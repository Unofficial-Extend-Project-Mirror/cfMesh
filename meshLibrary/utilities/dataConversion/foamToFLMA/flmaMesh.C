/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2007 Franjo Juretic
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "flmaMesh.H"
#include "IOmanip.H"
#include "tetMatcher.H"
#include "hexMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::flmaMesh::createPointsForDecomposition()
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    
    forAll(faces, faceI)
    {
        if( faces[faceI].size() < 5 )
            continue;
        
        const point p = faces[faceI].centre(points);
        faceCentreLabel_.insert(faceI, points.size()+additionalPoints_.size());
        additionalPoints_.append(p);
    }
    
    const labelList& owner = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    cellType_.setSize(cells.size());
    cellType_ = -1;
    
    tetMatcher tet;
    pyrMatcher pyr;
    hexMatcher hex;
    prismMatcher prism;
    
    forAll(cells, cellI)
    {
        if( tet.matchShape(true, faces, owner, cellI, cells[cellI]) )
        {
            cellType_[cellI] = 4;
        }
        else if( hex.matchShape(true, faces, owner, cellI, cells[cellI]) )
        {
            cellType_[cellI] = 5;
        }
        else if( prism.matchShape(true, faces, owner, cellI, cells[cellI]) )
        {
            cellType_[cellI] = 8;
        }
        else if( pyr.matchShape(true, faces, owner, cellI, cells[cellI]) )
        {
            cellType_[cellI] = 6;
        }
        else
        {
            FatalError << "Cell " << cellI << " is not a standard cell!!"
                << exit(FatalError);
            
            const labelList cp = cells[cellI].labels(faces);
            point p(vector::zero);
            forAll(cp, i)
                p += points[cp[i]];
            
            p /= cp.size();
            
            cellCentreLabel_.insert
            (
                cellI,
                points.size()+additionalPoints_.size()
            );
            additionalPoints_.append(p);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMeshGen
Foam::flmaMesh::flmaMesh(const polyMeshGen& mesh)
:
    mesh_(mesh),
    cellType_(),
    additionalPoints_(),
    faceCentreLabel_(),
    cellCentreLabel_()
{
    createPointsForDecomposition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flmaMesh::~flmaMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flmaMesh::writePoints
(
    Foam::OFstream& flmaGeometryFile
) const
{
    flmaGeometryFile << (mesh_.points().size()+additionalPoints_.size()) << nl;
    const pointFieldPMG& points = mesh_.points();
    forAll(points, pointI)
    {
        const point& p = points[pointI];
        flmaGeometryFile << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';
    }
    forAll(additionalPoints_, pointI)
    {
        const point& p = additionalPoints_[pointI];
        flmaGeometryFile << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';
    }
    
    flmaGeometryFile << nl;
}

void flmaMesh::writeCells
(
    OFstream& flmaGeometryFile
) const
{
    const labelList& owner = mesh_.owner();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    tetMatcher tet;
    pyrMatcher pyr;
    hexMatcher hex;
    prismMatcher prism;
    
    flmaGeometryFile << cells.size() << nl;
    
    forAll(cells, cellI)
    {
        if( tet.matchShape(false, faces, owner, cellI, cells[cellI]) )
        {
            const labelList& tetVrt = tet.vertLabels();
            flmaGeometryFile << tetVrt.size();
            forAll(tetVrt, i)
                flmaGeometryFile << ' ' << tetVrt[i];
            flmaGeometryFile << nl;
        }
        else if( hex.matchShape(false, faces, owner, cellI, cells[cellI]) )
        {
            const labelList& hexVrt = hex.vertLabels();
            flmaGeometryFile << hexVrt.size();
            forAll(hexVrt, i)
                flmaGeometryFile << ' ' << hexVrt[i];
            flmaGeometryFile << nl;
        }
        else if( prism.matchShape(false, faces, owner, cellI, cells[cellI]) )
        {
            const labelList& prismVrt = prism.vertLabels();
            flmaGeometryFile << prismVrt.size();
            forAll(prismVrt, i)
                flmaGeometryFile << ' ' << prismVrt[i];
            flmaGeometryFile << nl;
        }
        else if( pyr.matchShape(false, faces, owner, cellI, cells[cellI]) )
        {
            const labelList& pyrVrt = pyr.vertLabels();
            flmaGeometryFile << pyrVrt.size();
            forAll(pyrVrt, i)
                flmaGeometryFile << ' ' << pyrVrt[i];
            flmaGeometryFile << nl;
        }
    }
}

void Foam::flmaMesh::writeCellTypes
(
    OFstream& flmaGeometryFile
) const
{
    flmaGeometryFile << nl << cellType_.size() << nl;
    forAll(cellType_, cellI)
        flmaGeometryFile << cellType_[cellI] << nl;
    flmaGeometryFile << nl;
}

void Foam::flmaMesh::writeSelections
(
    Foam::OFstream& flmaGeometryFile
) const
{
    //- write patches as face selections
    const PtrList<writePatch>& patches = mesh_.boundaries();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    tetMatcher tet;
    pyrMatcher pyr;
    hexMatcher hex;
    prismMatcher prism;
    
    label nSubsets(0);
    
    nSubsets += patches.size();
    
    DynList<label> indices;
    mesh_.pointSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " point subsets" << endl;
    mesh_.faceSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " face subsets" << endl;
    mesh_.cellSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " cell subsets" << endl;
    
    flmaGeometryFile << nSubsets << nl;
    
    //- write patches as face selections
    forAll(patches, patchI)
    {
        const writePatch& patch = patches[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();
        
        flmaGeometryFile << patch.patchName() << nl;
        flmaGeometryFile << 3 << nl;
        flmaGeometryFile << 2 * patch.patchSize() << nl;
        
        for(label faceI=start;faceI<end;++faceI)
        {
            const label cellI = owner[faceI];
            const cell& c = cells[owner[faceI]];
            
            if( tet.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = tet.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( hex.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = hex.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( prism.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = prism.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( pyr.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = pyr.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
        }
        
        flmaGeometryFile << nl;
    }
    
    //- write node selections
    mesh_.pointSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelListPMG nodesInSubset;
        mesh_.pointsInSubset(indices[indexI], nodesInSubset);
        
        flmaGeometryFile << mesh_.pointSubsetName(indices[indexI]) << nl;
        flmaGeometryFile << 1 << nl;
        flmaGeometryFile << nodesInSubset.size() << nl;
        forAll(nodesInSubset, i)
            flmaGeometryFile << nodesInSubset[i] << ' ';
        flmaGeometryFile << nl;
    }
    
    //- write face selections
    mesh_.faceSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelListPMG facesInSubset;
        mesh_.facesInSubset(indices[indexI], facesInSubset);
        
        flmaGeometryFile << mesh_.faceSubsetName(indices[indexI]) << nl;
        flmaGeometryFile << 3 << nl;
        flmaGeometryFile << 2 * facesInSubset.size() << nl;
        forAll(facesInSubset, i)
        {
            const label faceI = facesInSubset[i];
            const label cellI = owner[faceI];
            const cell& c = cells[owner[faceI]];
            
            if( tet.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = tet.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( hex.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = hex.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( prism.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = prism.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
            else if( pyr.matchShape(false, faces, owner, cellI, c) )
            {
                const labelList& faceMap = pyr.faceMap();
                label dir(-1);
                forAll(faces, i)
                {
                    if( faceI == faceMap[i] )
                    {
                        dir = i;
                        break;
                    }
                }
                
                flmaGeometryFile << ' ' << cellI << ' ' << dir;
            }
        }
        
        flmaGeometryFile << nl;
    }
    
    //- write cell selections
    mesh_.cellSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelListPMG cellsInSubset;
        mesh_.cellsInSubset(indices[indexI], cellsInSubset);
        
        flmaGeometryFile << mesh_.cellSubsetName(indices[indexI]) << nl;
        flmaGeometryFile << 2 << nl;
        flmaGeometryFile << cellsInSubset.size() << nl;
        forAll(cellsInSubset, i)
            flmaGeometryFile << cellsInSubset[i] << ' ';
        flmaGeometryFile << nl;
    }
}


void flmaMesh::write
(
    OFstream& flmaGeometryFile
) const
{
    writePoints(flmaGeometryFile);

    writeCells(flmaGeometryFile);
    
    writeCellTypes(flmaGeometryFile);

    writeSelections(flmaGeometryFile);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
