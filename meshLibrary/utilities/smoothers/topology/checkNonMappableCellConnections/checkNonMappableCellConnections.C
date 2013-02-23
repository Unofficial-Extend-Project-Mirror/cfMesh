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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "checkNonMappableCellConnections.H"
#include "polyMeshGenModifier.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"

#include <omp.h>

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void checkNonMappableCellConnections::findCellTypes()
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    
    cellType_.setSize(cells.size());
    cellType_ = INTERNALCELL;
    
    //- find boundary cells
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();
        
        for(label faceI=start;faceI<end;++faceI)
            cellType_[owner[faceI]] = BNDCELL;
    }
    
    //- find boundary cells with all vertices at the boundary
    meshSurfaceEngine mse(mesh_);
    const labelList& bp = mse.bp();
    
    # pragma omp parallel for schedule(dynamic, 1000)
    for(label cellI=cells.size()-1;cellI>=0;--cellI)
    {
        if( cellType_[cellI] & INTERNALCELL )
            continue;
        
        const cell& c = cells[cellI];
        
        //- mark boundary cells with all vertices at the boundary
        const labelList cellPoints = c.labels(faces);
        bool allBoundary(true);
        forAll(cellPoints, cpI)
        {
            if( bp[cellPoints[cpI]] < 0 )
            {
                allBoundary = false;
                break;
            }
        }
        
        if( allBoundary )
        {
            cellType_[cellI] |= ALLBNDVERTEXCELL;
        }
        else
        {
            continue;
        }
        
        //- check if the internal faces are connected into a single group
        //- over their edges
        DynList<label> internalFaces;
        forAll(c, fI)
        {
            if( c[fI] < mesh_.nInternalFaces() )
            {
                internalFaces.append(c[fI]);
            }
            else if( mesh_.faceIsInProcPatch(c[fI]) != -1 )
            {
                internalFaces.append(c[fI]);
            }
        }
        
        Map<label> faceGroup(internalFaces.size());
        label nGroup(0);
        forAll(internalFaces, i)
        {
            if( faceGroup.found(internalFaces[i]) )
                continue;
            
            DynList<label> front;
            front.append(internalFaces[i]);
            faceGroup.insert(internalFaces[i], nGroup);
            
            while( front.size() )
            {
                const label fLabel = front.removeLastElement();
                
                forAll(internalFaces, j)
                {
                    const label nei = internalFaces[j];
                    if( faceGroup.found(nei) )
                        continue;
                    
                    if( help::shareAnEdge(faces[fLabel], faces[nei]) )
                    {
                        front.append(nei);
                        faceGroup.insert(nei, nGroup);
                    }
                }
            }
            
            ++nGroup;
        }

        if( nGroup > 1 )
            cellType_[cellI] |= INTERNALFACEGROUP;
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkNonMappableCellConnections::checkNonMappableCellConnections
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    cellType_()
{
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkNonMappableCellConnections::~checkNonMappableCellConnections()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkNonMappableCellConnections::findCells(labelHashSet& badCells)
{
    badCells.clear();
    
    //- classify cell types
    findCellTypes();
    
    //- select ALLBNDVERTEXCELL and INTERNALFACEGROUP cells
    //- with at least one INTERNALCELL neighbour
    //- these cells do not need to stay in the mesh
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    const PtrList<writeProcessorPatch>& procBoundaries = mesh_.procBoundaries();
    
    labelListList otherProcType;
    if( Pstream::parRun() )
    {
        //- exchange cell types at processor boundaries
        otherProcType.setSize(procBoundaries.size());
        
        //- send data to other processors
        forAll(procBoundaries, patchI)
        {
            label start = procBoundaries[patchI].patchStart();
            labelList patchCellType(procBoundaries[patchI].patchSize());
            
            forAll(patchCellType, faceI)
                patchCellType[faceI] = cellType_[owner[start++]];
            
            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                patchCellType.byteSize()
            );
            
            toOtherProc << patchCellType;
        }
        
        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            
            labelList& otherTypes = otherProcType[patchI];
            fromOtherProc >> otherTypes;
        }
    }
    
    # pragma omp parallel for schedule(dynamic)
    for(label cellI=cellType_.size()-1;cellI>=0;--cellI)
    {
        if( cellType_[cellI] & INTERNALFACEGROUP )
        {
            # pragma omp critical
            badCells.insert(cellI);
        }
        else if( cellType_[cellI] & (ALLBNDVERTEXCELL+INTERNALFACEGROUP) )
        {
            //- mark cells which have only one internal neighbour
            const cell& c = cells[cellI];
            
            bool hasInternalNeighbour(false);
            label nNeiCells(0);
            
            forAll(c, fI)
            {
                const label faceI = c[fI];
                
                if( faceI < nInternalFaces )
                {
                    ++nNeiCells;
                    
                    label nei = neighbour[c[fI]];
                    if( nei == cellI )
                        nei = owner[c[fI]];
                    
                    if( cellType_[nei] & INTERNALCELL )
                    {
                        hasInternalNeighbour = true;
                        break;
                    }
                }
                else if( mesh_.faceIsInProcPatch(faceI) != -1 )
                {
                    ++nNeiCells;
                    
                    const label patchI = mesh_.faceIsInProcPatch(faceI);
                    const label j = faceI - procBoundaries[patchI].patchStart();
                    
                    if( otherProcType[patchI][j] & INTERNALCELL )
                    {
                        hasInternalNeighbour = true;
                        break;
                    }
                }
            }
            
            if( hasInternalNeighbour || (nNeiCells == 1) )
            {
                # pragma omp critical
                badCells.insert(cellI);
            }
        }
    }
}

bool checkNonMappableCellConnections::removeCells()
{
    labelHashSet badCells;
    
    label nRemoved;
    bool changed(false);
    
    do
    {
        findCells(badCells);
        
        nRemoved = badCells.size();
        reduce(nRemoved, sumOp<label>());
        
        Info << "Found " << nRemoved << " non-mappable cells" << endl;
        
        if( nRemoved != 0 )
        {
            boolList removeCell(mesh_.cells().size(), false);
            forAllConstIter(labelHashSet, badCells, it)
                removeCell[it.key()] = true;
            
            polyMeshGenModifier(mesh_).removeCells(removeCell);
            
            changed = true;
        }
    } while( nRemoved );
    
    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
