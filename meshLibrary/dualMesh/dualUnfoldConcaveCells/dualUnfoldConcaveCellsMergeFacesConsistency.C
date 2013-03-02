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

#include "demandDrivenData.H"
#include "dualUnfoldConcaveCells.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"
#include "polyMeshGenAddressing.H"
#include "tetrahedron.H"

//#define DEBUGEdges

# ifdef DEBUGEdges
#include "cellSet.H"
#include "objectRegistry.H"
#include "Time.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualUnfoldConcaveCells::createNeighbouringBoundaryFaces
(
    const meshSurfaceEngine& mse
)
{
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelList& faceOwner = mse.faceOwners();
    const labelList& bp = mse.bp();
    
    //- find the patch for merged faces which belong to the given cell
    //- this here assumes that the original cell before edge extraction
    //- has only one boundary face!! This is the case if surface preparation
    //- is applied before edge extraction
    labelListPMG newPatchForCell(mesh_.cells().size(), -1);
    forAll(newBoundaryFaces_, faceI)
    {
        if( newPatchForCell[newBoundaryOwners_[faceI]] == -1 )
        {
            newPatchForCell[newBoundaryOwners_[faceI]] =
                newBoundaryPatches_[faceI];
        }
        else
        {
            FatalErrorIn
            (
                "void dualUnfoldConcaveCells::"
                "createNeighbouringBoundaryFaces(const meshSurfaceEngine&)"
            ) << "Cell " << newBoundaryOwners_[faceI]
                << " has more than one "
                << " boundary face! Cannot proceed!" << abort(FatalError);
        }
    }

    //- create a list of possible candidates to store
    labelListPMG front;
    forAll(bPoints, bpI)
        if( typeOfVertex_[bPoints[bpI]] & REMOVE )
        {
            front.append(bpI);
        }
        
    # ifdef DEBUGEdges
    Info << "Front points " << front << endl;
    # endif
    
    //- start creating new boundary faces
    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();
        
    while( front.size() != 0 )
    {
        labelListPMG newFront;
        
        forAll(front, fpI)
        {
            const label fLabel = front[fpI];
            
            label cellI(-1), patchOfTreated(-1);
            forAllRow(pointFaces, fLabel, pfI)
            {
                const label pfLabel = pointFaces(fLabel, pfI);
                
                if(
                    typeOfCell_[faceOwner[pfLabel]] ==
                    BOUNDARYCELL
                )
                {
                    cellI = faceOwner[pfLabel];

                }
                else if( newPatchForCell[faceOwner[pfLabel]] != -1 )
                {
                    patchOfTreated = newPatchForCell[faceOwner[pfLabel]];
                }
            }
                
            if( cellI == -1 )
                continue;
            
            typeOfCell_[cellI] |= TREATEDCELL;
            
            # ifdef DEBUGEdges
            Info << "Checking cell " << cellI << endl;
            Info << "Patch of treated " << patchOfTreated << endl;
            Info << "Front point " << bPoints[front[fpI]] << endl;
            # endif
            
            //- find boundary faces for the given cell
            DynList<face> bFaces(5);
            DynList<label> facePatch(5);
            forAll(cells[cellI], fI)
            {
                const label patch = mesh_.faceIsInPatch(cells[cellI][fI]);
                if( patch != -1 )
                {
                    bFaces.append(faces[cells[cellI][fI]]);
                    facePatch.append(patch);
                }
            }
            
            # ifdef DEBUGEdges
            Info << "Boundary faces " << bFaces << endl;
            Info << "Face patches " << facePatch << endl;
            # endif
            
            //- check if faces have to be merged
            bool merge(false);
            forAll(bFaces, bfI)
            {
                const face& f = bFaces[bfI];
                direction nPointsInFace(0);
                
                # ifdef DEBUGEdges
                DynList<label> removePoint(2);
                # endif
                
                forAll(f, pI)
                    if( !(typeOfVertex_[f[pI]] & REMOVE) )
                    {
                        ++nPointsInFace;
                    }
                    # ifdef DEBUGEdges
                    else
                    {
                        removePoint.append(f[pI]);
                    }
                    # endif
                    
                # ifdef DEBUGEdges
                Info << "Vertices removed from boundary face " << bfI
                    << " are " << removePoint << endl;
                # endif
                    
                if( nPointsInFace < 3 )
                {
                    merge = true;
                    break;
                }
            }
            
            if( (patchOfTreated == -1) || !facePatch.contains(patchOfTreated) )
                merge = true;
            
            if( merge )
            {
                # ifdef DEBUGEdges
                Info << "Merging boundary faces of cell" << endl;
                # endif
                //- merge boundary faces of the given cell
                newPatchForCell[cellI] = mergeBoundaryFacesOfCell(cellI);
                
                forAll(bFaces, bfI)
                {
                    const face& f = bFaces[bfI];
                    forAll(f, pI)
                        if( typeOfVertex_[f[pI]] & REMOVE )
                            newFront.append(bp[f[pI]]);
                }
            }
            else
            {
                # ifdef DEBUGEdges
                Info << "Changing boundary faces" << endl;
                # endif
                
                //- change boundary faces
                face triF(3);
                forAll(bFaces, bfI)
                    if( facePatch[bfI] != patchOfTreated )
                    {
                        const face& bf = bFaces[bfI];
                        //- remove edge vertex from face not it the same patch
                        //- as already treated cell
                        face newBf(bf.size() - 1);
                        label rpos(-1), pJ(0);
                        forAll(bf, pI)
                            if( typeOfVertex_[bf[pI]] & REMOVE )
                            {
                                rpos = pI;
                            }
                            else
                            {
                                newBf.newElmt(pJ++) = bf[pI];
                            }
                            
                        //- store shrinked boundary face
                        newBf.setSize(pJ);
                        newBoundaryFaces_.appendList(newBf);
                        newBoundaryOwners_.append(cellI);
                        newBoundaryPatches_.append(facePatch[bfI]);
                            
                        if( rpos != -1 )
                        {
                            triF[0] = bf.prevLabel(rpos);
                            triF[1] = bf[rpos];
                            triF[2] = bf.nextLabel(rpos);
                            
                            # ifdef DEBUGEdges
                            Info << "triF " << triF << endl;
                            # endif
                        }
                    }
                    
                forAll(bFaces, bfI)
                    if( facePatch[bfI] == patchOfTreated )
                    {
                        //- merge face with triF and store it
                        const face mf = help::mergeTwoFaces(bFaces[bfI], triF);
                        
                        # ifdef DEBUGEdges
                        Info << "Merged face " << mf << endl;
                        # endif
                        
                        face shrinkedFace(mf.size() -1);
                        direction i(0);
                        forAll(mf, pI)
                            if( !(typeOfVertex_[mf[pI]] & REMOVE) )
                                shrinkedFace[i++] = mf[pI];
                        
                        # ifdef DEBUGEdges
                        Info << "Storing face " << shrinkedFace << " into patch"
                            << patchOfTreated << endl;
                        # endif
                        
                        newBoundaryFaces_.appendList(shrinkedFace);
                        newBoundaryOwners_.append(cellI);
                        newBoundaryPatches_.append(patchOfTreated);
                    }
            }
        }
        
        front = newFront;
    }
}

void dualUnfoldConcaveCells::storeRemainingBoundaryFaces
(
    const meshSurfaceEngine& mse
)
{
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwner = mse.faceOwners();
    const labelList& facePatch = mse.boundaryFacePatches();
    
    forAll(bFaces, bfI)
        if( !(typeOfCell_[faceOwner[bfI]] & TREATEDCELL) )
        {
            newBoundaryFaces_.appendList(bFaces[bfI]);
            newBoundaryOwners_.append(faceOwner[bfI]);
            newBoundaryPatches_.append(facePatch[bfI]);
        }
}

void dualUnfoldConcaveCells::removeConcaveVerticesFromIntFaces
(
    const meshSurfaceEngine& mse
)
{
    const label nIntFaces = mesh_.nInternalFaces();
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    const VRWGraph& pointFaces = mesh_.addressingData().pointFaces();
    forAll(pointFaces, pI)
        if( (typeOfVertex_[pI] & REMOVE) )
        {
            # ifdef DEBUGEdges
            Info << "Removing vertex " << pI << " from internal faces" << endl;
            # endif
            
            forAllRow(pointFaces, pI, pfI)
            {
                const label pointFaceI = pointFaces(pI, pfI);
                if( pointFaceI < nIntFaces )
                {
                    const face& f = faces[pointFaceI];
                    
                    DynList<label> newF(f.size()-1);
                    forAll(f, pJ)
                        if( f[pJ] != pI )
                            newF.append(f[pJ]);
                        
                    if( newF.size() > 2 )
                    {
                        # ifdef DEBUGEdges
                        Info << "Internal face " << f << " is shrinked to"
                            << newF << endl;
                        # endif
                        newF.shrink();
                        faces[pointFaceI] = face(newF);
                    }
                }
            }
        }
        
    mesh_.clearAddressingData();
}

void dualUnfoldConcaveCells::checkAndRepairBoundary()
{
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    
    newBoundaryFaces_.setSize(boundaries.size());
    newBoundaryOwners_.setSize(boundaries.size());
    
    boolList touchedOwner(mesh_.cells().size(), false);
    
    bool changed(false);
    
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();
        
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];
            const label own = owner[faceI];
            
            bool merge(false);
            
            forAll(f, pI)
                if( typeOfVertex_[f[pI]] & REMOVE )
                {
                    merge = true;
                    break;
                }
                
            if( merge && (typeOfCell_[own] & TREATEDCELL) )
            {
                changed = true;
                if( touchedOwner[own] ) continue;
                    
                touchedOwner[own] = true;
                mergeBoundaryFacesOfCell(own);
            }
            else
            {
                newBoundaryFaces_.appendList(f);
                newBoundaryOwners_.append(own);
                newBoundaryPatches_.append(patchI);
            }
        }
    }
    
    if( changed )
        replaceBoundary();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
