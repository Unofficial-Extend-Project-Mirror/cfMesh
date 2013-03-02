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
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "correctEdgesBetweenPatches.H"

//#define DEBUGEdges

# ifdef DEBUGEdges
#include "writeMeshEnsight.H"
#include "primitiveMesh.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
void dualUnfoldConcaveCells::replaceBoundary()
{
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    wordList patchNames(boundaries.size());
    forAll(boundaries, patchI)
        patchNames[patchI] = boundaries[patchI].patchName();
    
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.replaceBoundary
    (
        patchNames,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
    
    newBoundaryFaces_.setSize(0);
    newBoundaryOwners_.setSize(0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMeshGen
dualUnfoldConcaveCells::dualUnfoldConcaveCells
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octree_(octree),
    typeOfCell_(mesh.cells().size(), INTERNALCELL),
    typeOfVertex_(mesh.points().size(), NONE),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_()
{
    mesh_.clearAddressingData();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dualUnfoldConcaveCells::~dualUnfoldConcaveCells()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualUnfoldConcaveCells::unfoldInvalidCells()
{
    const meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);
    
    if( findConcaveEdges(*msePtr) )
    {
        markVertexTypes(*msePtr);
        
        storeAndMergeBoundaryFaces(*msePtr);
        
        createNeighbouringBoundaryFaces(*msePtr);
        
        storeRemainingBoundaryFaces(*msePtr);
        
        removeConcaveVerticesFromIntFaces(*msePtr);
        
        deleteDemandDrivenData(msePtr);
        
        replaceBoundary();
        
        checkAndRepairBoundary();
        
        polyMeshGenModifier(mesh_).removeUnusedVertices();
        
        # ifdef DEBUGEdges
        mesh_.addressingData().checkMesh(true);
        mesh_.write();
        ::exit(EXIT_FAILURE);
        # endif
        
        correctEdgesBetweenPatches correctEdges(mesh_);
        meshSurfaceMapper(*msePtr, octree_).mapCornersAndEdges();
        
        # ifdef DEBUGEdges
        mesh_.addressingData().checkMesh();
        # endif
    }
    else
    {
        deleteDemandDrivenData(msePtr);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
