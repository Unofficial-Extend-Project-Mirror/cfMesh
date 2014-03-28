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
#include "meshOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * Private member functions * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& meshOptimizer::meshSurface() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

void meshOptimizer::clearSurface()
{
    deleteDemandDrivenData(msePtr_);
}

label meshOptimizer::findBadFaces
(
    labelHashSet& badFaces,
    const boolList& changedFace
) const
{
    badFaces.clear();

    polyMeshGenChecks::checkFacePyramids
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceFlatness
    (
        mesh_,
        false,
        0.8,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkCellPartTetrahedra
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

label meshOptimizer::findLowQualityFaces
(
    labelHashSet& badFaces,
    const boolList& changedFace
) const
{
    badFaces.clear();

    polyMeshGenChecks::checkFaceDotProduct
    (
        mesh_,
        false,
        70.0,
        &badFaces
    );

    polyMeshGenChecks::checkFaceSkewness
    (
        mesh_,
        false,
        2.0,
        &badFaces
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
meshOptimizer::meshOptimizer(polyMeshGen& mesh)
:
    mesh_(mesh),
    vertexLocation_(mesh.points().size(), INSIDE),
    msePtr_(NULL)
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bPoints = mse.boundaryPoints();

    //- mark boundary vertices
    forAll(bPoints, bpI)
        vertexLocation_[bPoints[bpI]] = BOUNDARY;

    //- mark edge vertices
    meshSurfacePartitioner mPart(mse);
    forAllConstIter(labelHashSet, mPart.edgePoints(), it)
        vertexLocation_[bPoints[it.key()]] = EDGE;

    //- mark corner vertices
    forAllConstIter(labelHashSet, mPart.corners(), it)
        vertexLocation_[bPoints[it.key()]] = CORNER;

    if( Pstream::parRun() )
    {
        const polyMeshGenAddressing& addresing = mesh_.addressingData();
        const VRWGraph& pointAtProcs = addresing.pointAtProcs();

        forAll(pointAtProcs, pointI)
            if( pointAtProcs.sizeOfRow(pointI) != 0 )
                vertexLocation_[pointI] |= PARALLELBOUNDARY;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOptimizer::~meshOptimizer()
{
    clearSurface();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
