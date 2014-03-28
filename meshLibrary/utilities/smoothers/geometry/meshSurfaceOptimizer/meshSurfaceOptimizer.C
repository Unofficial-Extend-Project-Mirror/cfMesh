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
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGenChecks.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::classifySurfaceVertices()
{
    const labelHashSet& corners = partitioner_.corners();
    const labelHashSet& edgePoints = partitioner_.edgePoints();

    //- set all vertices to partition
    vertexType_ = PARTITION;

    //- set corners
    forAllConstIter(labelHashSet, corners, it)
        vertexType_[it.key()] = CORNER;

    //- set edges
    forAllConstIter(labelHashSet, edgePoints, it)
        vertexType_[it.key()] = EDGE;

    if( Pstream::parRun() )
    {
        //- mark nodes at parallel boundaries
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            vertexType_[bpI] |= PROCBND;
        }
    }
}

label meshSurfaceOptimizer::findBadFaces
(
    labelHashSet& badFaces,
    boolList& changedFace
) const
{
    badFaces.clear();

    const polyMeshGen& mesh = surfaceEngine_.mesh();

    polyMeshGenChecks::checkFacePyramids
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkCellPartTetrahedra
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh surface and octree
meshSurfaceOptimizer::meshSurfaceOptimizer
(
    meshSurfaceEngine& surface,
    const meshOctree& octree
)
:
    surfaceEngine_(surface),
    meshOctree_(octree),
    vertexType_(surface.boundaryPoints().size()),
    trianglesPtr_(NULL),
    pointTrianglesPtr_(NULL),
    partitioner_(surface)
{
    classifySurfaceVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceOptimizer::~meshSurfaceOptimizer()
{
    deleteDemandDrivenData(trianglesPtr_);
    deleteDemandDrivenData(pointTrianglesPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
