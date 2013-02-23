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

#include "meshSurfaceEdgeExtractorFUN.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

#include <omp.h>

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshSurfaceEngine& meshSurfaceEdgeExtractorFUN::surfaceEngine()
{
    if( omp_in_parallel() )
        FatalErrorIn
        (
            "meshSurfaceEngine& meshSurfaceEdgeExtractorFUN::surfaceEngine()"
        ) << "Cannot create surface engine with a parallel region"
            << exit(FatalError);

    if( !surfaceEnginePtr_ )
        surfaceEnginePtr_ = new meshSurfaceEngine(mesh_);

    return *surfaceEnginePtr_;
}

void meshSurfaceEdgeExtractorFUN::clearOut()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and octree
meshSurfaceEdgeExtractorFUN::meshSurfaceEdgeExtractorFUN
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    meshOctree_(octree),
    surfaceEnginePtr_(NULL)
{
    if( Pstream::parRun() )
        FatalErrorIn
        (
            "meshSurfaceEdgeExtractorFUN::meshSurfaceEdgeExtractorFUN"
            "(polyMeshGen&, const meshOctree&)"
        ) << "Cannot run in parallel!" << exit(FatalError);

    distributeBoundaryFaces();

    reviseCorners();

    reviseEdges();

    //remapBoundaryPoints();

    createBasicFundamentalSheets();

    smoothMeshSurface();

    improveQualityOfFundamentalSheets();

    remapBoundaryPoints();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEdgeExtractorFUN::~meshSurfaceEdgeExtractorFUN()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
