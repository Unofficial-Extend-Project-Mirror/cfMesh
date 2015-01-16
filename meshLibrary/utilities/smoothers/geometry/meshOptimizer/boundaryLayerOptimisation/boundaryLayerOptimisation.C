/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryLayerOptimisation::boundaryLayerOptimisation(polyMeshGen& mesh)
:
    mesh_(mesh),
    meshSurfacePtr_(new meshSurfaceEngine(mesh)),
    deleteMeshSurface_(true),
    partitionerPtr_(NULL),
    hairEdges_(),
    hairEdgesAtBndPoint_(),
    hairEdgesNearHairEdge_(),
    isBndLayerBase_(),
    isExitFace_(),
    hairEdgeType_(),
    thinnedHairEdge_()
{
    calculateHairEdges();
}

boundaryLayerOptimisation::boundaryLayerOptimisation
(
    polyMeshGen& mesh,
    const meshSurfaceEngine& mse
)
:
    mesh_(mesh),
    meshSurfacePtr_(&mse),
    deleteMeshSurface_(false),
    partitionerPtr_(NULL),
    hairEdges_(),
    hairEdgesAtBndPoint_(),
    hairEdgesNearHairEdge_(),
    isBndLayerBase_(),
    isExitFace_(),
    hairEdgeType_(),
    thinnedHairEdge_()
{
    calculateHairEdges();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryLayerOptimisation::~boundaryLayerOptimisation()
{
    deleteDemandDrivenData(partitionerPtr_);

    if( deleteMeshSurface_ )
        deleteDemandDrivenData(meshSurfacePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const edgeLongList& boundaryLayerOptimisation::hairEdges() const
{
    return hairEdges_;
}

const VRWGraph& boundaryLayerOptimisation::hairEdgesAtBndPoint() const
{
    return hairEdgesAtBndPoint_;
}

const boolList& boundaryLayerOptimisation::isBaseFace() const
{
    return isBndLayerBase_;
}

const boolList& boundaryLayerOptimisation::isExitFace() const
{
    return isExitFace_;
}

void boundaryLayerOptimisation::optimiseLayer
(
    const dictionary& meshDict,
    boundaryLayerOptimisation& blOptimisation
)
{
    label nIterations(20);
    scalar thicknessTol(0.15);
    scalar featureSizeTol(0.3);

    if( meshDict.found("boundaryLayers") )
    {
        const dictionary& blDict = meshDict.subDict("boundaryLayers");

        if( blDict.found("boundaryLayerSmoother") )
        {
            const dictionary& smoothDict =
                blDict.subDict("boundaryLayerSmoother");

            if( smoothDict.found("nIterations") )
                nIterations = readLabel(smoothDict.lookup("nIterations"));

            if( smoothDict.found("thicknessAngle") )
            {
                const scalar angle =
                    readScalar(smoothDict.lookup("thicknessAngle"));

                thicknessTol = Foam::tan(angle * M_PI / 180.0);
            }

            if( smoothDict.found("featureSizeTolerance") )
            {
                const scalar featureSize =
                    readScalar(smoothDict.lookup("featureSizeTolerance"));

                if( featureSizeTol >= 1.0 || featureSizeTol < 0.0 )
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::optimiseLayer"
                        "(const dictionary&, boundaryLayerOptimisation&)"
                    ) << "Invalid feature size tolerance is out"
                      << " of a valid range 0 to 1" << exit(FatalError);

                featureSizeTol = featureSize;
            }
        }
    }

    //- optimise hair normals
    blOptimisation.optimiseLayer(nIterations, thicknessTol, featureSizeTol);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
