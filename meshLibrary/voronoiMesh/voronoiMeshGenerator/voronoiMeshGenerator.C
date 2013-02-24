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

#include "voronoiMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "objectRegistry.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "voronoiMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "triSurfacePatchManipulator.H"

#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "meshSurfaceEdgeExtractorFUN.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void voronoiMeshGenerator::createVoronoiMesh()
{
    //- create voronoi mesh from octree and Delaunay tets
    voronoiMeshExtractor vme(*octreePtr_, meshDict_, mesh_);

    vme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    do
    {
        surfaceMorpherCells* cmPtr = new surfaceMorpherCells(mesh_);
        cmPtr->morphMesh();
        deleteDemandDrivenData(cmPtr);
    }
    while( topologicalCleaner(mesh_).cleanTopology() );

/*    bool changed;
    do
    {
        changed = false;

        checkIrregularSurfaceConnections checkConnections(mesh_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;

        if( checkNonMappableCellConnections(mesh_).removeCells() )
            changed = true;

        if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
            changed = true;
    } while( changed );

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();
    */

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- map mesh surface on the geometry surface
    meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();
    mapper.mapVerticesOntoSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif

    //- untangle surface faces
    meshSurfaceOptimizer(*msePtr, *octreePtr_).untangleSurface();

    # ifdef DEBUG
    mesh_.write();
    ::exit(EXIT_FAILURE);
    # endif

    deleteDemandDrivenData(msePtr);
}

void voronoiMeshGenerator::mapEdgesAndCorners()
{
    //meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);
    meshSurfaceEdgeExtractorFUN(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer surfOptimiser(mse, *octreePtr_);
    surfOptimiser.optimizeSurface();
    surfOptimiser.untangleSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::decomposeConcaveCells()
{
    //decomposeCellsNearConcaveEdges dm(mesh_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::generateBoudaryLayers()
{
    boundaryLayers bl(mesh_);

    if( meshDict_.found("boundaryLayers") )
    {
        boundaryLayers bl(mesh_);

        const dictionary& bndLayers = meshDict_.subDict("boundaryLayers");

        if( bndLayers.found("nLayers") )
        {
            const label nLayers = readLabel(bndLayers.lookup("nLayers"));

            if( nLayers > 0 )
                bl.addLayerForAllPatches();
        }
        else if( bndLayers.found("patchBoundaryLayers") )
        {
            const dictionary& patchLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList createLayers = patchLayers.toc();

            forAll(createLayers, patchI)
                bl.addLayerForPatch(createLayers[patchI]);
        }
    }

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void voronoiMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.refineLayers();

        meshOptimizer optimizer(mesh_);

        optimizer.untangleMeshFV();
    }
}

void voronoiMeshGenerator::optimiseFinalMesh()
{
    //- final optimisation
    meshOptimizer optimizer(mesh_);

    optimizer.optimizeSurface(*octreePtr_);

    deleteDemandDrivenData(octreePtr_);

    optimizer.optimizeMeshFV();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void voronoiMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void voronoiMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void voronoiMeshGenerator::generateMesh()
{
    createVoronoiMesh();

    surfacePreparation();

    mapMeshToSurface();

    mapEdgesAndCorners();

    optimiseMeshSurface();

    generateBoudaryLayers();

    optimiseFinalMesh();

    refBoundaryLayers();

    renumberMesh();

    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Time
voronoiMeshGenerator::voronoiMeshGenerator(const Time& time)
:
    runTime_(time),
    surfacePtr_(NULL),
    octreePtr_(NULL),
    pointRegionsPtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            runTime_.system(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(time)
{
    if( true )
        checkMeshDict cmd(meshDict_);

    const fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurf(runTime_.path()/surfaceFile);

    if( surfacePtr_->featureEdges().size() != 0 )
    {
        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr_);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches(&meshDict_);

        //- delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr_);
        surfacePtr_ = surfaceWithPatches;
    }

    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator* octreeCreatorPtr =
        new meshOctreeCreator(*octreePtr_, meshDict_);
    octreeCreatorPtr->createOctreeBoxes();
    deleteDemandDrivenData(octreeCreatorPtr);

    generateMesh();
}

/*
voronoiMeshGenerator::voronoiMeshGenerator
(
    const Time& time,
    const volScalarField& localCellSize
)
:
    runTime_(time),
    surfacePtr_(NULL),
    octreePtr_(NULL),
    pointRegionsPtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(time)
{
    fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurface(runTime_.path()/surfaceFile);

    octreePtr_ = new meshOctree(*surfacePtr_);

    generateMesh();
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voronoiMeshGenerator::~voronoiMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
    deleteDemandDrivenData(pointRegionsPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
