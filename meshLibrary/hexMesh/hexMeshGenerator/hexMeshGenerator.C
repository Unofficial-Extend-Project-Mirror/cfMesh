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

#include "hexMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "hexMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorFUN.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "triSurfacePatchManipulator.H"

//#define DEBUG
//#define DEBUGflma

# ifdef DEBUG
#include "writeMeshEnsight.H"
#include "writeMeshFLMA.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void hexMeshGenerator::generateOctree()
{
    if( !octreePtr_ )
        octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator creator(*octreePtr_, meshDict_);
    creator.activateHexRefinement();
    creator.createOctreeBoxes();
}

void hexMeshGenerator::generateDualMesh()
{
    //- create polyMesh from octree boxes
    hexMeshExtractor dme(*octreePtr_, meshDict_, mesh_);
    dme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "hexMesh");
    # else
    writeMeshEnsight(mesh_, "hexMesh");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void hexMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells which create non-manifold surfaces
    //- and the cells which connect two clusters of cells
    //- and make sure all cells make a single cluster of cells
    bool changed;
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

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "afterTopoCleaning");
    # else
    writeMeshEnsight(mesh_, "afterTopoCleaning");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void hexMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- map mesh surface on the geometry surface
    meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();
    mapper.mapVerticesOntoSurface();

    # ifdef DEBUG
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "afterMapping");
    # else
    writeMeshEnsight(mesh_, "afterMapping");
    # endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif

    //- untangle surface faces
    meshSurfaceOptimizer(*msePtr, *octreePtr_).preOptimizeSurface();

    # ifdef DEBUG
    //meshOptimizer(*octreePtr_, mesh_).preOptimize();
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "afterSurfaceSmoothing");
    # else
    writeMeshEnsight(mesh_, "afterSurfaceSmoothing");
    # endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif

    deleteDemandDrivenData(msePtr);

    //- extract edges and corners
    meshSurfaceEdgeExtractorFUN(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "withEdges");
    # else
    writeMeshEnsight(mesh_, "withEdges");
    #endif
    //::exit(EXIT_FAILURE);
    # endif
}

void hexMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "optSurfaceWithEdges");
    # else
    writeMeshEnsight(mesh_, "optSurfaceWithEdges");
    #endif
    # endif
}

void hexMeshGenerator::generateBoundaryLayers()
{
    boundaryLayers bl(mesh_);

    if( meshDict_.found("boundaryLayers") )
    {
        wordList createLayers;

        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& dict = meshDict_.subDict("boundaryLayers");
            createLayers = dict.toc();
        }
        else
        {
            wordList bndLayers(meshDict_.lookup("boundaryLayers"));
            createLayers.transfer(bndLayers);
        }

        forAll(createLayers, patchI)
            bl.addLayerForPatch(createLayers[patchI]);
    }
    else
    {
        //bl.createOTopologyLayers();
        bl.addLayerForAllPatches();
    }

    # ifdef DEBUG
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_, "meshWithBndLayer");
    # else
    writeMeshEnsight(mesh_, "meshWithBndLayer");
    # endif
    mesh_.write();
    # endif
}

void hexMeshGenerator::optimiseFinalMesh()
{
    //- final optimisation
    meshOptimizer optimizer(mesh_);

    optimizer.optimizeSurface(*octreePtr_);

    deleteDemandDrivenData(octreePtr_);

    optimizer.optimizeMeshFV();

    # ifdef DEBUG
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_,"optimisedMesh");
    # else
    writeMeshEnsight(mesh_, "optimisedMesh");
    #endif
    # endif
}

void hexMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_,"renamedPatchesMesh");
    # else
    writeMeshEnsight(mesh_, "renamedPatchesMesh");
    #endif
    # endif
}

void hexMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    # ifdef DEBUGflma
    writeMeshFLMA(mesh_,"renumberedMesh");
    # else
    writeMeshEnsight(mesh_, "renumberedMesh");
    #endif
    # endif
}

void hexMeshGenerator::generateMesh()
{
    generateDualMesh();

    surfacePreparation();

    mapMeshToSurface();

    optimiseMeshSurface();

    generateBoundaryLayers();

    optimiseFinalMesh();

    renumberMesh();

    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
hexMeshGenerator::hexMeshGenerator
(
    const Time& runTime
)
:
    runTime_(runTime),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(runTime)
{
    if( Pstream::parRun() )
    {
        FatalError << "Cannot run in parallel" << exit(FatalError);
    }

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

    generateOctree();

    generateMesh();
}

/*
hexMeshGenerator::hexMeshGenerator
(
    const Time& runTime,
    const volScalarField& localCellSize
)
:
    runTime_(runTime),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(time)
{
    fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurface(db_.path()/surfaceFile);

    generateOctree();

    generateMesh();
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hexMeshGenerator::~hexMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void hexMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
