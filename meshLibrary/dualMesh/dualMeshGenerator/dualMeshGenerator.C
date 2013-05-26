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

#include "dualMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "dualMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "dualUnfoldConcaveCells.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"

//#define DEBUG
//#define DEBUGfpma

# ifdef DEBUG
#include "writeMeshEnsight.H"
#include "writeMeshFPMA.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void dualMeshGenerator::generateDualMesh()
{
    //- create polyMesh from octree boxes
    dualMeshExtractor dme(*octreePtr_, meshDict_, mesh_);
    dme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "dualMesh");
    # else
    writeMeshEnsight(mesh_, "dualMesh");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    do
    {
        surfaceMorpherCells* cmPtr = new surfaceMorpherCells(mesh_);
        cmPtr->morphMesh();
        deleteDemandDrivenData(cmPtr);
    } while( topologicalCleaner(mesh_).cleanTopology() );

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterTopoCleaning");
    # else
    writeMeshEnsight(mesh_, "afterTopoCleaning");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- map mesh surface on the geometry surface
    meshSurfaceMapper(*msePtr, *octreePtr_).mapVerticesOntoSurface();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterMapping");
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
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterSurfaceSmoothing");
    # else
    writeMeshEnsight(mesh_, "afterSurfaceSmoothing");
    # endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif

    deleteDemandDrivenData(msePtr);

    //- extract edges and corners
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "withEdges");
    # else
    writeMeshEnsight(mesh_, "withEdges");
    #endif
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "optSurfaceWithEdges");
    # else
    writeMeshEnsight(mesh_, "optSurfaceWithEdges");
    #endif
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::checkConcaveEdges()
{
    //- optimize surface to get rid of nearby vertices
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);
    meshSurfaceOptimizer(*msePtr, *octreePtr_).optimizeSurface();
    deleteDemandDrivenData(msePtr);

    //- repair mesh near concave edges
    dualUnfoldConcaveCells(mesh_, *octreePtr_).unfoldInvalidCells();

    # ifdef DEBUG
    mesh_.write();
    //meshOptimizer(*octreePtr_, mesh_).preOptimize();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "correctedEdges");
    # else
    writeMeshEnsight(mesh_, "correctedEdges");
    #endif
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::generateBoudaryLayers()
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
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "meshWithBndLayer");
    # else
    writeMeshEnsight(mesh_, "meshWithBndLayer");
    # endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif
}

void dualMeshGenerator::optimiseFinalMesh()
{
    //- final optimisation
    meshOptimizer optimizer(mesh_);

    optimizer.optimizeSurface(*octreePtr_);

    deleteDemandDrivenData(octreePtr_);

    optimizer.optimizeMeshFV();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_,"optimisedMesh");
    # else
    writeMeshEnsight(mesh_, "optimisedMesh");
    #endif
    # endif
}

void dualMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_,"renamedPatchesMesh");
    # else
    writeMeshEnsight(mesh_, "renamedPatchesMesh");
    #endif
    # endif
}

void dualMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_,"renumberedMesh");
    # else
    writeMeshEnsight(mesh_, "renumberedMesh");
    #endif
    # endif
}

void dualMeshGenerator::generateMesh()
{
    generateDualMesh();

    surfacePreparation();

    mapMeshToSurface();

    optimiseMeshSurface();

    checkConcaveEdges();

    generateBoudaryLayers();

    optimiseFinalMesh();

    renumberMesh();

    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
dualMeshGenerator::dualMeshGenerator
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
    if( true )
        checkMeshDict cmd(meshDict_);

    const fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurf(runTime_.path()/surfaceFile);

//     if( meshDict_.found("subsetFileName") )
//     {
//         const fileName subsetFileName = meshDict_.lookup("subsetFileName");
//         surfacePtr_->readFaceSubsets(runTime_.path()/subsetFileName);
//     }

    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}

/*
dualMeshGenerator::dualMeshGenerator
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

    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dualMeshGenerator::~dualMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualMeshGenerator::writeMesh() const
{
    # ifdef DEBUG
    mesh_.addressingData().checkCellsZipUp(true);
    mesh_.addressingData().checkCellVolumes(true);
    # endif

    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
