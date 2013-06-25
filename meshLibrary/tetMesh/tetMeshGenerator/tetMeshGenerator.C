/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright(C) 2005-2007 Franjo Juretic
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or(at your
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

#include "tetMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "meshOctreeAutomaticRefinement.H"
#include "tetMeshExtractorOctree.H"
#include "delaunayTessellation.H"
#include "subdivisionTessellation.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "triSurfacePatchManipulator.H"

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

void tetMeshGenerator::createTetMesh()
{
    //- create tet Mesh from octree and Delaunay tets
    tetMeshExtractorOctree tme(*octreePtr_, meshDict_, mesh_);

    tme.createMesh();

# ifdef DEBUG
    mesh_.write();
# ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "tetMesh");
# else
    writeMeshEnsight(mesh_, "tetMesh");
# endif
    ::exit(EXIT_FAILURE);
# endif
}

void tetMeshGenerator::surfacePreparation()
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

void tetMeshGenerator::mapMeshToSurface()
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
# ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterSurfaceSmoothing");
# else
    writeMeshEnsight(mesh_, "afterSurfaceSmoothing");
# endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
# endif
    deleteDemandDrivenData(msePtr);
}

void tetMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);

# ifdef DEBUG
    mesh_.write();
    //meshOptimizer(*octreePtr_, mesh_).preOptimize();
# ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "withEdges");
# else
    writeMeshEnsight(mesh_, "withEdges");
#endif
    //::exit(EXIT_FAILURE);
# endif
}

void tetMeshGenerator::optimiseMeshSurface()
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

void tetMeshGenerator::generateBoudaryLayers()
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

# ifdef DEBUG
    writeMeshEnsight(mesh_, "meshWithBndLayer");
    mesh_.write();
    //::exit(EXIT_FAILURE);
# endif
}

void tetMeshGenerator::optimiseFinalMesh()
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

void tetMeshGenerator::replaceBoundaries()
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

void tetMeshGenerator::renumberMesh()
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

void tetMeshGenerator::generateMesh()
{
    createTetMesh();

    surfacePreparation();

    mapMeshToSurface();

    mapEdgesAndCorners();

    optimiseMeshSurface();

    generateBoudaryLayers();

    optimiseFinalMesh();

    renumberMesh();

    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Time
tetMeshGenerator::tetMeshGenerator
(
    const Time& time
)
:
    runTime_(time),
    surfacePtr_(NULL),
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
    octreePtr_(NULL),
    mesh_(time)
{
    if( true )
    {
        checkMeshDict cmd(meshDict_);
    }

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
tetMeshGenerator::tetMeshGenerator
(
    const Time& time,
    const volScalarField& localCellSize
)
:
    runTime_(time),
    surfacePtr_(NULL),
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
    octreePtr_(NULL),
    mesh_(time)
{
    fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurface(runTime_.path()/surfaceFile);

    octreePtr_ = new meshOctree(*surfacePtr_);

    generateMesh();
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshGenerator::~tetMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
