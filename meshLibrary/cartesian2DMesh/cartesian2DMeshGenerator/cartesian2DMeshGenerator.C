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

#include "cartesian2DMeshGenerator.H"
#include "polyMeshGen2DEngine.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "objectRegistry.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEdgeExtractor2D.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"

//#define DEBUG
//#define DEBUGfpma
//#define DEBUGEnsight

# ifdef DEBUG
#include "writeMeshEnsight.H"
#include "writeMeshFPMA.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void cartesian2DMeshGenerator::createCartesianMesh()
{
    //- create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);

    if( meshDict_.found("decomposePolyhedraIntoTetsAndPyrs") )
    {
        if( readBool(meshDict_.lookup("decomposePolyhedraIntoTetsAndPyrs")) )
            cme.decomposeSplitHexes();
    }

    cme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "cartesianMesh");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "cartesianMesh");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void cartesian2DMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
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

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterTopoCleaning");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "afterTopoCleaning");
    # endif
    //::exit(EXIT_FAILURE);
    # endif
}

void cartesian2DMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- pre-map mesh surface
    meshSurfaceMapper2D mapper(*msePtr, *octreePtr_);

    mapper.adjustZCoordinates();

    mapper.preMapVertices();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "preMappedMesh");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "preMappedMesh");
    # endif
    mesh_.write();
    //::exit(EXIT_FAILURE);
    # endif

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "afterMapping");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "afterMapping");
    # endif
    mesh_.write();
    ::exit(0);
    # endif

    deleteDemandDrivenData(msePtr);
}

void cartesian2DMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractor2D(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    ::exit(0);
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "withEdges");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "withEdges");
    #endif
    ::exit(0);
    # endif
}

void cartesian2DMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer optimizer(mse, *octreePtr_);
    optimizer.optimizeSurface2D();
    optimizer.untangleSurface2D();

    # ifdef DEBUG
    mesh_.write();
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "optSurfaceWithEdges");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "optSurfaceWithEdges");
    #endif
    ::exit(0);
    # endif
}

void cartesian2DMeshGenerator::generateBoundaryLayers()
{
    boundaryLayers bl(mesh_);

    bl.activate2DMode();

    bl.addLayerForAllPatches();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "meshWithBndLayer");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "meshWithBndLayer");
    # endif
    mesh_.write();
    //::exit(0);
    # endif
}

void cartesian2DMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.activate2DMode();

        refLayers.refineLayers();

        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer optimizer(mse, *octreePtr_);

        optimizer.optimizeSurface2D();
        optimizer.untangleSurface2D();
    }
}

void cartesian2DMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_,"renamedPatchesMesh");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "renamedPatchesMesh");
    #endif
    # endif
}

void cartesian2DMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    # ifdef DEBUGfpma
    writeMeshFPMA(mesh_,"renumberedMesh");
    # elif DEBUGEnsight
    writeMeshEnsight(mesh_, "renumberedMesh");
    #endif
    # endif
}

void cartesian2DMeshGenerator::generateMesh()
{
    createCartesianMesh();

    surfacePreparation();

    mapMeshToSurface();

    mapEdgesAndCorners();

    optimiseMeshSurface();

    generateBoundaryLayers();

    optimiseMeshSurface();

    refBoundaryLayers();

    renumberMesh();

    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
cartesian2DMeshGenerator::cartesian2DMeshGenerator(const Time& time)
:
    db_(time),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.system(),
            db_,
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

    fileName surfaceFile = meshDict_.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    surfacePtr_ = new triSurf(db_.path()/surfaceFile);

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

    octreePtr_ = new meshOctree(*surfacePtr_, true);

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}

/*
cartesian2DMeshGenerator::cartesian2DMeshGenerator
(
    const objectRegistry& time,
    const volScalarField& localCellSize
)
:
    db_(time),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.time().constant(),
            db_,
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

    generateMesh();
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesian2DMeshGenerator::~cartesian2DMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesian2DMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
