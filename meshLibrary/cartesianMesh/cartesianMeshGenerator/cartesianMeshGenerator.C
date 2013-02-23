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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "objectRegistry.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
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
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"

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
    
void cartesianMeshGenerator::createCartesianMesh()
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
	# else
	writeMeshEnsight(mesh_, "cartesianMesh");
	# endif
	//::exit(EXIT_FAILURE);
	# endif
}
	
void cartesianMeshGenerator::surfacePreparation()
{
	//- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
	//- It also checks topology of cells after morphing is performed
/*	do
	{
		surfaceMorpherCells* cmPtr = new surfaceMorpherCells(mesh_);
		cmPtr->morphMesh();
		deleteDemandDrivenData(cmPtr);
	} while( topologicalCleaner(mesh_).cleanTopology() );
    */
    
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
	# else
	writeMeshEnsight(mesh_, "afterTopoCleaning");
	# endif
	//::exit(EXIT_FAILURE);
	# endif
}
		
void cartesianMeshGenerator::mapMeshToSurface()
{
	//- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);
    
    //- pre-map mesh surface
    meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();
    
    # ifdef DEBUG
	# ifdef DEBUGfpma
    writeMeshFPMA(mesh_, "preMappedMesh");
    # else
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

void cartesianMeshGenerator::mapEdgesAndCorners()
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

void cartesianMeshGenerator::optimiseMeshSurface()
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
		
void cartesianMeshGenerator::generateBoudaryLayers()
{
	boundaryLayers bl(mesh_);
    
    if( meshDict_.found("boundaryLayers") )
    {
        wordList createLayers(meshDict_.lookup("boundaryLayers"));
        
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
		
void cartesianMeshGenerator::optimiseFinalMesh()
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

void cartesianMeshGenerator::replaceBoundaries()
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

void cartesianMeshGenerator::renumberMesh()
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
	
void cartesianMeshGenerator::generateMesh()
{
	createCartesianMesh();
	
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

// Construct from objectRegistry
cartesianMeshGenerator::cartesianMeshGenerator(const Time& time)
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
	
	if( meshDict_.found("subsetFileName") )
	{
		fileName subsetFileName = meshDict_.lookup("subsetFileName");
        if( Pstream::parRun() )
			subsetFileName = ".."/subsetFileName;
		surfacePtr_->readFaceSubsets(db_.path()/subsetFileName);
	}

    octreePtr_ = new meshOctree(*surfacePtr_);
	
	meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}

/*
cartesianMeshGenerator::cartesianMeshGenerator
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

cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshGenerator::writeMesh() const
{	
	mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
