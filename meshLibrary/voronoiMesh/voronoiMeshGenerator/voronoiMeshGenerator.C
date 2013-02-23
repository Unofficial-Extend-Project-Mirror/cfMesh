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

#include "voronoiMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "objectRegistry.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "voronoiMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "decomposeCellsNearConcaveEdges.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
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

void voronoiMeshGenerator::createVoronoiMesh()
{	
	//- create voronoi mesh from octree and Delaunay tets
	voronoiMeshExtractor vme(*octreePtr_, meshDict_, mesh_);

	vme.createMesh();
	
	# ifdef DEBUG
	mesh_.write();
	# ifdef DEBUGfpma
	writeMeshFPMA(mesh_, "voronoiMesh");
	# else
	writeMeshEnsight(mesh_, "voronoiMesh");
	# endif
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
		
void voronoiMeshGenerator::mapMeshToSurface()
{
	//- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

	//- map mesh surface on the geometry surface		
	meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();
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
	::exit(EXIT_FAILURE);
	# endif
    deleteDemandDrivenData(msePtr);
}

void voronoiMeshGenerator::mapEdgesAndCorners()
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

void voronoiMeshGenerator::optimiseMeshSurface()
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

void voronoiMeshGenerator::decomposeConcaveCells()
{
	decomposeCellsNearConcaveEdges dm(mesh_);
	
	# ifdef DEBUG
	mesh_.write();
	//meshOptimizer(*octreePtr_, mesh_).preOptimize();
	# ifdef DEBUGfpma
	writeMeshFPMA(mesh_, "decomposedConcaveCells");
	# else
	writeMeshEnsight(mesh_, "decomposedConcaveCells");
	#endif
	//::exit(EXIT_FAILURE);
	# endif
}
		
void voronoiMeshGenerator::generateBoudaryLayers()
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
	writeMeshEnsight(mesh_, "meshWithBndLayer");
	mesh_.write();
	//::exit(EXIT_FAILURE);
	# endif
}
		
void voronoiMeshGenerator::optimiseFinalMesh()
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

void voronoiMeshGenerator::replaceBoundaries()
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

void voronoiMeshGenerator::renumberMesh()
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
	
void voronoiMeshGenerator::generateMesh()
{
	createVoronoiMesh();
	
	surfacePreparation();
	
	mapMeshToSurface();
	
	mapEdgesAndCorners();
	
	optimiseMeshSurface();
	
	decomposeConcaveCells();

	generateBoudaryLayers();
	
	optimiseFinalMesh();
	
	renumberMesh();
    
    replaceBoundaries();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Time
voronoiMeshGenerator::voronoiMeshGenerator
(
    const Time& time
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
	
	if( meshDict_.found("subsetFileName") )
	{
		const fileName subsetFileName = meshDict_.lookup("subsetFileName");
		surfacePtr_->readFaceSubsets(runTime_.path()/subsetFileName);
	}

    octreePtr_ = new meshOctree(*surfacePtr_);
	
	meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

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
