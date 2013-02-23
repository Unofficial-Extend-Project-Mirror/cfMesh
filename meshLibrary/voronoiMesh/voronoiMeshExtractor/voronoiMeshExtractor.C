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

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

// #define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::clearOut()
{
	deleteDemandDrivenData(pointEdgesPtr_);
	deleteDemandDrivenData(edgesPtr_);
	deleteDemandDrivenData(edgeTetsPtr_);
	deleteDemandDrivenData(boundaryEdgePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
voronoiMeshExtractor::voronoiMeshExtractor
(
    const meshOctree& octree,
	const IOdictionary& meshDict,
	polyMeshGen& mesh
)
:
    tetCreator_(octree, meshDict),
	mesh_(mesh),
	pointEdgesPtr_(NULL),
	edgesPtr_(NULL),
	edgeTetsPtr_(NULL),
	boundaryEdgePtr_(NULL)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voronoiMeshExtractor::~voronoiMeshExtractor()
{
	clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createMesh()
{
	Info << "Extracting voronoi mesh" << endl;
	
	//- copy tet points into the mesh
	createPoints();

	//- create the mesh
    createPolyMesh();
	
	polyMeshGenModifier(mesh_).reorderBoundaryFaces();
	polyMeshGenModifier(mesh_).removeUnusedVertices();
	
	Info << "Mesh has :" << nl
		<< mesh_.points().size() << " vertices " << nl
		<< mesh_.faces().size() << " faces" << nl
		<< mesh_.cells().size() << " cells" << endl;
	
    Info << "Finished extracting voronoi mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
