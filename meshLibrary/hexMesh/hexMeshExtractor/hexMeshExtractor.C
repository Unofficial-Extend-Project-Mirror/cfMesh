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

#include "hexMeshExtractor.H"
#include "meshOctree.H"

//#define DEBUGDual

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
void hexMeshExtractor::clearOut()
{
	deleteDemandDrivenData(centreNodeLabelPtr_);
    deleteDemandDrivenData(subVerticesPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
hexMeshExtractor::hexMeshExtractor
(
    const meshOctree& octree,
    const IOdictionary& dict,
    polyMeshGen& mesh
)
:
    octreeAddressing_(octree, dict),
    mesh_(mesh),
    centreNodeLabelPtr_(NULL),
    subVerticesPtr_(NULL),
    octreeVertexType_()
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hexMeshExtractor::~hexMeshExtractor()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void hexMeshExtractor::createMesh()
{
    Info << "Extracting hex mesh" << endl;
    
    classifyOctreePoints();
	
    createPoints();

    createHexMesh();
	
    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).reorderBoundaryFaces();
	
    Info << "Mesh has :" << nl
	 << mesh_.points().size() << " vertices " << nl
	 << mesh_.faces().size() << " faces" << nl
	 << mesh_.cells().size() << " cells" << endl;
	
    Info << "Finished extracting hex template" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
