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
