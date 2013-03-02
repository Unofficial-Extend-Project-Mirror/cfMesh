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

#include "dualMeshExtractor.H"
#include "meshOctree.H"

//#define DEBUGDual

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const direction dualMeshExtractor::faceFlip_[6][4] =
    {
        {0, 4, 6, 2},
        {5, 1, 3, 7},
        {0, 2, 3, 1},
        {4, 5, 7, 6},
        {2, 6, 7, 3},
        {0, 1, 5, 4}
    };
    
void dualMeshExtractor::clearOut()
{
    deleteDemandDrivenData(centreNodeLabelPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
dualMeshExtractor::dualMeshExtractor
(
    const meshOctree& octree,
    const IOdictionary& dict,
    polyMeshGen& mesh
)
:
    octreeAddressing_(octree, dict, true),
    mesh_(mesh),
    centreNodeLabelPtr_(NULL)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dualMeshExtractor::~dualMeshExtractor()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualMeshExtractor::createMesh()
{
    Info << "Extracting polyMesh" << endl;
    
    createPoints();

    createPolyMesh();
    
    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).reorderBoundaryFaces();
    
    Info << "Mesh has :" << nl
        << mesh_.points().size() << " vertices " << nl
        << mesh_.faces().size() << " faces" << nl
        << mesh_.cells().size() << " cells" << endl;
    
    # ifdef DEBUGDual
    Info << "Points start at address " << long(&mesh_.points()) << endl;
    Info << "Faces start at address " << long(&mesh_.faces()) << endl;
    Info << "Cells start at address " << long(&mesh_.cells()) << endl;
    # endif
    
    Info << "Finished extracting polyMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
