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

#include "tetMeshExtractor.H"
#include "tetTessellation.H"
#include "meshOctree.H"
#include "triSurf.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
tetMeshExtractor::tetMeshExtractor
(
    const tetTessellation& tessellation,
    const meshOctree& octree,
    polyMeshGen& mesh
)
        :
        tessellation_ ( tessellation ),
        octree_ ( octree ),
        mesh_ ( mesh ),
        useBoundaryTets_ ( false ),
        patchesForBoundaryTet_ ( 0 ),
        useElement_ ( tessellation.elmts().size(), false )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshExtractor::~tetMeshExtractor()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshExtractor::useTetsIntersectingBoundary()
{
    useBoundaryTets_ = true;
}

void tetMeshExtractor::useTetsIntersectingBoundaryPatches
(
    const wordList& patchNames
)
{
    patchesForBoundaryTet_.setSize ( patchNames.size() );

    const triSurf& ts = octree_.surface();

    forAll ( patchNames, nameI )
    forAll ( ts.patches(), patchI )
    if ( ts.patches() [patchI].name() == patchNames[nameI] )
    {
        patchesForBoundaryTet_[nameI] = patchI;
    }
}

void tetMeshExtractor::createMesh()
{
    Info << "Extracting tetMesh" << endl;

    //- selecting which elements will be used as mesh cells
    selectElements();

    //- create points and pointLeaves addressing
    createPoints();

    //- create the mesh
    createPolyMesh();

    polyMeshGenModifier ( mesh_ ).reorderBoundaryFaces();
    polyMeshGenModifier ( mesh_ ).removeUnusedVertices();

    Info << "Mesh has :" << nl
    << mesh_.points().size() << " vertices " << nl
    << mesh_.faces().size() << " faces" << nl
    << mesh_.cells().size() << " cells" << endl;

    Info << "Finished extracting tetMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
