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

#include "boundaryFacesGenerator.H"
#include "boundaryOutwardOrientation.H"

//#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private member functions * * * * * * * * * * * //

void boundaryFacesGenerator::createBoundaryFaces()
{
    //- set sizes of boundary cells
    polyMeshGenModifier meshModifier(mesh_);
    cellListPMG& cells = meshModifier.cellsAccess();
    forAll(boundaryCell_, cellI)
        if( boundaryCell_[cellI] )
            cells[cellI].setSize(nFacesInCell_[cellI]);

    //- create boundary faces
    forAll(boundaryCell_, cellI)
        if( boundaryCell_[cellI] )
            createBoundaryFacesForCell(cellI);
}

void boundaryFacesGenerator::storeBoundaryFaces()
{
    wordList patchNames;
    if( hasFacesInDefaultPatch_ )
    {
        patchNames.setSize(surface_.patches().size() + 1);
        patchNames[surface_.patches().size()] = "defaultFaces";
    }
    else
    {
        patchNames.setSize(surface_.patches().size());
    }

    forAll(surface_.patches(), patchI)
        patchNames[patchI] = surface_.patches()[patchI].name();

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        boundaryFaces_,
        boundaryOwners_,
        boundaryPatches_
    );
    polyMeshGenModifier(mesh_).removeUnusedVertices();
}

void boundaryFacesGenerator::createSurfaceCornersAndPatches()
{
    const VRWGraph& pointFaces = surface_.pointFacets();
    forAll(pointFaces, pI)
    {
        DynList<label> patches;
        forAllRow(pointFaces, pI, pfI)
            patches.appendIfNotIn(surface_[pointFaces(pI, pfI)].region());

        //- corners are points in at least three regions
        if( patches.size() > 2 )
        {
            surfaceCorners_.append(pI);
            cornersPatches_.append(patches);
        }
    }

    Info << "The surface has " << surfaceCorners_.size() << " corners" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

boundaryFacesGenerator::boundaryFacesGenerator
(
    polyMeshGen& mesh,
    List<direction>& nFacesInCell,
    label& nPoints,
    const label nIntFaces,
    const boolList& boundaryCell,
    VRWGraph& pointPatches,
    const triSurf& surface
)
:
    surface_(surface),
    mesh_(mesh),
    nFacesInCell_(nFacesInCell),
    boundaryCell_(boundaryCell),
    pointRegions_(pointPatches),
    nPoints_(nPoints),
    nIntFaces_(nIntFaces),
    boundaryFaces_(),
    boundaryOwners_(),
    boundaryPatches_(),
    hasFacesInDefaultPatch_(false),
    defaultPatchID_(surface.patches().size()),
    surfaceCorners_(surface.patches().size()),
    cornersPatches_(surface.patches().size())
{
    createSurfaceCornersAndPatches();

    createBoundaryFaces();

    storeBoundaryFaces();

    boundaryOutwardOrientation boo(mesh);
}

// Destructor
boundaryFacesGenerator::~boundaryFacesGenerator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
