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

#include "correctEdgesBetweenPatches.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "decomposeCells.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& correctEdgesBetweenPatches::meshSurface() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

//- delete mesh surface
void correctEdgesBetweenPatches::clearMeshSurface()
{
    deleteDemandDrivenData(msePtr_);
}

void correctEdgesBetweenPatches::replaceBoundary()
{
    clearMeshSurface();

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
}

void correctEdgesBetweenPatches::decomposeCorrectedCells()
{
    if( decompose_ )
    {
        clearMeshSurface();

        decomposeCells dc(mesh_);
        dc.decomposeMesh(decomposeCell_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, octree, regions for boundary vertices
correctEdgesBetweenPatches::correctEdgesBetweenPatches
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    msePtr_(NULL),
    patchNames_(mesh.boundaries().size()),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_(),
    decomposeCell_(mesh_.cells().size(), false),
    decompose_(false)
{
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
        patchNames_[patchI] = boundaries[patchI].patchName();

    decomposeProblematicFaces();

    patchCorrection();

    decomposeCorrectedCells();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

correctEdgesBetweenPatches::~correctEdgesBetweenPatches()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
