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

#include "triangulateBoundaryLayerBaseFaces.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triangulateBoundaryLayerBaseFaces::triangulateBoundaryLayerBaseFaces
(
    polyMeshGen& mesh,
    const VRWGraph& layerCellsInColumn
)
:
    mesh_(mesh),
    layerCellsInColumn_(layerCellsInColumn),
    invertedCell_(mesh_.cells().size(), false),
    ownerInColumn_(),
    neiInColumn_(),
    faceType_(),
    faceCentreLabel_()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triangulateBoundaryLayerBaseFaces::~triangulateBoundaryLayerBaseFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateBoundaryLayerBaseFaces::setBadCells(const labelHashSet& s)
{
    forAllConstIter(labelHashSet, s, it)
        invertedCell_[it.key()] = true;
}

void triangulateBoundaryLayerBaseFaces::triangulateLayers()
{
    classifyMeshFaces();

    if( !createNewPoints() )
    {
        WarningIn
        (
            "void triangulateBoundaryLayerBaseFaces::triangulateLayers()"
        ) << "Boundary layer is not modified" << endl;

        return;
    }

    createNewFacesAndCells();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
