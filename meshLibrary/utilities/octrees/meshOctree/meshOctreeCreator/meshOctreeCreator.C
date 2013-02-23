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

#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "IOdictionary.H"
#include "boundBox.H"
#include "demandDrivenData.H"


// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshOctreeCreator::meshOctreeCreator(meshOctree& mo)
:
    octree_(mo),
    meshDictPtr_(NULL),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0))
{}

meshOctreeCreator::meshOctreeCreator
(
    meshOctree& mo,
    const IOdictionary& dict
)
:
    octree_(mo),
    meshDictPtr_(&dict),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0))
{}

/*
meshOctreeCreator::meshOctreeCreator
(
    meshOctree& mo,
    const IOdictionary& dict,
    const volScalarField& localCellSize
)
:
    octree_(mo),
    meshDict_(dict),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0))
{
    FatalErrorIn
    (
        "meshOctreeCreator::meshOctreeCreator"
        "("
        "meshOctree& mo,"
        "const IOdictionary& dict,"
        "const volScalarField& localCellSize"
        ")"
    ) << "Not implemented!" << exit(FatalError);

    Info << "Constructing octree" << endl;

    Info << "Finished constructing octree" << endl;
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctreeCreator::~meshOctreeCreator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCreator::activateHexRefinement()
{
    hexRefinement_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
