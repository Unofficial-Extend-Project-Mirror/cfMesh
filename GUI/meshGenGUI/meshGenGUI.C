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

#include "meshGenGUI.H"
#include "triSurface.H"
#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void meshGenGUI::createSurface(const fileName& name) const
{
	surfacePtr_ = new triSurface(name);
}

const triSurface& meshGenGUI::surface() const
{
	if( !surfacePtr_ )
		createSurface(surfaceFileName());
	
	return *surfacePtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
meshGenGUI::meshGenGUI
(
    const objectRegistry& reg
)
:
	meshDict_
	(
		IOobject
        (
            "meshDict",
            reg.time().system(),
            reg,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
	surfacePtr_(NULL)
{
	Info << meshDict_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshGenGUI::~meshGenGUI()
{
    deleteDemandDrivenData(surfacePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
