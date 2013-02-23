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

#include "meshGenGTK.H"
#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
meshGenGTK::meshGenGTK
(
    label argc, char* argv[], const objectRegistry& reg
)
:
	meshGui_(reg),
	generalPageFramePtr_(NULL),
	mainWindowPtr_(NULL),
	localRefinementFramePtr_(NULL),
	keepCellsIntersectingPatchesFramePtr_(NULL),
	bndLayersFramePtr_(NULL),
	localSettingsMainFramePtr_(NULL),
	
	boxRefinementFramePtr_(NULL),
	lineRefinementFramePtr_(NULL),
	coneRefinementFramePtr_(NULL),
	sphereRefinementFramePtr_(NULL),
	objectRefinementMainFramePtr_(NULL),
	
	renameBoundaryFramePtr_(NULL)
{
	gtk_init(&argc, &argv);
	
	createLocalSettingsMainWindowPage();
	createGeneralPage();
	createObjectRefinementMainWindowPage();
	createRenameBoundaryMainWindowPage();
	createMainWindowPage();
	
	gtk_main();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshGenGTK::~meshGenGTK()
{
	meshGui_.writeDict();
    gtk_main_quit();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
