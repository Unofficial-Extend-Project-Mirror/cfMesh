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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGTK::createLocalSettingsMainWindowPage()
{
	//- create pages within this page
	createLocalRefinementWindowPage();
	createKeepCellsIntersectingPatchesWindowPage();
	createBndLayersWindowPage();
	
	localSettingsMainFramePtr_ = gtk_frame_new(NULL);
	
	GtkWidget* notebook1 = gtk_notebook_new ();

	gtk_container_add(GTK_CONTAINER(localSettingsMainFramePtr_), notebook1);
	
	//- add page for local refinement
	GtkWidget* label_localRefinement = gtk_label_new("Local refinement");
	gtk_widget_show(label_localRefinement);
	
	gtk_notebook_append_page
	(
		GTK_NOTEBOOK(notebook1),
		localRefinementFramePtr_,
		label_localRefinement
	);
	
	//- add page for keepCellsIntersectingPatches
	GtkWidget* label_keepBoxesIntersectingPatches =
		gtk_label_new("Keep octree boxes intersecting patches");
	gtk_widget_show(label_keepBoxesIntersectingPatches);

	gtk_notebook_append_page
	(
		GTK_NOTEBOOK(notebook1),
		keepCellsIntersectingPatchesFramePtr_,
		label_keepBoxesIntersectingPatches
	);
	
	//- add page for bnd layers
	GtkWidget* label_bndLayers =
		gtk_label_new("Boundary layer for patches");
	gtk_widget_show(label_bndLayers);
	
	gtk_notebook_append_page
	(
		GTK_NOTEBOOK(notebook1),
		bndLayersFramePtr_,
		label_bndLayers
	);
	
	//- show the notebook
	gtk_widget_show(notebook1);
	
	/* Store pointers to all widgets, for use by lookup_widget(). */
	GLADE_HOOKUP_OBJECT_NO_REF
	(
		localSettingsMainFramePtr_,
		localSettingsMainFramePtr_,
		"localSettingsMainFramePtr_"
	);
	GLADE_HOOKUP_OBJECT
	(
		localSettingsMainFramePtr_,
		notebook1,
		"notebook1"
	);
	GLADE_HOOKUP_OBJECT
	(
		localSettingsMainFramePtr_,
		label_localRefinement,
		"label_localRefinement"
	);
	GLADE_HOOKUP_OBJECT
	(
		localSettingsMainFramePtr_,
		label_keepBoxesIntersectingPatches,
		"label_keepBoxesIntersectingPatches"
	);
	GLADE_HOOKUP_OBJECT
	(
		localSettingsMainFramePtr_,
		label_bndLayers,
		"label_bndLayers"
	);

	gtk_widget_show(localSettingsMainFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
