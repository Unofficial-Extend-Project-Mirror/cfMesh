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
#include "helperFunctions.H"
#include "gtkHelpers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// callback functions
	
extern "C"
{
	
static void useButtonCallback(GtkWidget* widget, gpointer data)
{
	if( (gchar*)data == "UseBndSize" )
	{
		if( guiPtr->boundaryCellSizeEntryExist() ) 
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 0);
			guiPtr->removeBoundaryCellSize();
		}
		else
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 1);
			guiPtr->setBoundaryCellSize
			(
				guiPtr->maxCellSize()
			);
		}
	}
	else if( (gchar*)data == "UseAutoRef" )
	{
		if( guiPtr->minCellSizeEntryExist() )
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 0);
			guiPtr->removeMinCellSize();
		}
		else
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 1);
			guiPtr->setMinCellSize(guiPtr->maxCellSize());
		}
	}
	else if( (gchar*)data == "UseKeepCellsIntersectingBoundary" )
	{
		if( guiPtr->keepCellsIntersectingBoundaryEntryExist() )
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 0);
			guiPtr->removeKeepCellsIntersectingBoundary();
			
			if( guiPtr->checkForGluedMeshEntryExist() )
			{
				GtkWidget* checkButtonGluedMesh =
					(GtkWidget*)g_object_get_data
					(
						G_OBJECT(widget),
						"checkbutton_checkForGluedMesh"
					);
				
				gtk_toggle_button_set_active
				(
					GTK_TOGGLE_BUTTON(checkButtonGluedMesh),
					0
				);
				guiPtr->removeCheckForGluedMesh();
			}
		}
		else
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 1);
			guiPtr->setKeepCellsIntersectingBoundary();
		}
	}
	else if( (gchar*)(data) == "UseCheckForGluedMesh" )
	{
		if( !guiPtr->keepCellsIntersectingBoundary() )
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 0);
			return;
		}
		
		if( guiPtr->checkForGluedMeshEntryExist() )
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 0);
			guiPtr->removeCheckForGluedMesh();
		}
		else
		{
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), 1);
			guiPtr->setCheckForGluedMesh();
		}
	}
}

static void enterSurfaceFile(GtkWidget* widget, GtkEntry* entry)
{
	const fileName file = gtk_entry_get_text(entry);
	guiPtr->setSurfaceFileName(file);
}

static void enterMaxCellSize(GtkWidget* widget, GtkEntry* entry)
{
	const word s = gtk_entry_get_text(GTK_ENTRY(widget));
	const scalar size = help::textToScalar(s);
	guiPtr->setMaxCellSize(size);
}

static void enterBndCellSize(GtkWidget* widget, GtkEntry* entry)
{
	if( !guiPtr->boundaryCellSizeEntryExist() )
		return;
	const word s = gtk_entry_get_text(GTK_ENTRY(widget));
	const scalar size = help::textToScalar(s);
	guiPtr->setBoundaryCellSize(size);
}

static void enterMinCellSize(GtkWidget* widget, GtkEntry* entry)
{
	if( !guiPtr->minCellSizeEntryExist() )
		return;
	const word s = gtk_entry_get_text(GTK_ENTRY(widget));
	const scalar size = help::textToScalar(s);
	guiPtr->setMinCellSize(size);
}

static void editableEntry(GtkWidget *checkbutton, GtkWidget *entry)
{
	gtk_editable_set_editable
	(
		GTK_EDITABLE(entry),
		GTK_TOGGLE_BUTTON(checkbutton)->active
	);
	
	if( !GTK_TOGGLE_BUTTON(checkbutton)->active )
		gtk_entry_set_text(GTK_ENTRY(entry), "");
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGTK::createGeneralPage()
{
	guiPtr = &meshGui_;
	
	//- create a window for a general page
	generalPageFramePtr_ = gtk_frame_new(NULL);
	
	GtkWidget* table1 = gtk_table_new (11, 5, FALSE);
	gtk_widget_show(table1);
	gtk_container_add(GTK_CONTAINER(generalPageFramePtr_), table1);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the surface file
	GtkWidget* label_surfaceFile = gtk_label_new("Surface file");
	gtk_widget_show(label_surfaceFile);
	gtk_table_attach(GTK_TABLE(table1), label_surfaceFile, 1, 2, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC (label_surfaceFile), 0, 0.5);

	GtkWidget* label_Space = gtk_label_new("");
	gtk_widget_show(label_Space);
	gtk_table_attach(GTK_TABLE(table1), label_Space, 2, 3, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_widget_set_size_request(label_Space, 10, -1);
	gtk_misc_set_alignment(GTK_MISC(label_Space), 0, 0.5);
	
	//- create entry box for the surface file
	GtkWidget* entry_surfaceFile = gtk_entry_new();
	gtk_widget_show(entry_surfaceFile);
	gtk_table_attach(GTK_TABLE(table1), entry_surfaceFile, 3, 4, 0, 1,
                    (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	gtk_entry_set_text
	(
		GTK_ENTRY(entry_surfaceFile),
		meshGui_.surfaceFileName().c_str()
	);
	
	gtk_editable_set_editable
	(
		GTK_EDITABLE(entry_surfaceFile),
		1
	);
	
	g_signal_connect(G_OBJECT(entry_surfaceFile), "changed",
		      G_CALLBACK(enterSurfaceFile),
		      (gpointer)entry_surfaceFile);

/*	GtkWidget* image1 =
		gtk_image_new_from_stock("gtk-open", GTK_ICON_SIZE_BUTTON);
	gtk_widget_show(image1);
	gtk_table_attach(GTK_TABLE (table1), image1, 4, 5, 0, 1,
                    (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions)(GTK_FILL), 0, 0);
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- box for maxCellSize
	GtkWidget* label_maxCellSize = gtk_label_new("Max cell size [m]");
	gtk_widget_show(label_maxCellSize);
	gtk_table_attach(GTK_TABLE(table1), label_maxCellSize, 1, 2, 1, 2,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_maxCellSize), 0, 0.5);

	//- entry for maxCellSize
	GtkWidget* entry_maxCellSize = gtk_entry_new();
	gtk_widget_show(entry_maxCellSize);
	gtk_table_attach(GTK_TABLE(table1), entry_maxCellSize, 3, 5, 1, 2,
                    (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	gtk_entry_set_text
	(
		GTK_ENTRY(entry_maxCellSize),
		help::scalarToText(meshGui_.maxCellSize()).c_str()
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_maxCellSize), "changed",
		G_CALLBACK(enterMaxCellSize),
		(gpointer)entry_maxCellSize
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- horizontal separator
	GtkWidget* hseparator0 = gtk_hseparator_new();
	gtk_widget_show(hseparator0);
	gtk_table_attach(GTK_TABLE(table1), hseparator0, 0, 5, 2, 3,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- boundary cell size

	//- checkbutton for boundaryCellSize
	GtkWidget* checkbutton_boundaryCellSize =
		gtk_check_button_new_with_mnemonic
		(
			"Use different size for boundary cells"
		);
	gtk_widget_show(checkbutton_boundaryCellSize);
	gtk_table_attach(GTK_TABLE(table1),
					checkbutton_boundaryCellSize, 0, 5, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 1, 0);

	
	//- add callbacks for the button
	if( meshGui_.boundaryCellSizeEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_boundaryCellSize),
			1
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_boundaryCellSize),
			0
		);
	}
	
	g_signal_connect
	(
		gpointer(checkbutton_boundaryCellSize),
		"clicked",
		G_CALLBACK(useButtonCallback),
		(gpointer)"UseBndSize"
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //	
	//- create a box for boundary cell size
	GtkWidget* label_boundaryCellSize = gtk_label_new("Boundary cell size [m]");
	gtk_widget_show(label_boundaryCellSize);
	gtk_table_attach(GTK_TABLE(table1), label_boundaryCellSize, 1, 2, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boundaryCellSize), 0, 0.5);
	
	GtkWidget* entry_boundaryCellSize = gtk_entry_new();
	gtk_widget_show(entry_boundaryCellSize);
	gtk_table_attach(GTK_TABLE (table1), entry_boundaryCellSize, 3, 5, 4, 5,
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions) (0), 1, 0);
	gtk_entry_set_invisible_char (GTK_ENTRY (entry_boundaryCellSize), 9679);
	
	gtk_entry_set_visibility
	(
		GTK_ENTRY(entry_boundaryCellSize),
		meshGui_.boundaryCellSizeEntryExist()
	);
	
	if( meshGui_.boundaryCellSizeEntryExist() )
	{
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_boundaryCellSize),
			help::scalarToText(meshGui_.boundaryCellSize()).c_str()
		);
	}
	else
	{
		gtk_editable_set_editable(GTK_EDITABLE(entry_boundaryCellSize), 0);
	}
	
	g_signal_connect
	(
		G_OBJECT(checkbutton_boundaryCellSize),
		"toggled",
		G_CALLBACK(editableEntry),
		(gpointer)entry_boundaryCellSize
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boundaryCellSize),
		"changed",
		G_CALLBACK(enterBndCellSize),
		(gpointer)entry_boundaryCellSize
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- horizontal separator
	GtkWidget* hseparator1 = gtk_hseparator_new();
	gtk_widget_show(hseparator1);
	gtk_table_attach(GTK_TABLE(table1), hseparator1, 0, 5, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a box for automatic refinement
	
	//- checkbox for automatic refinement
	GtkWidget* checkbutton_automaticRefinement =
		gtk_check_button_new_with_mnemonic("Use automatic refinement");
	gtk_widget_show(checkbutton_automaticRefinement);
	gtk_table_attach(GTK_TABLE(table1), checkbutton_automaticRefinement,
					0, 5, 6, 7,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	if( meshGui_.minCellSizeEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_automaticRefinement), 1
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_automaticRefinement), 0
		);
	}
	
	g_signal_connect
	(
		gpointer(checkbutton_automaticRefinement),
		"clicked",
		G_CALLBACK(useButtonCallback),
		(gpointer)"UseAutoRef"
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //	
	//- create a box for min cell size
	GtkWidget* label_minCellSize = gtk_label_new("Min cell size [m]");
	gtk_widget_show(label_minCellSize);
	gtk_table_attach(GTK_TABLE(table1), label_minCellSize, 1, 2, 7, 8,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (0), 0, 0);
	gtk_misc_set_alignment (GTK_MISC (label_minCellSize), 0, 0.5);
	
	GtkWidget* entry_minCellSize = gtk_entry_new();
	gtk_widget_show(entry_minCellSize);
	gtk_table_attach(GTK_TABLE (table1), entry_minCellSize, 3, 5, 7, 8,
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions) (0), 1, 0);
	
	if( meshGui_.minCellSizeEntryExist() )
	{
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_minCellSize),
			help::scalarToText(meshGui_.minCellSize()).c_str()
		);
	}
	else
	{
		gtk_editable_set_editable(GTK_EDITABLE(entry_minCellSize), 0);
	}
	
	g_signal_connect
	(
		G_OBJECT(checkbutton_automaticRefinement),
		"toggled",
		G_CALLBACK(editableEntry),
		(gpointer)entry_minCellSize
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_minCellSize),
		"changed",
		G_CALLBACK(enterMinCellSize),
		(gpointer)entry_minCellSize
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- horizontal separator
	GtkWidget* hseparator2 = gtk_hseparator_new();
	gtk_widget_show(hseparator2);
	gtk_table_attach(GTK_TABLE(table1), hseparator2, 0, 5, 8, 9,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- activate keepCellsIntersectingBoundary option
	GtkWidget* checkbutton_keepCellsIntersectingBoundary =
		gtk_check_button_new_with_mnemonic
		(
			"Keep octree boxes intersecting the boundary"
		);
	gtk_widget_show(checkbutton_keepCellsIntersectingBoundary);
	gtk_table_attach(GTK_TABLE(table1),
					checkbutton_keepCellsIntersectingBoundary,
					0, 5, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	if( meshGui_.keepCellsIntersectingBoundaryEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_keepCellsIntersectingBoundary)
			,1
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_keepCellsIntersectingBoundary)
			,0
		);
	}
	
	g_signal_connect
	(
		gpointer(checkbutton_keepCellsIntersectingBoundary),
		"clicked",
		G_CALLBACK(useButtonCallback),
		(gpointer)"UseKeepCellsIntersectingBoundary"
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- activate checkForGluedMesh
	GtkWidget* checkbutton_checkForGluedMesh =
		gtk_check_button_new_with_mnemonic("Remove invalid mesh connections");
	gtk_widget_show(checkbutton_checkForGluedMesh);
	gtk_table_attach(GTK_TABLE(table1), checkbutton_checkForGluedMesh,
					0, 5, 10, 11,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	if( !meshGui_.keepCellsIntersectingBoundary() )
		meshGui_.removeCheckForGluedMesh();
	
	if( meshGui_.checkForGluedMesh() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_checkForGluedMesh),
			1
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_checkForGluedMesh),
			0
		);
	}
	
	g_signal_connect
	(
		gpointer(checkbutton_checkForGluedMesh),
		"clicked",
		G_CALLBACK(useButtonCallback),
		(gpointer)"UseCheckForGluedMesh"
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	/* Store pointers to all widgets, for use by lookup_widget(). */
	GLADE_HOOKUP_OBJECT_NO_REF
	(
		generalPageFramePtr_,
		generalPageFramePtr_,
		"window_generalPage_"
	);
	GLADE_HOOKUP_OBJECT
	(
		checkbutton_keepCellsIntersectingBoundary,
		checkbutton_checkForGluedMesh,
		"checkbutton_checkForGluedMesh"
	);
	
	gtk_widget_show(generalPageFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
