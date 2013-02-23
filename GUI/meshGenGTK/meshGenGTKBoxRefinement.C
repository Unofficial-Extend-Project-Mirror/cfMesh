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
#include "boxRefinement.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- callback functions
	
extern "C"
{
	
static void showDataForBoxRefinement
(
	GtkWidget* widget,
	GtkWidget* boxRefinementFramePtr
)
{
	GtkWidget* comboboxentry_boxRefinement_createdBoxes =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"comboboxentry_boxRefinement_createdBoxes"
		);
	GtkWidget* entry_boxRefinement_selectedCellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCellSize"
		);
	GtkWidget* entry_boxRefinement_selectedCentreX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreX"
		);
	GtkWidget* entry_boxRefinement_selectedCentreY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreY"
		);
	GtkWidget* entry_boxRefinement_selectedCentreZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreZ"
		);
	GtkWidget* entry_boxRefinement_selectedLengthX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthX"
		);
	GtkWidget* entry_boxRefinement_selectedLengthY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthY"
		);
	GtkWidget* entry_boxRefinement_selectedLengthZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthZ"
		);
	
	//- find the object in the dictionary
	const word boxName =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_boxRefinement_createdBoxes)->entry
			)
		);
	
	const PtrList<entry> refs = guiPtr->objectRefinements();
	
	forAll(refs, refI)
		if( refs[refI].keyword() == boxName )
		{
			const dictionary dict = refs[refI].dict();
			
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedCellSize),
				help::scalarToText
				(
					scalar(readScalar(dict.lookup("cellSize")))
				).c_str()
			);
			const vector centre(dict.lookup("centre"));
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedCentreX),
				help::scalarToText(centre.x()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedCentreY),
				help::scalarToText(centre.y()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedCentreZ),
				help::scalarToText(centre.z()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedLengthX),
				help::scalarToText(readScalar(dict.lookup("lengthX"))).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedLengthY),
				help::scalarToText(readScalar(dict.lookup("lengthY"))).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_boxRefinement_selectedLengthZ),
				help::scalarToText(readScalar(dict.lookup("lengthZ"))).c_str()
			);
			
			return;
		}
		
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedCellSize), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedCentreX), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedCentreY), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedCentreZ), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedLengthX), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedLengthY), "");
	gtk_entry_set_text(GTK_ENTRY(entry_boxRefinement_selectedLengthZ), "");
}

static void resetValuesForSelectedBox
(
	GtkWidget* widget,
	GtkWidget* boxRefinementFramePtr
)
{
	GtkWidget* comboboxentry_boxRefinement_createdBoxes =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"comboboxentry_boxRefinement_createdBoxes"
		);
	GtkWidget* entry_boxRefinement_selectedCellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCellSize"
		);
	GtkWidget* entry_boxRefinement_selectedCentreX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreX"
		);
	GtkWidget* entry_boxRefinement_selectedCentreY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreY"
		);
	GtkWidget* entry_boxRefinement_selectedCentreZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedCentreZ"
		);
	GtkWidget* entry_boxRefinement_selectedLengthX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthX"
		);
	GtkWidget* entry_boxRefinement_selectedLengthY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthY"
		);
	GtkWidget* entry_boxRefinement_selectedLengthZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_selectedLengthZ"
		);
		
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO
				(
					comboboxentry_boxRefinement_createdBoxes
				)->entry
			)
		);
	
	if( name == "" )
		return;
	
	const PtrList<entry> refs = guiPtr->objectRefinements();
	dictionary dict;
	forAll(refs, refI)
		if( refs[refI].keyword() == name )
			dict = refs[refI].dict();
	guiPtr->removeObjectRefinement(name);
	if( widget == entry_boxRefinement_selectedCellSize )
	{
		const word size =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedCellSize));
		dict.remove("cellSize");
		dict.add("cellSize", help::textToScalar(size));
	}
	else if( widget == entry_boxRefinement_selectedCentreX )
	{
		point centre(dict.lookup("centre"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedCentreX));
		dict.remove("centre");
		centre.x() = help::textToScalar(c);
		dict.add("centre", centre);
	}
	else if( widget == entry_boxRefinement_selectedCentreY )
	{
		point centre(dict.lookup("centre"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedCentreY));
		dict.remove("centre");
		centre.y() = help::textToScalar(c);
		dict.add("centre", centre);
	}
	else if( widget == entry_boxRefinement_selectedCentreZ )
	{
		point centre(dict.lookup("centre"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedCentreZ));
		dict.remove("centre");
		centre.z() = help::textToScalar(c);
		dict.add("centre", centre);
	}
	else if( widget == entry_boxRefinement_selectedLengthX )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedLengthX));
		dict.remove("lengthX");
		dict.add("lengthX", help::textToScalar(l));
	}
	else if( widget == entry_boxRefinement_selectedLengthY )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedLengthY));
		dict.remove("lengthY");
		dict.add("lengthY", help::textToScalar(l));
	}
	else if( widget == entry_boxRefinement_selectedLengthZ )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_boxRefinement_selectedLengthZ));
		dict.remove("lengthZ");
		dict.add("lengthZ", help::textToScalar(l));
	}
	
	boxRefinement br(name, dict);
	guiPtr->addObjectRefinement(br);
}

static void updateBoxRefinements(GtkWidget* boxRefinementFramePtr)
{
	GtkWidget* comboboxentry_boxRefinement_createdBoxes =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"comboboxentry_boxRefinement_createdBoxes"
		);
	
	const PtrList<entry> refs = guiPtr->objectRefinements();
	label nBoxes(0);
	GList* createdBoxes = 0;
	forAll(refs, refI)
	{
		const word& key = refs[refI].keyword();
		const dictionary& dict = refs[refI].dict();
		const word type(dict.lookup("type"));
		
		if( type == "box" )
		{
			++nBoxes;
			boxRefinement br(key, dict);
			createdBoxes =
				g_list_append
				(
					createdBoxes,
					const_cast<char*>(br.name().c_str())
				);
		}
	}
		
	gtk_combo_set_popdown_strings
	(
		GTK_COMBO
		(
			comboboxentry_boxRefinement_createdBoxes
		), createdBoxes
	);
	
	if( nBoxes == 0 )
	{
		gtk_entry_set_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_boxRefinement_createdBoxes)->entry
			),
			""
		);
	}
		
	g_list_free(createdBoxes);
	
	showDataForBoxRefinement(NULL, boxRefinementFramePtr);
}

static void addBoxRefinement
(
	GtkWidget* button,
	GtkWidget* boxRefinementFramePtr
)
{
	GtkWidget* name =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxName"
		);
	GtkWidget* cellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_cellSize"
		);
	GtkWidget* centreX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_centreX"
		);
	GtkWidget* centreY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_centreY"
		);
	GtkWidget* centreZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_centreZ"
		);
	GtkWidget* lengthX =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_lengthX"
		);
	GtkWidget* lengthY =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_lengthY"
		);
	GtkWidget* lengthZ =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"entry_boxRefinement_lengthZ"
		);
	
	//- create all necessary settings
	const word boxName = gtk_entry_get_text(GTK_ENTRY(name));
	gtk_entry_set_text(GTK_ENTRY(name), "");
	const scalar boxCellSize =
		help::textToScalar(gtk_entry_get_text(GTK_ENTRY(cellSize)));
	gtk_entry_set_text(GTK_ENTRY(cellSize), "");
	point centre;
	centre.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreX)));
	gtk_entry_set_text(GTK_ENTRY(centreX), "");
	centre.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreY)));
	gtk_entry_set_text(GTK_ENTRY(centreY), "");
	centre.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreZ)));
	gtk_entry_set_text(GTK_ENTRY(centreZ), "");
	const scalar boxLengthX =
		help::textToScalar(gtk_entry_get_text(GTK_ENTRY(lengthX)));
	gtk_entry_set_text(GTK_ENTRY(lengthX), "");
	const scalar boxLengthY =
		help::textToScalar(gtk_entry_get_text(GTK_ENTRY(lengthY)));
	gtk_entry_set_text(GTK_ENTRY(lengthY), "");
	const scalar boxLengthZ =
		help::textToScalar(gtk_entry_get_text(GTK_ENTRY(lengthZ)));
	gtk_entry_set_text(GTK_ENTRY(lengthZ), "");
		
	boxRefinement box
	(
		boxName,
		boxCellSize,
		centre,
		boxLengthX,
		boxLengthY,
		boxLengthZ
	);
	
	guiPtr->addObjectRefinement(box);
	
	updateBoxRefinements(boxRefinementFramePtr);
}

static void removeBoxRefinement
(
	GtkWidget* button,
	GtkWidget* boxRefinementFramePtr
)
{
	GtkWidget* createdBoxes =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(boxRefinementFramePtr),
			"comboboxentry_boxRefinement_createdBoxes"
		);
	
	const word boxName =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(createdBoxes)->entry
			)
		);

	guiPtr->removeObjectRefinement(boxName);

	updateBoxRefinements(boxRefinementFramePtr);
}

}

void meshGenGTK::createBoxRefinementWindowPage()
{
	guiPtr = &meshGui_;
	
	boxRefinementFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a table
	GtkWidget* table1 = gtk_table_new(5, 11, FALSE);
	gtk_widget_show(table1);
	gtk_container_add
	(
		GTK_CONTAINER(boxRefinementFramePtr_),
		table1
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for box name
	GtkWidget* label_boxName = gtk_label_new("Box name");
	gtk_widget_show(label_boxName);
	gtk_table_attach(GTK_TABLE(table1), label_boxName, 0, 1, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boxName), 0, 0.5);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create an entry for box name
	GtkWidget* entry_boxName = gtk_entry_new();
	gtk_widget_show(entry_boxName);
	
	gtk_table_attach(GTK_TABLE(table1), entry_boxName, 1, 4, 0, 1,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for desired cell size
	GtkWidget* label_boxRefinement_cellSize = gtk_label_new("Cell size [m]");
	gtk_widget_show(label_boxRefinement_cellSize);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_cellSize,
					0, 1, 1, 2,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boxRefinement_cellSize), 0, 0.5);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for desired cell size
	GtkWidget* entry_boxRefinement_cellSize = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_cellSize);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_cellSize,
					1, 4, 1, 2,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create labels for axes
	GtkWidget* label_axes = gtk_label_new("X");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					1, 2, 2, 3,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
	
	label_axes = gtk_label_new("Y");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					2, 3, 2, 3,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
	
	label_axes = gtk_label_new("Z");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					3, 4, 2, 3,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for box centre
	GtkWidget* label_boxRefinement_centre = gtk_label_new("Centre");
	gtk_widget_show(label_boxRefinement_centre);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_centre,
					0, 1, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boxRefinement_centre), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for box centre
	GtkWidget* entry_boxRefinement_centreX = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_centreX);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_centreX,
					1, 2, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_centreY = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_centreY);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_centreY,
					2, 3, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_centreZ = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_centreZ);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_centreZ,
					3, 4, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for side lengths
	GtkWidget* label_boxRefinement_lengths = gtk_label_new("Side lengths");
	gtk_widget_show(label_boxRefinement_lengths);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_lengths,
					0, 1, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boxRefinement_lengths), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for side lengths
	GtkWidget* entry_boxRefinement_lengthX = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_lengthX);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_lengthX,
					1, 2, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_lengthY = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_lengthY);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_lengthY,
					2, 3, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_lengthZ = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_lengthZ);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_lengthZ,
					3, 4, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the button for adding boxes
	GtkWidget* button1 = gtk_button_new();
	gtk_widget_show(button1);
	gtk_table_attach
	(
		GTK_TABLE(table1),
		button1,
		4, 5, 0, 5,
		(GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 5
	);
	
	g_signal_connect
	(
		G_OBJECT(button1),
		"clicked",
		G_CALLBACK(addBoxRefinement),
		boxRefinementFramePtr_
	);

	GtkWidget* alignment1 = gtk_alignment_new(0.5, 0.5, 0, 0);
	gtk_widget_show(alignment1);
	gtk_container_add(GTK_CONTAINER(button1), alignment1);
	
	GtkWidget* hbox1 = gtk_hbox_new(FALSE, 2);
	gtk_widget_show(hbox1);
	gtk_container_add(GTK_CONTAINER(alignment1), hbox1);
	
	GtkWidget* image1 =
		gtk_image_new_from_stock("gtk-add", GTK_ICON_SIZE_BUTTON);
	gtk_widget_show(image1);
	gtk_box_pack_start(GTK_BOX(hbox1), image1, FALSE, FALSE, 0);
	
	GtkWidget* label3 = gtk_label_new_with_mnemonic("Add");
	gtk_widget_show(label3);
	gtk_box_pack_start(GTK_BOX(hbox1), label3, FALSE, FALSE, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- horizontal separator
	GtkWidget* hseparator0 = gtk_hseparator_new();
	gtk_widget_show(hseparator0);
	gtk_table_attach(GTK_TABLE(table1), hseparator0, 0, 5, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for existing box objects
	GtkWidget* label_boxRefinement_existing = gtk_label_new("Existing boxes");
	gtk_widget_show(label_boxRefinement_existing);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_existing,
					0, 1, 6, 7,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- combobox for existing box objects
	GtkWidget* comboboxentry_boxRefinement_createdBoxes =
		gtk_combo_new();
	gtk_widget_show(comboboxentry_boxRefinement_createdBoxes);
	gtk_table_attach(GTK_TABLE(table1),comboboxentry_boxRefinement_createdBoxes,
					1, 2, 6, 7,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
	gtk_editable_set_editable
	(
		GTK_EDITABLE
		(
			GTK_COMBO(comboboxentry_boxRefinement_createdBoxes)->entry
		),
		0
	);
	
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for selected cell size
	GtkWidget* label_boxRefinement_selectedCellSize =
		gtk_label_new("Cell size [m]");
	gtk_widget_show(label_boxRefinement_selectedCellSize);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_selectedCellSize,
					0, 1, 7, 8,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment
	(
		GTK_MISC(label_boxRefinement_selectedCellSize), 0, 0.5
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for selected cell size
	GtkWidget* entry_boxRefinement_selectedCellSize = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedCellSize);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedCellSize,
					1, 4, 7, 8,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create labels for axes
	label_axes = gtk_label_new("X");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					1, 2, 8, 9,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
	
	label_axes = gtk_label_new("Y");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					2, 3, 8, 9,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
	
	label_axes = gtk_label_new("Z");
	gtk_widget_show(label_axes);
	gtk_table_attach(GTK_TABLE(table1), label_axes,
					3, 4, 8, 9,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for box centre
	GtkWidget* label_boxRefinement_selectedCentre = gtk_label_new("Centre");
	gtk_widget_show(label_boxRefinement_selectedCentre);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_selectedCentre,
					0, 1, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_boxRefinement_selectedCentre), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for selected box centre
	GtkWidget* entry_boxRefinement_selectedCentreX = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedCentreX);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedCentreX,
					1, 2, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_selectedCentreY = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedCentreY);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedCentreY,
					2, 3, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_selectedCentreZ = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedCentreZ);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedCentreZ,
					3, 4, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for selected side lengths
	GtkWidget* label_boxRefinement_selectedLengths =
		gtk_label_new("Side lengths");
	gtk_widget_show(label_boxRefinement_selectedLengths);
	gtk_table_attach(GTK_TABLE(table1), label_boxRefinement_selectedLengths,
					0, 1, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment
	(
		GTK_MISC(label_boxRefinement_selectedLengths), 0, 0.5
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for side lengths
	GtkWidget* entry_boxRefinement_selectedLengthX = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedLengthX);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedLengthX,
					1, 2, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_selectedLengthY = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedLengthY);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedLengthY,
					2, 3, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_boxRefinement_selectedLengthZ = gtk_entry_new();
	gtk_widget_show(entry_boxRefinement_selectedLengthZ);
	gtk_table_attach(GTK_TABLE(table1), entry_boxRefinement_selectedLengthZ,
					3, 4, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the button for removing patches
	GtkWidget* button2 = gtk_button_new();
	gtk_widget_show(button2);
	gtk_table_attach
	(
		GTK_TABLE(table1),
		button2,
		4, 5, 6, 11,
		(GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 5
	);
	
	g_signal_connect
	(
		G_OBJECT(button2),
		"clicked",
		G_CALLBACK(removeBoxRefinement),
		boxRefinementFramePtr_
	);
	
	GtkWidget* alignment2 = gtk_alignment_new(0.5, 0.5, 0, 0);
	gtk_widget_show(alignment2);
	gtk_container_add(GTK_CONTAINER(button2), alignment2);
	
	GtkWidget* hbox2 = gtk_hbox_new(FALSE, 2);
	gtk_widget_show(hbox2);
	gtk_container_add(GTK_CONTAINER(alignment2), hbox2);
	
	GtkWidget* image2 =
		gtk_image_new_from_stock("gtk-remove", GTK_ICON_SIZE_BUTTON);
	gtk_widget_show(image2);
	gtk_box_pack_start(GTK_BOX (hbox2), image2, FALSE, FALSE, 0);
	
	GtkWidget* label5 = gtk_label_new_with_mnemonic("Remove");
	gtk_widget_show(label5);
	gtk_box_pack_start(GTK_BOX (hbox2), label5, FALSE, FALSE, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	/* Store pointers to all widgets, for use by lookup_widget(). */
	GLADE_HOOKUP_OBJECT_NO_REF
	(
		boxRefinementFramePtr_,
		boxRefinementFramePtr_,
		"boxRefinementFramePtr_");
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		table1,
		"table1"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxName,
		"entry_boxName"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_cellSize,
		"entry_boxRefinement_cellSize"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_centreX,
		"entry_boxRefinement_centreX"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_centreY,
		"entry_boxRefinement_centreY"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_centreZ,
		"entry_boxRefinement_centreZ"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_lengthX,
		"entry_boxRefinement_lengthX"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_lengthY,
		"entry_boxRefinement_lengthY"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_lengthZ,
		"entry_boxRefinement_lengthZ"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		comboboxentry_boxRefinement_createdBoxes,
		"comboboxentry_boxRefinement_createdBoxes"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedCellSize,
		"entry_boxRefinement_selectedCellSize"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedCentreX,
		"entry_boxRefinement_selectedCentreX"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedCentreY,
		"entry_boxRefinement_selectedCentreY"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedCentreZ,
		"entry_boxRefinement_selectedCentreZ"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedLengthX,
		"entry_boxRefinement_selectedLengthX"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedLengthY,
		"entry_boxRefinement_selectedLengthY"
	);
	GLADE_HOOKUP_OBJECT
	(
		boxRefinementFramePtr_,
		entry_boxRefinement_selectedLengthZ,
		"entry_boxRefinement_selectedLengthZ"
	);
	
	updateBoxRefinements(boxRefinementFramePtr_);
	
	g_signal_connect
	(
		G_OBJECT
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_boxRefinement_createdBoxes)->entry
			)
		),
		"changed",
		G_CALLBACK(showDataForBoxRefinement),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedCellSize),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedCentreX),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedCentreY),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedCentreZ),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedLengthX),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedLengthY),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_boxRefinement_selectedLengthZ),
		"changed",
		G_CALLBACK(resetValuesForSelectedBox),
		boxRefinementFramePtr_
	);
	
	gtk_widget_show(boxRefinementFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
