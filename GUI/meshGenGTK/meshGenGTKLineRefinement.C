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
#include "lineRefinement.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- callback functions

extern "C"
{

static void showDataForLineRefinement
(
	GtkWidget* widget,
	GtkWidget* lineRefinementFramePtr
)
{
	GtkWidget* comboboxentry_lineRefinement_createdLines =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"comboboxentry_lineRefinement_createdLines"
		);
	GtkWidget* entry_lineRefinement_selectedCellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedCellSize"
		);
	GtkWidget* entry_lineRefinement_selectedP0X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0X"
		);
	GtkWidget* entry_lineRefinement_selectedP0Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0Y"
		);
	GtkWidget* entry_lineRefinement_selectedP0Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0Z"
		);
	GtkWidget* entry_lineRefinement_selectedP1X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1X"
		);
	GtkWidget* entry_lineRefinement_selectedP1Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1Y"
		);
	GtkWidget* entry_lineRefinement_selectedP1Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1Z"
		);
	
	//- find the object in the dictionary
	const word lineName =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_lineRefinement_createdLines)->entry
			)
		);
	
	const PtrList<entry> refs = guiPtr->objectRefinements();
	
	forAll(refs, refI)
		if( refs[refI].keyword() == lineName )
		{
			const dictionary dict = refs[refI].dict();
			
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedCellSize),
				help::scalarToText
				(
					scalar(readScalar(dict.lookup("cellSize")))
				).c_str()
			);
			const vector p0(dict.lookup("p0"));
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP0X),
				help::scalarToText(p0.x()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP0Y),
				help::scalarToText(p0.y()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP0Z),
				help::scalarToText(p0.z()).c_str()
			);
			
			const vector p1(dict.lookup("p1"));
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP1X),
				help::scalarToText(p1.x()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP1Y),
				help::scalarToText(p1.y()).c_str()
			);
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_lineRefinement_selectedP1Z),
				help::scalarToText(p1.z()).c_str()
			);
			
			return;
		}
		
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedCellSize), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP0X), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP0Y), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP0Z), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP1X), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP1Y), "");
	gtk_entry_set_text(GTK_ENTRY(entry_lineRefinement_selectedP1Z), "");
}

static void resetValuesForSelectedLine
(
	GtkWidget* widget,
	GtkWidget* lineRefinementFramePtr
)
{
	GtkWidget* comboboxentry_lineRefinement_createdLines =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"comboboxentry_lineRefinement_createdLines"
		);
	GtkWidget* entry_lineRefinement_selectedCellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedCellSize"
		);
	GtkWidget* entry_lineRefinement_selectedP0X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0X"
		);
	GtkWidget* entry_lineRefinement_selectedP0Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0Y"
		);
	GtkWidget* entry_lineRefinement_selectedP0Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP0Z"
		);
	GtkWidget* entry_lineRefinement_selectedP1X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1X"
		);
	GtkWidget* entry_lineRefinement_selectedP1Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1Y"
		);
	GtkWidget* entry_lineRefinement_selectedP1Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_selectedP1Z"
		);
		
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO
				(
					comboboxentry_lineRefinement_createdLines
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
	if( widget == entry_lineRefinement_selectedCellSize )
	{
		const word size =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedCellSize));
		dict.remove("cellSize");
		dict.add("cellSize", help::textToScalar(size));
	}
	else if( widget == entry_lineRefinement_selectedP0X )
	{
		point p0(dict.lookup("p0"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP0X));
		dict.remove("p0");
		p0.x() = help::textToScalar(c);
		dict.add("p0", p0);
	}
	else if( widget == entry_lineRefinement_selectedP0Y )
	{
		point p0(dict.lookup("p0"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP0Y));
		dict.remove("p0");
		p0.y() = help::textToScalar(c);
		dict.add("p0", p0);
	}
	else if( widget == entry_lineRefinement_selectedP0Z )
	{
		point p0(dict.lookup("p0"));
		const word c =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP0Z));
		dict.remove("p0");
		p0.z() = help::textToScalar(c);
		dict.add("p0", p0);
	}
	else if( widget == entry_lineRefinement_selectedP1X )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP1X));
		point p1(dict.lookup("p1"));
		dict.remove("p1");
		p1.x() = help::textToScalar(l);
		dict.add("p1", p1);
	}
	else if( widget == entry_lineRefinement_selectedP1Y )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP1Y));
		point p1(dict.lookup("p1"));
		dict.remove("p1");
		p1.y() = help::textToScalar(l);
		dict.add("p1", p1);
	}
	else if( widget == entry_lineRefinement_selectedP1Z )
	{
		const word l =
			gtk_entry_get_text(GTK_ENTRY(entry_lineRefinement_selectedP1Z));
		point p1(dict.lookup("p1"));
		dict.remove("p1");
		p1.z() = help::textToScalar(l);
		dict.add("p1", p1);
	}
	
	lineRefinement br(name, dict);
	guiPtr->addObjectRefinement(br);
}

static void updateLineRefinements(GtkWidget* lineRefinementFramePtr)
{
	GtkWidget* comboboxentry_lineRefinement_createdLines =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"comboboxentry_lineRefinement_createdLines"
		);
	
	const PtrList<entry> refs = guiPtr->objectRefinements();
	label nLines(0);
	GList* createdLines = 0;
	forAll(refs, refI)
	{
		const word& key = refs[refI].keyword();
		const dictionary& dict = refs[refI].dict();
		const word type(dict.lookup("type"));
		
		if( type == "line" )
		{
			++nLines;
			lineRefinement line(key, dict);
			createdLines =
				g_list_append
				(
					createdLines,
					const_cast<char*>(line.name().c_str())
				);
		}
	}
		
	gtk_combo_set_popdown_strings
	(
		GTK_COMBO
		(
			comboboxentry_lineRefinement_createdLines
		), createdLines
	);
	
	if( nLines == 0 )
	{
		gtk_entry_set_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_lineRefinement_createdLines)->entry
			),
			""
		);
	}
		
	g_list_free(createdLines);
	
	showDataForLineRefinement(NULL, lineRefinementFramePtr);
}

static void addLineRefinement
(
	GtkWidget* button,
	GtkWidget* lineRefinementFramePtr
)
{
	GtkWidget* name =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineName"
		);
	GtkWidget* cellSize =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_cellSize"
		);
	GtkWidget* p0X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p0X"
		);
	GtkWidget* p0Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p0Y"
		);
	GtkWidget* p0Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p0Z"
		);
	GtkWidget* p1X =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p1X"
		);
	GtkWidget* p1Y =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p1Y"
		);
	GtkWidget* p1Z =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"entry_lineRefinement_p1Z"
		);
	
	//- create all necessary settings
	const word lineName = gtk_entry_get_text(GTK_ENTRY(name));
	gtk_entry_set_text(GTK_ENTRY(name), "");
	const scalar lineCellSize =
		help::textToScalar(gtk_entry_get_text(GTK_ENTRY(cellSize)));
	gtk_entry_set_text(GTK_ENTRY(cellSize), "");
	point p0, p1;
	p0.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0X)));
	gtk_entry_set_text(GTK_ENTRY(p0X), "");
	p0.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0Y)));
	gtk_entry_set_text(GTK_ENTRY(p0Y), "");
	p0.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0Z)));
	gtk_entry_set_text(GTK_ENTRY(p0Z), "");
	p1.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1X)));
	gtk_entry_set_text(GTK_ENTRY(p1X), "");
	p1.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1Y)));
	gtk_entry_set_text(GTK_ENTRY(p1Y), "");
	p1.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1Z)));
	gtk_entry_set_text(GTK_ENTRY(p1Z), "");
		
	lineRefinement line
	(
		lineName,
		lineCellSize,
		p0,
		p1
	);
	
	guiPtr->addObjectRefinement(line);
	
	updateLineRefinements(lineRefinementFramePtr);
}

static void removeLineRefinement
(
	GtkWidget* button,
	GtkWidget* lineRefinementFramePtr
)
{
	GtkWidget* createdLines =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(lineRefinementFramePtr),
			"comboboxentry_lineRefinement_createdLines"
		);
	
	const word lineName =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(createdLines)->entry
			)
		);

	guiPtr->removeObjectRefinement(lineName);

	updateLineRefinements(lineRefinementFramePtr);
}

}

void meshGenGTK::createLineRefinementWindowPage()
{
	guiPtr = &meshGui_;
	
	lineRefinementFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a table
	GtkWidget* table1 = gtk_table_new(5, 11, FALSE);
	gtk_widget_show(table1);
	gtk_container_add
	(
		GTK_CONTAINER(lineRefinementFramePtr_),
		table1
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for box name
	GtkWidget* label_lineName = gtk_label_new("Line name");
	gtk_widget_show(label_lineName);
	gtk_table_attach(GTK_TABLE(table1), label_lineName, 0, 1, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_lineName), 0, 0.5);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create an entry for box name
	GtkWidget* entry_lineName = gtk_entry_new();
	gtk_widget_show(entry_lineName);
	
	gtk_table_attach(GTK_TABLE(table1), entry_lineName, 1, 4, 0, 1,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for desired cell size
	GtkWidget* label_lineRefinement_cellSize = gtk_label_new("Cell size [m]");
	gtk_widget_show(label_lineRefinement_cellSize);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_cellSize,
					0, 1, 1, 2,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_lineRefinement_cellSize), 0, 0.5);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for desired cell size
	GtkWidget* entry_lineRefinement_cellSize = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_cellSize);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_cellSize,
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
	//- label for p0
	GtkWidget* label_lineRefinement_p0 = gtk_label_new("Starting point");
	gtk_widget_show(label_lineRefinement_p0);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_p0,
					0, 1, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_lineRefinement_p0), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for p0
	GtkWidget* entry_lineRefinement_p0X = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p0X);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p0X,
					1, 2, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_p0Y = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p0Y);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p0Y,
					2, 3, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_p0Z = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p0Z);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p0Z,
					3, 4, 3, 4,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for p1
	GtkWidget* label_lineRefinement_p1 = gtk_label_new("End point");
	gtk_widget_show(label_lineRefinement_p1);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_p1,
					0, 1, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_lineRefinement_p1), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for p1
	GtkWidget* entry_lineRefinement_p1X = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p1X);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p1X,
					1, 2, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_p1Y = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p1Y);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p1Y,
					2, 3, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_p1Z = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_p1Z);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_p1Z,
					3, 4, 4, 5,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the button for adding lines
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
		G_CALLBACK(addLineRefinement),
		lineRefinementFramePtr_
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
	//- label for existing line objects
	GtkWidget* label_lineRefinement_existing = gtk_label_new("Existing lines");
	gtk_widget_show(label_lineRefinement_existing);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_existing,
					0, 1, 6, 7,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- combobox for existing box objects
	GtkWidget* comboboxentry_lineRefinement_createdLines =
		gtk_combo_new();
	gtk_widget_show(comboboxentry_lineRefinement_createdLines);
	gtk_table_attach
	(
		GTK_TABLE(table1),comboboxentry_lineRefinement_createdLines,
		1, 2, 6, 7,
		(GtkAttachOptions)(GTK_FILL),
		(GtkAttachOptions)(0), 0, 0
	);
					
	gtk_editable_set_editable
	(
		GTK_EDITABLE
		(
			GTK_COMBO(comboboxentry_lineRefinement_createdLines)->entry
		),
		0
	);
	
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for selected cell size
	GtkWidget* label_lineRefinement_selectedCellSize =
		gtk_label_new("Cell size [m]");
	gtk_widget_show(label_lineRefinement_selectedCellSize);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_selectedCellSize,
					0, 1, 7, 8,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment
	(
		GTK_MISC(label_lineRefinement_selectedCellSize), 0, 0.5
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a label for selected cell size
	GtkWidget* entry_lineRefinement_selectedCellSize = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedCellSize);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedCellSize,
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
	//- label for p0
	GtkWidget* label_lineRefinement_selectedP0 =
		gtk_label_new("Starting point");
	gtk_widget_show(label_lineRefinement_selectedP0);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_selectedP0,
					0, 1, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label_lineRefinement_selectedP0), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for selected p0
	GtkWidget* entry_lineRefinement_selectedP0X = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP0X);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP0X,
					1, 2, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_selectedP0Y = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP0Y);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP0Y,
					2, 3, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_selectedP0Z = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP0Z);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP0Z,
					3, 4, 9, 10,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- label for selected p1
	GtkWidget* label_lineRefinement_selectedP1 =
		gtk_label_new("End point");
	gtk_widget_show(label_lineRefinement_selectedP1);
	gtk_table_attach(GTK_TABLE(table1), label_lineRefinement_selectedP1,
					0, 1, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	gtk_misc_set_alignment
	(
		GTK_MISC(label_lineRefinement_selectedP1), 0, 0.5
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- entries for p1
	GtkWidget* entry_lineRefinement_selectedP1X = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP1X);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP1X,
					1, 2, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_selectedP1Y = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP1Y);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP1Y,
					2, 3, 10, 11,
					(GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	GtkWidget* entry_lineRefinement_selectedP1Z = gtk_entry_new();
	gtk_widget_show(entry_lineRefinement_selectedP1Z);
	gtk_table_attach(GTK_TABLE(table1), entry_lineRefinement_selectedP1Z,
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
		G_CALLBACK(removeLineRefinement),
		lineRefinementFramePtr_
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
		lineRefinementFramePtr_,
		lineRefinementFramePtr_,
		"lineRefinementFramePtr_");
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		table1,
		"table1"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineName,
		"entry_lineName"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_cellSize,
		"entry_lineRefinement_cellSize"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p0X,
		"entry_lineRefinement_p0X"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p0Y,
		"entry_lineRefinement_p0Y"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p0Z,
		"entry_lineRefinement_p0Z"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p1X,
		"entry_lineRefinement_p1X"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p1Y,
		"entry_lineRefinement_p1Y"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_p1Z,
		"entry_lineRefinement_p1Z"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		comboboxentry_lineRefinement_createdLines,
		"comboboxentry_lineRefinement_createdLines"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedCellSize,
		"entry_lineRefinement_selectedCellSize"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP0X,
		"entry_lineRefinement_selectedP0X"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP0Y,
		"entry_lineRefinement_selectedP0Y"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP0Z,
		"entry_lineRefinement_selectedP0Z"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP1X,
		"entry_lineRefinement_selectedP1X"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP1Y,
		"entry_lineRefinement_selectedP1Y"
	);
	GLADE_HOOKUP_OBJECT
	(
		lineRefinementFramePtr_,
		entry_lineRefinement_selectedP1Z,
		"entry_lineRefinement_selectedP1Z"
	);
	
	updateLineRefinements(lineRefinementFramePtr_);
	
	g_signal_connect
	(
		G_OBJECT
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_lineRefinement_createdLines)->entry
			)
		),
		"changed",
		G_CALLBACK(showDataForLineRefinement),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedCellSize),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP0X),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP0Y),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP0Z),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP1X),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP1Y),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_lineRefinement_selectedP1Z),
		"changed",
		G_CALLBACK(resetValuesForSelectedLine),
		lineRefinementFramePtr_
	);
	
	gtk_widget_show(lineRefinementFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
