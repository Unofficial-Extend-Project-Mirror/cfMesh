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
#include "triSurface.H"
#include "helperFunctions.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
extern "C"
{

static void resetNewPatchNamesForSelectedPatch
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* comboboxentry_renameBoundary_selectedPatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_selectedPatches"
		);
	GtkWidget* entry_renameBoundary_selectedPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_selectedPatchName"
		);
	GtkWidget* entry_renameBoundary_selectedPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_selectedPatchType"
		);
	
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_renameBoundary_selectedPatches)->entry
			)
		);
	
	const PtrList<entry> newPatchNames = guiPtr->newPatchNames();
	forAll(newPatchNames, nameI)
		if( newPatchNames[nameI].keyword() == name )
		{
			const dictionary dict = newPatchNames[nameI].dict();
			
			dictionary newDict;
			if( widget == entry_renameBoundary_selectedPatchName )
			{
				const word newName =
					gtk_entry_get_text
					(
						GTK_ENTRY(entry_renameBoundary_selectedPatchName)
					);
				word type("patch");
				if( dict.found("type") )
				{
					type = word(dict.lookup("type"));
				}
				newDict.add("newName", newName);
				newDict.add("type", type);
			}
			else
			{
				const word newName(dict.lookup("newName"));
				const word type =
					gtk_entry_get_text
					(
						GTK_ENTRY(entry_renameBoundary_selectedPatchType)
					);
				newDict.add("newName", newName);
				newDict.add("type", type);
			}
			
			guiPtr->removePatchName(name);
			guiPtr->addNewPatchName(name, newDict);
		}
}

static void updateNewPatchNamesForSelectedPatch
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* comboboxentry_renameBoundary_selectedPatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_selectedPatches"
		);
	GtkWidget* entry_renameBoundary_selectedPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_selectedPatchName"
		);
	GtkWidget* entry_renameBoundary_selectedPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_selectedPatchType"
		);
	
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_renameBoundary_selectedPatches)->entry
			)
		);
		
	if( name == "" )
	{
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_selectedPatchName), ""
		);
		
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_selectedPatchType), ""
		);
		
		return;
	}
	
	const PtrList<entry> newPatchNames = guiPtr->newPatchNames();
	forAll(newPatchNames, nameI)
		if( newPatchNames[nameI].keyword() == name )
		{
			const dictionary dict = newPatchNames[nameI].dict();
			
			const word newName(dict.lookup("newName"));
			word type("patch");
			if( dict.found("type") )
			{
				type = word(dict.lookup("type"));
			}
			
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_renameBoundary_selectedPatchName),
				newName.c_str()
			);
			
			gtk_entry_set_text
			(
				GTK_ENTRY(entry_renameBoundary_selectedPatchType),
				type.c_str()
			);
		}
}
	
static void updateNewPatchNames
(
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* comboboxentry_renameBoundary_availablePatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_availablePatches"
		);
	GtkWidget* comboboxentry_renameBoundary_selectedPatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_selectedPatches"
		);
	
	if( guiPtr->newPatchNamesEntryExist() )
	{
		wordHashSet selectedPatches;
		PtrList<entry> newPatchNames = guiPtr->newPatchNames();
		forAll(newPatchNames, pI)
			selectedPatches.insert(newPatchNames[pI].keyword());
		
		label nAvailable(0);
		GList* avPatches = 0;
		const triSurface& surf = guiPtr->surface();
		forAll(surf.patches(), patchI)
		{
			if( selectedPatches.found(surf.patches()[patchI].name()) )
				continue;
			
			++nAvailable;
			const char* name = surf.patches()[patchI].name().c_str();
			avPatches = g_list_append(avPatches, const_cast<char*>(name));
		}
		
		if( nAvailable == 0 )
			gtk_entry_set_text
			(
				GTK_ENTRY
				(
					GTK_COMBO
					(
						comboboxentry_renameBoundary_availablePatches
					)->entry
				),
				""
			);
		
		gtk_combo_set_popdown_strings
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_availablePatches
			), avPatches
		);
		
		g_list_free(avPatches);
		
		GList* selPatches = 0;
		forAll(newPatchNames, nameI)
		{
			const char* name = newPatchNames[nameI].keyword().c_str();
			selPatches = g_list_append(selPatches, const_cast<char*>(name));
		}
		
		gtk_combo_set_popdown_strings
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_selectedPatches
			), selPatches
		);
		
		g_list_free(selPatches);
		
		if( nAvailable == newPatchNames.size() )
		{
			gtk_entry_set_text
			(
				GTK_ENTRY
				(
					GTK_COMBO
					(
						comboboxentry_renameBoundary_selectedPatches
					)->entry
				),
				""
			);
		}
	}
	else
	{
		GList* avPatches = 0;
		
		const triSurface& surf = guiPtr->surface();
		forAll(surf.patches(), patchI)
		{
			const char* name = surf.patches()[patchI].name().c_str();
			avPatches = g_list_append(avPatches, const_cast<char*>(name));
		}
		
		gtk_combo_set_popdown_strings
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_availablePatches
			), avPatches
		);
		
		g_list_free(avPatches);
		
		gtk_combo_set_popdown_strings
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_selectedPatches
			), avPatches
		);
		
		gtk_entry_set_text
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_renameBoundary_selectedPatches)->entry
			),
			""
		);
	}
	
	updateNewPatchNamesForSelectedPatch(NULL, renameBoundaryFramePtr);
}

static void setDefaultValues(GtkWidget* renameBoundaryFramePtr)
{
	GtkWidget* checkbutton_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchName"
		);
	GtkWidget* entry_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchName"
		);
	GtkWidget* checkbutton_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchType"
		);
	GtkWidget* entry_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchType"
		);
	
	if( guiPtr->defaultPatchNameEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchName),
			1
		);
		
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchName), 1
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchName),
			guiPtr->defaultPatchName().c_str()
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchName),
			0
		);
		
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchName), 0
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchName),
			""
		);
	}
	
	if( guiPtr->defaultPatchTypeEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchType),
			1
		);
		
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchType), 1
		);
		
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchType),
			guiPtr->defaultPatchType().c_str()
		);
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchType),
			0
		);
		
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchType), 0
		);
		
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchType),
			""
		);
	}
}

static void useDefaultPatchName
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* checkbutton_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchName"
		);
	GtkWidget* entry_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchName"
		);
	
	if( guiPtr->defaultPatchNameEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchName), 0
		);
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchName), 0
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchName), ""
		);
		guiPtr->removeDefaultPatchName();
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchName), 1
		);
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchName), 1
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchName), "defaultFaces"
		);
		guiPtr->setDefaultPatchName("defaultFaces");
	}
}

static void setDefaultPatchName
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* checkbutton_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchName"
		);
	GtkWidget* entry_renameBoundary_defaultPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchName"
		);
	
	if(
		gtk_toggle_button_get_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchName)
		)
	)
	{
		const word name =
			gtk_entry_get_text
			(
				GTK_ENTRY(entry_renameBoundary_defaultPatchName)
			);
		
		guiPtr->setDefaultPatchName(name);
	}
}

static void useDefaultPatchType
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* checkbutton_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchType"
		);
	GtkWidget* entry_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchType"
		);
	
	if( guiPtr->defaultPatchTypeEntryExist() )
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchType), 0
		);
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchType), 0
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchType), ""
		);
		guiPtr->removeDefaultPatchType();
	}
	else
	{
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchType), 1
		);
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_defaultPatchType), 1
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_defaultPatchType), "patch"
		);
		guiPtr->setDefaultPatchType("patch");
	}
}

static void setDefaultPatchType
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* checkbutton_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_defaultPatchType"
		);
	GtkWidget* entry_renameBoundary_defaultPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_defaultPatchType"
		);
	
	if(
		gtk_toggle_button_get_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_defaultPatchType)
		)
	)
	{
		const word pType =
			gtk_entry_get_text
			(
				GTK_ENTRY(entry_renameBoundary_defaultPatchType)
			);
		
		guiPtr->setDefaultPatchType(pType);
	}
}

static void useNewPatchType
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* entry_renameBoundary_newPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_newPatchType"
		);
	
	Info << "Here " << gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)) << endl;
	if( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)) )
	{
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_newPatchType), 1
		);
	}
	else
	{
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_newPatchType), ""
		);
		gtk_editable_set_editable
		(
			GTK_EDITABLE(entry_renameBoundary_newPatchType), 0
		);
	}
}

static void addNewPatchName
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* comboboxentry_renameBoundary_availablePatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_availablePatches"
		);
	GtkWidget* entry_renameBoundary_newPatchName =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_newPatchName"
		);
	GtkWidget* checkbutton_renameBoundary_newPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"checkbutton_renameBoundary_newPatchType"
		);
	GtkWidget* entry_renameBoundary_newPatchType =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"entry_renameBoundary_newPatchType"
		);
		
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO
				(
					comboboxentry_renameBoundary_availablePatches
				)->entry
			)
		);
		
	if( name == "" )
		return;

	const word newName =
		gtk_entry_get_text
		(
			GTK_ENTRY(entry_renameBoundary_newPatchName)
		);
	gtk_entry_set_text(GTK_ENTRY(entry_renameBoundary_newPatchName), "");
	
	word newType("patch");
	if(
		gtk_toggle_button_get_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_newPatchType)
		)
	)
	{
		newType =
			gtk_entry_get_text
			(
				GTK_ENTRY(entry_renameBoundary_newPatchType)
			);
	
		gtk_toggle_button_set_active
		(
			GTK_TOGGLE_BUTTON(checkbutton_renameBoundary_newPatchType), 0
		);
		gtk_entry_set_text
		(
			GTK_ENTRY(entry_renameBoundary_newPatchType), ""
		);
	}
	
	dictionary dict;
	dict.add("newName", newName);
	dict.add("type", newType);
	guiPtr->addNewPatchName(name, dict);
	
	updateNewPatchNames(renameBoundaryFramePtr);
}

static void removeNewPatchName
(
	GtkWidget* widget,
	GtkWidget* renameBoundaryFramePtr
)
{
	GtkWidget* comboboxentry_renameBoundary_selectedPatches =
		(GtkWidget*)g_object_get_data
		(
			G_OBJECT(renameBoundaryFramePtr),
			"comboboxentry_renameBoundary_selectedPatches"
		);
	
	const word name =
		gtk_entry_get_text
		(
			GTK_ENTRY
			(
				GTK_COMBO
				(
					comboboxentry_renameBoundary_selectedPatches
				)->entry
			)
		);
		
	if( name == "" )
		return;
	
	guiPtr->removePatchName(name);
	
	updateNewPatchNames(renameBoundaryFramePtr);
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGTK::createRenameBoundaryMainWindowPage()
{
	guiPtr = &meshGui_;
	
	renameBoundaryFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a table
	GtkWidget* table1 = gtk_table_new(3, 10, FALSE);
	gtk_widget_show(table1);
	gtk_container_add
	(
		GTK_CONTAINER(renameBoundaryFramePtr_),
		table1
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a check button for default patch name
	GtkWidget* checkbutton_renameBoundary_defaultPatchName =
		gtk_check_button_new_with_mnemonic
		(
			"Default patch name"
		);
	gtk_widget_show(checkbutton_renameBoundary_defaultPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					checkbutton_renameBoundary_defaultPatchName, 0, 1, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create the entry for default patch name
	GtkWidget* entry_renameBoundary_defaultPatchName = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_defaultPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_defaultPatchName, 1, 2, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create a check button for default patch type
	GtkWidget* checkbutton_renameBoundary_defaultPatchType =
		gtk_check_button_new_with_mnemonic
		(
			"Default patch type"
		);
	gtk_widget_show(checkbutton_renameBoundary_defaultPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					checkbutton_renameBoundary_defaultPatchType, 0, 1, 1, 2,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- create the entry for default patch type
	GtkWidget* entry_renameBoundary_defaultPatchType = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_defaultPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_defaultPatchType, 1, 2, 1, 2,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a horizontal separator
	GtkWidget* hseparator = gtk_hseparator_new();
	gtk_widget_show(hseparator);
	gtk_table_attach(GTK_TABLE(table1), hseparator, 0, 3, 2, 3,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a label for available patches
	GtkWidget* label_renameBoundary_availablePatches =
		gtk_label_new("Available patches");
	gtk_widget_show(label_renameBoundary_availablePatches);
	
	gtk_table_attach(GTK_TABLE(table1),
					label_renameBoundary_availablePatches, 0, 1, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a combo box for available patches
	GtkWidget* comboboxentry_renameBoundary_availablePatches =
		gtk_combo_new();
	gtk_widget_show(comboboxentry_renameBoundary_availablePatches);
	
	gtk_table_attach(GTK_TABLE(table1),
					comboboxentry_renameBoundary_availablePatches, 1, 2, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
	gtk_editable_set_editable
	(
		GTK_EDITABLE
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_availablePatches
			)->entry
		),
		0
	);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a label for new patch name
	GtkWidget* label_renameBoundary_newPatchName =
		gtk_label_new("New name");
	gtk_widget_show(label_renameBoundary_newPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					label_renameBoundary_newPatchName, 0, 1, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add the entry for new patch name
	GtkWidget* entry_renameBoundary_newPatchName = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_newPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_newPatchName, 1, 2, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a check button for new patch type
	GtkWidget* checkbutton_renameBoundary_newPatchType =
		gtk_check_button_new_with_mnemonic
		(
			"New type"
		);
	gtk_widget_show(checkbutton_renameBoundary_newPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					checkbutton_renameBoundary_newPatchType, 0, 1, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add the entry for new patch type
	GtkWidget* entry_renameBoundary_newPatchType = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_newPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_newPatchType, 1, 2, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the button for adding new patch names
	GtkWidget* button1 = gtk_button_new();
	gtk_widget_show(button1);
	gtk_table_attach
	(
		GTK_TABLE(table1),
		button1,
		2, 3, 3, 6,
		(GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 3
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
	//- add a horizontal separator
	hseparator = gtk_hseparator_new();
	gtk_widget_show(hseparator);
	gtk_table_attach(GTK_TABLE(table1), hseparator, 0, 3, 6, 7,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
					0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a label for selected patches
	GtkWidget* label_renameBoundary_selectedPatches =
		gtk_label_new("Selected patches");
	gtk_widget_show(label_renameBoundary_selectedPatches);
	
	gtk_table_attach(GTK_TABLE(table1),
					label_renameBoundary_selectedPatches, 0, 1, 7, 8,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a combo box for selected patches
	GtkWidget* comboboxentry_renameBoundary_selectedPatches =
		gtk_combo_new();
	gtk_widget_show(comboboxentry_renameBoundary_selectedPatches);
	
	gtk_table_attach(GTK_TABLE(table1),
					comboboxentry_renameBoundary_selectedPatches, 1, 2, 7, 8,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
	
	gtk_editable_set_editable
	(
		GTK_EDITABLE
		(
			GTK_COMBO
			(
				comboboxentry_renameBoundary_selectedPatches
			)->entry
		),
		0
	);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a label for new patch name
	GtkWidget* label_renameBoundary_selectedPatchName =
		gtk_label_new("New name");
	gtk_widget_show(label_renameBoundary_selectedPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					label_renameBoundary_selectedPatchName, 0, 1, 8, 9,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add the entry for selected patch name
	GtkWidget* entry_renameBoundary_selectedPatchName = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_selectedPatchName);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_selectedPatchName, 1, 2, 8, 9,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
					
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add a label for selected patch type
	GtkWidget* label_renameBoundary_selectedPatchType =
		gtk_label_new
		(
			"New type"
		);
	gtk_widget_show(label_renameBoundary_selectedPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					label_renameBoundary_selectedPatchType, 0, 1, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- add the entry for selected patch type
	GtkWidget* entry_renameBoundary_selectedPatchType = gtk_entry_new();
	gtk_widget_show(entry_renameBoundary_selectedPatchType);
	
	gtk_table_attach(GTK_TABLE(table1),
					entry_renameBoundary_selectedPatchType, 1, 2, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//- set the button for removing new patch names
	GtkWidget* button2 = gtk_button_new();
	gtk_widget_show(button2);
	gtk_table_attach
	(
		GTK_TABLE(table1),
		button2,
		2, 3, 7, 10,
		(GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 5
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

	GLADE_HOOKUP_OBJECT_NO_REF
	(
		renameBoundaryFramePtr_,
		renameBoundaryFramePtr_,
		"renameBoundaryFramePtr_"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		checkbutton_renameBoundary_defaultPatchName,
		"checkbutton_renameBoundary_defaultPatchName"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_defaultPatchName,
		"entry_renameBoundary_defaultPatchName"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		checkbutton_renameBoundary_defaultPatchType,
		"checkbutton_renameBoundary_defaultPatchType"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_defaultPatchType,
		"entry_renameBoundary_defaultPatchType"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		comboboxentry_renameBoundary_availablePatches,
		"comboboxentry_renameBoundary_availablePatches"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_newPatchName,
		"entry_renameBoundary_newPatchName"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		checkbutton_renameBoundary_newPatchType,
		"checkbutton_renameBoundary_newPatchType"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_newPatchType,
		"entry_renameBoundary_newPatchType"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		comboboxentry_renameBoundary_selectedPatches,
		"comboboxentry_renameBoundary_selectedPatches"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_selectedPatchName,
		"entry_renameBoundary_selectedPatchName"
	);
	GLADE_HOOKUP_OBJECT
	(
		renameBoundaryFramePtr_,
		entry_renameBoundary_selectedPatchType,
		"entry_renameBoundary_selectedPatchType"
	);
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	setDefaultValues(renameBoundaryFramePtr_);
	updateNewPatchNames(renameBoundaryFramePtr_);
	
	g_signal_connect
	(
		G_OBJECT
		(
			GTK_ENTRY
			(
				GTK_COMBO(comboboxentry_renameBoundary_selectedPatches)->entry
			)
		),
		"changed",
		G_CALLBACK(updateNewPatchNamesForSelectedPatch),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(checkbutton_renameBoundary_defaultPatchName),
		"clicked",
		G_CALLBACK(useDefaultPatchName),
		renameBoundaryFramePtr_
	);

	g_signal_connect
	(
		G_OBJECT(entry_renameBoundary_defaultPatchName),
		"changed",
		G_CALLBACK(setDefaultPatchName),
		renameBoundaryFramePtr_
	);

	g_signal_connect
	(
		G_OBJECT(checkbutton_renameBoundary_defaultPatchType),
		"clicked",
		G_CALLBACK(useDefaultPatchType),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_renameBoundary_defaultPatchType),
		"changed",
		G_CALLBACK(setDefaultPatchType),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(checkbutton_renameBoundary_newPatchType),
		"clicked",
		G_CALLBACK(useNewPatchType),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_renameBoundary_selectedPatchName),
		"changed",
		G_CALLBACK(resetNewPatchNamesForSelectedPatch),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(entry_renameBoundary_selectedPatchType),
		"changed",
		G_CALLBACK(resetNewPatchNamesForSelectedPatch),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(button1),
		"clicked",
		G_CALLBACK(addNewPatchName),
		renameBoundaryFramePtr_
	);
	
	g_signal_connect
	(
		G_OBJECT(button2),
		"clicked",
		G_CALLBACK(removeNewPatchName),
		renameBoundaryFramePtr_
	);
	
	gtk_widget_show(renameBoundaryFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
