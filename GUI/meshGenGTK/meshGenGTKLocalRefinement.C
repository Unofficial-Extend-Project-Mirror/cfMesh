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
#include "patchRefinementList.H"
#include "HashSet.H"
#include "triSurface.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
extern "C"
{
    
static void updateSizeForSelectedPatch
(
    GtkWidget* widget,
    GtkWidget* entry_localRefinement_selectedCellSize
)
{
    //- display the size for a patch selected in combo box
    const word name = 
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                widget
            )
        );
    
    if( name == "" )
    {
        gtk_entry_set_text
        (
            GTK_ENTRY(entry_localRefinement_selectedCellSize),
            ""
        );
    }

    patchRefinementList refPatches = guiPtr->patchCellSize();
    forAll(refPatches, patchI)
        if( refPatches[patchI].patchName() == name )
        {
            const word s = help::scalarToText(refPatches[patchI].cellSize());
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_localRefinement_selectedCellSize),
                s.c_str()
            );
        }
}

static void resetSizeForSelectedPatch
(
    GtkWidget* widget,
    GtkWidget* comboboxentry_localRefinement_selectedPatches
)
{
    const word name =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_selectedPatches
                )->entry
            )
        );
    
    if( name == "" )
        return;
    
    guiPtr->removePatchCellSize(name);
    
    const word text =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                widget
            )
        );
    const scalar size = help::textToScalar(text);

    patchRefinement pr(name, Foam::min(size, guiPtr->maxCellSize()));

    guiPtr->addPatchCellSize(pr);
}
    
static void updateLocalRefinementPatches(GtkWidget* localRefinementFramePtr)
{
    GtkWidget* comboboxentry_localRefinement_availablePatches =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "comboboxentry_localRefinement_availablePatches"
        );
    GtkWidget* comboboxentry_localRefinement_selectedPatches =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "comboboxentry_localRefinement_selectedPatches"
        );
    //- update available patches
    patchRefinementList alreadyUsed = guiPtr->patchCellSize();

    wordHashSet usedNames;
    forAll(alreadyUsed, nameI)
        usedNames.insert(alreadyUsed[nameI].patchName());

    GList* avPatches = 0;
    label nAvailable(0);
    const triSurface& surf = guiPtr->surface();
    forAll(surf.patches(), patchI)
    {
        if( usedNames.found(surf.patches()[patchI].name()) )
            continue;
        const char* name = surf.patches()[patchI].name().c_str();
        avPatches = g_list_append(avPatches, const_cast<char*>(name));
        ++nAvailable;
    }

    gtk_combo_set_popdown_strings
    (
        GTK_COMBO
        (
            comboboxentry_localRefinement_availablePatches
        ), avPatches
    );
    
    if( nAvailable == 0 )
    {
        gtk_entry_set_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_availablePatches
                )->entry
            ), ""
        );
    }
    
    //- update used patches
    GList* selPatches = 0;
    forAll(alreadyUsed, patchI)
    {
        const char* name = alreadyUsed[patchI].patchName().c_str();
        selPatches = g_list_append(selPatches, const_cast<char*>(name));
    }
    
    gtk_combo_set_popdown_strings
    (
        GTK_COMBO
        (
            comboboxentry_localRefinement_selectedPatches
        ), selPatches
    );
    
    if( alreadyUsed.size() == 0 )
    {
        gtk_entry_set_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_selectedPatches
                )->entry
            ), ""
        );
    }

    g_list_free(avPatches);
    g_list_free(selPatches);
    
    GtkWidget* entry_localRefinement_selectedCellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "entry_localRefinement_selectedCellSize"
        );
    
    updateSizeForSelectedPatch
    (
        GTK_COMBO
        (
            comboboxentry_localRefinement_selectedPatches
        )->entry,
        entry_localRefinement_selectedCellSize
    );
}
    
static void addPatchCellSize
(
    GtkWidget* widget,
    GtkWidget* localRefinementFramePtr
)
{
    GtkWidget* comboboxentry_localRefinement_availablePatches =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "comboboxentry_localRefinement_availablePatches"
        );
    GtkWidget* entry_localRefinement_cellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "entry_localRefinement_cellSize"
        );
    
    const word name =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_availablePatches
                )->entry
            )
        );
    
    if( name == "" )
        return;
    
    const word text =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                entry_localRefinement_cellSize
            )
        );
    const scalar size = help::textToScalar(text);

    patchRefinement pr(name, Foam::min(size, guiPtr->maxCellSize()));

    guiPtr->addPatchCellSize(pr);

    updateLocalRefinementPatches(localRefinementFramePtr);
}

static void removePatchRefinement
(
    GtkWidget* widget,
    GtkWidget* localRefinementFramePtr
)
{
    GtkWidget* comboboxentry_localRefinement_selectedPatches =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(localRefinementFramePtr),
            "comboboxentry_localRefinement_selectedPatches"
        );
    const word name =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_selectedPatches
                )->entry
            )
        );
    
    if( name == "" )
        return;
    
    guiPtr->removePatchCellSize(name);
    
    updateLocalRefinementPatches(localRefinementFramePtr);
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGTK::createLocalRefinementWindowPage()
{
    guiPtr = &meshGui_;
    
    localRefinementFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a table
    GtkWidget* table1 = gtk_table_new(2, 4, FALSE);
    gtk_widget_show(table1);
    gtk_container_add
    (
        GTK_CONTAINER(localRefinementFramePtr_),
        table1
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- set label for available patches
    GtkWidget* label_localRefinement_availablePatches =
        gtk_label_new("Available patches");
    gtk_widget_show(label_localRefinement_availablePatches);
    gtk_table_attach
    (
        GTK_TABLE(table1),
        label_localRefinement_availablePatches,
        0, 1, 0, 1,
        (GtkAttachOptions)(GTK_FILL),
        (GtkAttachOptions) (0), 0, 0
    );
    gtk_misc_set_alignment
    (
        GTK_MISC(label_localRefinement_availablePatches),
        0,
        0.5
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create the combo box for available patches
    GtkWidget* comboboxentry_localRefinement_availablePatches =
        gtk_combo_new();
    gtk_widget_show
    (
        comboboxentry_localRefinement_availablePatches
    );
    gtk_table_attach
    (
        GTK_TABLE(table1),
        comboboxentry_localRefinement_availablePatches,
        1, 2, 0, 1,
        (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
        (GtkAttachOptions)(GTK_FILL), 0, 0
    );
    
    gtk_editable_set_editable
    (
        GTK_EDITABLE
        (
            GTK_COMBO
            (
                comboboxentry_localRefinement_availablePatches
            )->entry
        ),
        0
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for cell size
    GtkWidget* entry_localRefinement_cellSize = gtk_entry_new();
    gtk_widget_show(entry_localRefinement_cellSize);
    
    gtk_table_attach
    (
        GTK_TABLE(table1),
        entry_localRefinement_cellSize,
        2, 3, 0, 1,
        (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
        (GtkAttachOptions)(GTK_FILL), 0, 0
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- set the button for adding patches
    GtkWidget* button1 = gtk_button_new();
    gtk_widget_show(button1);
    gtk_table_attach
    (
        GTK_TABLE(table1),
        button1,
        3, 4, 0, 1,
        (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 0
    );
    
    g_signal_connect
    (
        G_OBJECT(button1),
        "clicked",
        G_CALLBACK(addPatchCellSize),
        localRefinementFramePtr_
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
    //- set label for selected patches
    GtkWidget* label_localRefinement_selectedPatches =
        gtk_label_new("Selected patches");
    gtk_widget_show(label_localRefinement_selectedPatches);
    gtk_table_attach
    (
        GTK_TABLE(table1),
        label_localRefinement_selectedPatches,
        0, 1, 1, 2,
        (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 0
    );
    gtk_misc_set_alignment
    (
        GTK_MISC(label_localRefinement_selectedPatches),
        0,
        0.5
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for selected cell size
    GtkWidget* entry_localRefinement_selectedCellSize = gtk_entry_new();
    gtk_widget_show(entry_localRefinement_selectedCellSize);
    
    gtk_table_attach
    (
        GTK_TABLE(table1),
        entry_localRefinement_selectedCellSize,
        2, 3, 1, 2,
        (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),
        (GtkAttachOptions)(GTK_FILL), 0, 0
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create the combo box for selected patches
    GtkWidget* comboboxentry_localRefinement_selectedPatches = gtk_combo_new();
    gtk_widget_show(comboboxentry_localRefinement_selectedPatches);

    gtk_table_attach
    (
        GTK_TABLE(table1),
        comboboxentry_localRefinement_selectedPatches,
        1, 2, 1, 2,
        (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0
    );
    
    gtk_editable_set_editable
    (
        GTK_EDITABLE
        (
            GTK_COMBO
            (
                comboboxentry_localRefinement_selectedPatches
            )->entry
        ),
        0
    );
    
    //- change the displayed value to match the selected patch
    g_signal_connect
    (
        G_OBJECT
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_localRefinement_selectedPatches
                )->entry
            )
        ),
        "changed",
        G_CALLBACK(updateSizeForSelectedPatch),
        entry_localRefinement_selectedCellSize
    );
    
    //- change the selected size in case the user decides it
    g_signal_connect
    (
        G_OBJECT(entry_localRefinement_selectedCellSize),
        "changed",
        G_CALLBACK(resetSizeForSelectedPatch),
        comboboxentry_localRefinement_selectedPatches
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- set the button for removing patches
    GtkWidget* button2 = gtk_button_new();
    gtk_widget_show(button2);
    gtk_table_attach
    (
        GTK_TABLE(table1),
        button2,
        3, 4, 1, 2,
        (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 0
    );
    
    g_signal_connect
    (
        G_OBJECT(button2),
        "clicked",
        G_CALLBACK(removePatchRefinement),
        localRefinementFramePtr_
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
        localRefinementFramePtr_,
        localRefinementFramePtr_,
        "Local refinement frame"
    );
    GLADE_HOOKUP_OBJECT
    (
        localRefinementFramePtr_,
        comboboxentry_localRefinement_availablePatches,
        "comboboxentry_localRefinement_availablePatches"
    );
    GLADE_HOOKUP_OBJECT
    (
        localRefinementFramePtr_,
        entry_localRefinement_cellSize,
        "entry_localRefinement_cellSize"
    );
    GLADE_HOOKUP_OBJECT
    (
        localRefinementFramePtr_,
        comboboxentry_localRefinement_selectedPatches,
        "comboboxentry_localRefinement_selectedPatches"
    );
    GLADE_HOOKUP_OBJECT
    (
        localRefinementFramePtr_,
        entry_localRefinement_selectedCellSize,
        "entry_localRefinement_selectedCellSize"
    );
    
    updateLocalRefinementPatches(localRefinementFramePtr_);
    
    gtk_widget_show(localRefinementFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
