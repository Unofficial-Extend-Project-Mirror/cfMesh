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
    
// Callback functions go here
void saveMeshDictCallback(void* data)
{
    guiPtr->writeDict();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGTK::createMainWindowPage()
{
    guiPtr = &meshGui_;
    
    //- create the main window
    GtkWidget* mainWindowPtr_ = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW (mainWindowPtr_), "Mesh generation GUI");
    g_signal_connect
    (
        G_OBJECT(mainWindowPtr_),
        "destroy",
        G_CALLBACK(gtk_main_quit),
        NULL
    );
    
    GtkWidget* vbox1 = gtk_vbox_new(FALSE, 0);
    gtk_widget_show(vbox1);
    gtk_container_add(GTK_CONTAINER (mainWindowPtr_), vbox1);
    
    //- create a toolbar (for reading and saving meshDict)
    GtkWidget* toolbar1 = gtk_toolbar_new();
    gtk_widget_show(toolbar1);
    gtk_box_pack_start(GTK_BOX(vbox1), toolbar1, FALSE, FALSE, 0);
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar1), GTK_TOOLBAR_ICONS);
    //GtkIconSize tmp_toolbar_icon_size =
        //gtk_toolbar_get_icon_size(GTK_TOOLBAR(toolbar1));
    
    GtkWidget* toolitem1 = (GtkWidget*)gtk_tool_item_new();
    gtk_widget_show(toolitem1);
    gtk_container_add(GTK_CONTAINER(toolbar1), toolitem1);
    
    GtkWidget* button_mainWindow_save =
        gtk_button_new_with_mnemonic("Save meshDict");
    gtk_widget_show(button_mainWindow_save);
    gtk_container_add(GTK_CONTAINER(toolitem1), button_mainWindow_save);
    
    g_signal_connect
    (
        G_OBJECT(button_mainWindow_save),
        "clicked",
        G_CALLBACK(saveMeshDictCallback),
        static_cast<void*>(&meshGui_)
    );

    //- create exit button
    GtkWidget* toolitem2 = (GtkWidget*)gtk_tool_item_new();
    gtk_widget_show(toolitem2);
    gtk_container_add(GTK_CONTAINER(toolbar1), toolitem2);
    
    GtkWidget* button_mainWindow_exit = gtk_button_new_with_mnemonic("Exit");
    gtk_widget_show(button_mainWindow_exit);
    gtk_container_add(GTK_CONTAINER(toolitem2), button_mainWindow_exit);
    
    g_signal_connect
    (
        G_OBJECT(button_mainWindow_exit),
        "clicked",
        G_CALLBACK(gtk_main_quit),
        NULL
    );
    
    //- create a notebook which allows the user to switch between
    //- the local and global settings
    GtkWidget* notebook1 = gtk_notebook_new();
    
    gtk_box_pack_start(GTK_BOX(vbox1), notebook1, TRUE, TRUE, 0);
    
    //- append the page for general settings
    GtkWidget* label_generalSettings = gtk_label_new("General settings");
    gtk_widget_show(label_generalSettings);

    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        generalPageFramePtr_,
        label_generalSettings
    );
    
    //- add local settings    
    GtkWidget* label_localSettings = gtk_label_new("Local settings");
    gtk_widget_show(label_localSettings);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        localSettingsMainFramePtr_,
        label_localSettings
    );
    
    //- add object refinements
    GtkWidget* label_objectRefinements = gtk_label_new("Object refinements");
    gtk_widget_show(label_objectRefinements);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        objectRefinementMainFramePtr_,
        label_objectRefinements
    );
    
    //- add patch renaming
    GtkWidget* label_renameBoundary = gtk_label_new("Rename boundary patches");
    gtk_widget_show(label_renameBoundary);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        renameBoundaryFramePtr_,
        label_renameBoundary
    );
    
    gtk_widget_show(notebook1);
    
    /* Store pointers to all widgets, for use by lookup_widget(). */
    GLADE_HOOKUP_OBJECT_NO_REF
    (
        mainWindowPtr_,
        mainWindowPtr_,
        "Mesh generation GUI"
    );

    GLADE_HOOKUP_OBJECT
    (
        mainWindowPtr_,
        button_mainWindow_save,
        "button_mainWindow_save"
    );
    GLADE_HOOKUP_OBJECT
    (
        mainWindowPtr_,
        button_mainWindow_exit,
        "button_mainWindow_exit"
    );
    
    gtk_widget_show(mainWindowPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
