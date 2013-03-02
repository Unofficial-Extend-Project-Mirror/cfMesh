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
#include "sphereRefinement.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- callback functions

extern "C"
{

static void showDataForSphereRefinement
(
    GtkWidget* widget,
    GtkWidget* sphereRefinementFramePtr
)
{
    GtkWidget* comboboxentry_sphereRefinement_createdSpheres =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "comboboxentry_sphereRefinement_createdSpheres"
        );
    GtkWidget* entry_sphereRefinement_selectedCellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCellSize"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreX =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreX"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreY =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreY"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreZ =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreZ"
        );
    GtkWidget* entry_sphereRefinement_selectedRadius =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedRadius"
        );
    
    //- find the object in the dictionary
    const word sphereName =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_sphereRefinement_createdSpheres)->entry
            )
        );
    
    const PtrList<entry> refs = guiPtr->objectRefinements();
    
    forAll(refs, refI)
        if( refs[refI].keyword() == sphereName )
        {
            const dictionary dict = refs[refI].dict();
            
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_sphereRefinement_selectedCellSize),
                help::scalarToText
                (
                    scalar(readScalar(dict.lookup("cellSize")))
                ).c_str()
            );
            const vector centre(dict.lookup("centre"));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_sphereRefinement_selectedCentreX),
                help::scalarToText(centre.x()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_sphereRefinement_selectedCentreY),
                help::scalarToText(centre.y()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_sphereRefinement_selectedCentreZ),
                help::scalarToText(centre.z()).c_str()
            );
            
            const scalar radius(readScalar(dict.lookup("radius")));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_sphereRefinement_selectedRadius),
                help::scalarToText(radius).c_str()
            );
            
            return;
        }
        
    gtk_entry_set_text(GTK_ENTRY(entry_sphereRefinement_selectedCellSize), "");
    gtk_entry_set_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreX), "");
    gtk_entry_set_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreY), "");
    gtk_entry_set_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreZ), "");
    gtk_entry_set_text(GTK_ENTRY(entry_sphereRefinement_selectedRadius), "");
}

static void resetValuesForSelectedLine
(
    GtkWidget* widget,
    GtkWidget* sphereRefinementFramePtr
)
{
    GtkWidget* comboboxentry_sphereRefinement_createdSpheres =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "comboboxentry_sphereRefinement_createdSpheres"
        );
    GtkWidget* entry_sphereRefinement_selectedCellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCellSize"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreX =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreX"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreY =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreY"
        );
    GtkWidget* entry_sphereRefinement_selectedCentreZ =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedCentreZ"
        );
    GtkWidget* entry_sphereRefinement_selectedRadius =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_selectedRadius"
        );
        
    const word name =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_sphereRefinement_createdSpheres
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
    if( widget == entry_sphereRefinement_selectedCellSize )
    {
        const word size =
            gtk_entry_get_text(GTK_ENTRY(entry_sphereRefinement_selectedCellSize));
        dict.remove("cellSize");
        dict.add("cellSize", help::textToScalar(size));
    }
    else if( widget == entry_sphereRefinement_selectedCentreX )
    {
        point centre(dict.lookup("centre"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreX));
        dict.remove("centre");
        centre.x() = help::textToScalar(c);
        dict.add("centre", centre);
    }
    else if( widget == entry_sphereRefinement_selectedCentreY )
    {
        point centre(dict.lookup("centre"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreY));
        dict.remove("centre");
        centre.y() = help::textToScalar(c);
        dict.add("centre", centre);
    }
    else if( widget == entry_sphereRefinement_selectedCentreZ )
    {
        point centre(dict.lookup("centre"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_sphereRefinement_selectedCentreZ));
        dict.remove("centre");
        centre.z() = help::textToScalar(c);
        dict.add("centre", centre);
    }
    else if( widget == entry_sphereRefinement_selectedRadius )
    {
        const word l =
            gtk_entry_get_text(GTK_ENTRY(entry_sphereRefinement_selectedRadius));
        dict.remove("radius");
        scalar radius = help::textToScalar(l);
        dict.add("radius", radius);
    }
    
    sphereRefinement sr(name, dict);
    guiPtr->addObjectRefinement(sr);
}

static void updateSphereRefinements(GtkWidget* sphereRefinementFramePtr)
{
    GtkWidget* comboboxentry_sphereRefinement_createdSpheres =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "comboboxentry_sphereRefinement_createdSpheres"
        );
    
    const PtrList<entry> refs = guiPtr->objectRefinements();
    label nSpheres(0);
    GList* createdSpheres = 0;
    forAll(refs, refI)
    {
        const word& key = refs[refI].keyword();
        const dictionary& dict = refs[refI].dict();
        const word type(dict.lookup("type"));
        
        if( type == "sphere" )
        {
            ++nSpheres;
            sphereRefinement sphere(key, dict);
            createdSpheres =
                g_list_append
                (
                    createdSpheres,
                    const_cast<char*>(sphere.name().c_str())
                );
        }
    }
        
    gtk_combo_set_popdown_strings
    (
        GTK_COMBO
        (
            comboboxentry_sphereRefinement_createdSpheres
        ), createdSpheres
    );
    
    if( nSpheres == 0 )
    {
        gtk_entry_set_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_sphereRefinement_createdSpheres)->entry
            ),
            ""
        );
    }
        
    g_list_free(createdSpheres);
    
    showDataForSphereRefinement(NULL, sphereRefinementFramePtr);
}

static void addSphereRefinement
(
    GtkWidget* button,
    GtkWidget* sphereRefinementFramePtr
)
{
    GtkWidget* name =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereName"
        );
    GtkWidget* cellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_cellSize"
        );
    GtkWidget* centreX =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_centreX"
        );
    GtkWidget* centreY =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_centreY"
        );
    GtkWidget* centreZ =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_centreZ"
        );
    GtkWidget* radius =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "entry_sphereRefinement_radius"
        );
    
    //- create all necessary settings
    const word sphereName = gtk_entry_get_text(GTK_ENTRY(name));
    gtk_entry_set_text(GTK_ENTRY(name), "");
    const scalar sphereCellSize =
        help::textToScalar(gtk_entry_get_text(GTK_ENTRY(cellSize)));
    gtk_entry_set_text(GTK_ENTRY(cellSize), "");
    point centre, p1;
    centre.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreX)));
    gtk_entry_set_text(GTK_ENTRY(centreX), "");
    centre.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreY)));
    gtk_entry_set_text(GTK_ENTRY(centreY), "");
    centre.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(centreZ)));
    gtk_entry_set_text(GTK_ENTRY(centreZ), "");
    const scalar selectedRadius =
        help::textToScalar(gtk_entry_get_text(GTK_ENTRY(radius)));
    gtk_entry_set_text(GTK_ENTRY(radius), "");
        
    sphereRefinement sphere
    (
        sphereName,
        sphereCellSize,
        centre,
        selectedRadius
    );
    
    guiPtr->addObjectRefinement(sphere);
    
    updateSphereRefinements(sphereRefinementFramePtr);
}

static void removeSphereRefinement
(
    GtkWidget* button,
    GtkWidget* sphereRefinementFramePtr
)
{
    GtkWidget* createdSpheres =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(sphereRefinementFramePtr),
            "comboboxentry_sphereRefinement_createdSpheres"
        );
    
    const word sphereName =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(createdSpheres)->entry
            )
        );

    guiPtr->removeObjectRefinement(sphereName);

    updateSphereRefinements(sphereRefinementFramePtr);
}

}

void meshGenGTK::createSphereRefinementWindowPage()
{
    guiPtr = &meshGui_;
    
    sphereRefinementFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a table
    GtkWidget* table1 = gtk_table_new(5, 11, FALSE);
    gtk_widget_show(table1);
    gtk_container_add
    (
        GTK_CONTAINER(sphereRefinementFramePtr_),
        table1
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for sphere name
    GtkWidget* label_sphereName = gtk_label_new("Sphere name");
    gtk_widget_show(label_sphereName);
    gtk_table_attach(GTK_TABLE(table1), label_sphereName, 0, 1, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_sphereName), 0, 0.5);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create an entry for box name
    GtkWidget* entry_sphereName = gtk_entry_new();
    gtk_widget_show(entry_sphereName);
    
    gtk_table_attach(GTK_TABLE(table1), entry_sphereName, 1, 4, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for desired cell size
    GtkWidget* label_sphereRefinement_cellSize = gtk_label_new("Cell size [m]");
    gtk_widget_show(label_sphereRefinement_cellSize);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_cellSize,
                    0, 1, 1, 2,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_sphereRefinement_cellSize), 0, 0.5);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for desired cell size
    GtkWidget* entry_sphereRefinement_cellSize = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_cellSize);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_cellSize,
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
    //- label for centre
    GtkWidget* label_sphereRefinement_centre = gtk_label_new("Centre");
    gtk_widget_show(label_sphereRefinement_centre);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_centre,
                    0, 1, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_sphereRefinement_centre), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for centre
    GtkWidget* entry_sphereRefinement_centreX = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_centreX);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_centreX,
                    1, 2, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_sphereRefinement_centreY = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_centreY);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_centreY,
                    2, 3, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_sphereRefinement_centreZ = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_centreZ);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_centreZ,
                    3, 4, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for radius
    GtkWidget* label_sphereRefinement_p1 = gtk_label_new("Radius [m]");
    gtk_widget_show(label_sphereRefinement_p1);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_p1,
                    0, 1, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_sphereRefinement_p1), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for radius
    GtkWidget* entry_sphereRefinement_radius = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_radius);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_radius,
                    1, 4, 4, 5,
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
        G_CALLBACK(addSphereRefinement),
        sphereRefinementFramePtr_
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
    GtkWidget* label_sphereRefinement_existing =
        gtk_label_new("Existing spheres");
    gtk_widget_show(label_sphereRefinement_existing);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_existing,
                    0, 1, 6, 7,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- combobox for existing box objects
    GtkWidget* comboboxentry_sphereRefinement_createdSpheres =
        gtk_combo_new();
    gtk_widget_show(comboboxentry_sphereRefinement_createdSpheres);
    gtk_table_attach
    (
        GTK_TABLE(table1),comboboxentry_sphereRefinement_createdSpheres,
        1, 2, 6, 7,
        (GtkAttachOptions)(GTK_FILL),
        (GtkAttachOptions)(0), 0, 0
    );
                    
    gtk_editable_set_editable
    (
        GTK_EDITABLE
        (
            GTK_COMBO(comboboxentry_sphereRefinement_createdSpheres)->entry
        ),
        0
    );
    
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for selected cell size
    GtkWidget* label_sphereRefinement_selectedCellSize =
        gtk_label_new("Cell size [m]");
    gtk_widget_show(label_sphereRefinement_selectedCellSize);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_selectedCellSize,
                    0, 1, 7, 8,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_sphereRefinement_selectedCellSize), 0, 0.5
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for selected cell size
    GtkWidget* entry_sphereRefinement_selectedCellSize = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_selectedCellSize);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_selectedCellSize,
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
    //- label for centre
    GtkWidget* label_sphereRefinement_selectedCentre =
        gtk_label_new("Centre");
    gtk_widget_show(label_sphereRefinement_selectedCentre);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_selectedCentre,
                    0, 1, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_sphereRefinement_selectedCentre), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for selected centre
    GtkWidget* entry_sphereRefinement_selectedCentreX = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_selectedCentreX);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_selectedCentreX,
                    1, 2, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_sphereRefinement_selectedCentreY = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_selectedCentreY);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_selectedCentreY,
                    2, 3, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_sphereRefinement_selectedCentreZ = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_selectedCentreZ);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_selectedCentreZ,
                    3, 4, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for selected radius
    GtkWidget* label_sphereRefinement_selectedRadius =
        gtk_label_new("Radius [m]");
    gtk_widget_show(label_sphereRefinement_selectedRadius);
    gtk_table_attach(GTK_TABLE(table1), label_sphereRefinement_selectedRadius,
                    0, 1, 10, 11,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_sphereRefinement_selectedRadius), 0, 0.5
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for radius
    GtkWidget* entry_sphereRefinement_selectedRadius = gtk_entry_new();
    gtk_widget_show(entry_sphereRefinement_selectedRadius);
    gtk_table_attach(GTK_TABLE(table1), entry_sphereRefinement_selectedRadius,
                    1, 4, 10, 11,
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
        G_CALLBACK(removeSphereRefinement),
        sphereRefinementFramePtr_
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
        sphereRefinementFramePtr_,
        sphereRefinementFramePtr_,
        "sphereRefinementFramePtr_");
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        table1,
        "table1"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereName,
        "entry_sphereName"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_cellSize,
        "entry_sphereRefinement_cellSize"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_centreX,
        "entry_sphereRefinement_centreX"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_centreY,
        "entry_sphereRefinement_centreY"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_centreZ,
        "entry_sphereRefinement_centreZ"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_radius,
        "entry_sphereRefinement_radius"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        comboboxentry_sphereRefinement_createdSpheres,
        "comboboxentry_sphereRefinement_createdSpheres"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_selectedCellSize,
        "entry_sphereRefinement_selectedCellSize"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_selectedCentreX,
        "entry_sphereRefinement_selectedCentreX"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_selectedCentreY,
        "entry_sphereRefinement_selectedCentreY"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_selectedCentreZ,
        "entry_sphereRefinement_selectedCentreZ"
    );
    GLADE_HOOKUP_OBJECT
    (
        sphereRefinementFramePtr_,
        entry_sphereRefinement_selectedRadius,
        "entry_sphereRefinement_selectedRadius"
    );
    
    updateSphereRefinements(sphereRefinementFramePtr_);
    
    g_signal_connect
    (
        G_OBJECT
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_sphereRefinement_createdSpheres)->entry
            )
        ),
        "changed",
        G_CALLBACK(showDataForSphereRefinement),
        sphereRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_sphereRefinement_selectedCellSize),
        "changed",
        G_CALLBACK(resetValuesForSelectedLine),
        sphereRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_sphereRefinement_selectedCentreX),
        "changed",
        G_CALLBACK(resetValuesForSelectedLine),
        sphereRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_sphereRefinement_selectedCentreY),
        "changed",
        G_CALLBACK(resetValuesForSelectedLine),
        sphereRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_sphereRefinement_selectedCentreZ),
        "changed",
        G_CALLBACK(resetValuesForSelectedLine),
        sphereRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_sphereRefinement_selectedRadius),
        "changed",
        G_CALLBACK(resetValuesForSelectedLine),
        sphereRefinementFramePtr_
    );
    
    gtk_widget_show(sphereRefinementFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
