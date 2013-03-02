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
#include "coneRefinement.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- callback functions

extern "C"
{

static void showDataForConeRefinement
(
    GtkWidget* widget,
    GtkWidget* coneRefinementFramePtr
)
{
    GtkWidget* comboboxentry_coneRefinement_createdCones =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "comboboxentry_coneRefinement_createdCones"
        );
    GtkWidget* entry_coneRefinement_selectedCellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedCellSize"
        );
    GtkWidget* entry_coneRefinement_selectedP0X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0X"
        );
    GtkWidget* entry_coneRefinement_selectedP0Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0Y"
        );
    GtkWidget* entry_coneRefinement_selectedP0Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0Z"
        );
    GtkWidget* entry_coneRefinement_selectedRadius0 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedRadius0"
        );
    GtkWidget* entry_coneRefinement_selectedP1X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1X"
        );
    GtkWidget* entry_coneRefinement_selectedP1Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1Y"
        );
    GtkWidget* entry_coneRefinement_selectedP1Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1Z"
        );
    GtkWidget* entry_coneRefinement_selectedRadius1 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedRadius1"
        );
    
    //- find the object in the dictionary
    const word coneName =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_coneRefinement_createdCones)->entry
            )
        );
    
    const PtrList<entry> refs = guiPtr->objectRefinements();
    
    forAll(refs, refI)
        if( refs[refI].keyword() == coneName )
        {
            const dictionary dict = refs[refI].dict();
            
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedCellSize),
                help::scalarToText
                (
                    scalar(readScalar(dict.lookup("cellSize")))
                ).c_str()
            );
            const vector p0(dict.lookup("p0"));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP0X),
                help::scalarToText(p0.x()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP0Y),
                help::scalarToText(p0.y()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP0Z),
                help::scalarToText(p0.z()).c_str()
            );
            const scalar r0(readScalar(dict.lookup("radius0")));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedRadius0),
                help::scalarToText(r0).c_str()
            );
            
            const vector p1(dict.lookup("p1"));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP1X),
                help::scalarToText(p1.x()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP1Y),
                help::scalarToText(p1.y()).c_str()
            );
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedP1Z),
                help::scalarToText(p1.z()).c_str()
            );
            const scalar r1(readScalar(dict.lookup("radius1")));
            gtk_entry_set_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedRadius1),
                help::scalarToText(r1).c_str()
            );
            
            return;
        }
        
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedCellSize), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP0X), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP0Y), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP0Z), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedRadius0), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP1X), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP1Y), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedP1Z), "");
    gtk_entry_set_text(GTK_ENTRY(entry_coneRefinement_selectedRadius1), "");
}

static void resetValuesForSelectedCone
(
    GtkWidget* widget,
    GtkWidget* coneRefinementFramePtr
)
{
    GtkWidget* comboboxentry_coneRefinement_createdCones =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "comboboxentry_coneRefinement_createdCones"
        );
    GtkWidget* entry_coneRefinement_selectedCellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedCellSize"
        );
    GtkWidget* entry_coneRefinement_selectedP0X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0X"
        );
    GtkWidget* entry_coneRefinement_selectedP0Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0Y"
        );
    GtkWidget* entry_coneRefinement_selectedP0Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP0Z"
        );
    GtkWidget* entry_coneRefinement_selectedRadius0 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedRadius0"
        );
    GtkWidget* entry_coneRefinement_selectedP1X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1X"
        );
    GtkWidget* entry_coneRefinement_selectedP1Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1Y"
        );
    GtkWidget* entry_coneRefinement_selectedP1Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedP1Z"
        );
    GtkWidget* entry_coneRefinement_selectedRadius1 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_selectedRadius1"
        );
        
    const word name =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO
                (
                    comboboxentry_coneRefinement_createdCones
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
    if( widget == entry_coneRefinement_selectedCellSize )
    {
        const word size =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedCellSize));
        dict.remove("cellSize");
        dict.add("cellSize", help::textToScalar(size));
    }
    else if( widget == entry_coneRefinement_selectedP0X )
    {
        point p0(dict.lookup("p0"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP0X));
        dict.remove("p0");
        p0.x() = help::textToScalar(c);
        dict.add("p0", p0);
    }
    else if( widget == entry_coneRefinement_selectedP0Y )
    {
        point p0(dict.lookup("p0"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP0Y));
        dict.remove("p0");
        p0.y() = help::textToScalar(c);
        dict.add("p0", p0);
    }
    else if( widget == entry_coneRefinement_selectedP0Z )
    {
        point p0(dict.lookup("p0"));
        const word c =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP0Z));
        dict.remove("p0");
        p0.z() = help::textToScalar(c);
        dict.add("p0", p0);
    }
    else if( widget == entry_coneRefinement_selectedRadius0 )
    {
        scalar r0(readScalar(dict.lookup("radius0")));
        const word c =
            gtk_entry_get_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedRadius0)
            );
        dict.remove("radius0");
        r0 = help::textToScalar(c);
        dict.add("radius0", r0);
    }
    else if( widget == entry_coneRefinement_selectedP1X )
    {
        const word l =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP1X));
        point p1(dict.lookup("p1"));
        dict.remove("p1");
        p1.x() = help::textToScalar(l);
        dict.add("p1", p1);
    }
    else if( widget == entry_coneRefinement_selectedP1Y )
    {
        const word l =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP1Y));
        point p1(dict.lookup("p1"));
        dict.remove("p1");
        p1.y() = help::textToScalar(l);
        dict.add("p1", p1);
    }
    else if( widget == entry_coneRefinement_selectedP1Z )
    {
        const word l =
            gtk_entry_get_text(GTK_ENTRY(entry_coneRefinement_selectedP1Z));
        point p1(dict.lookup("p1"));
        dict.remove("p1");
        p1.z() = help::textToScalar(l);
        dict.add("p1", p1);
    }
    else if( widget == entry_coneRefinement_selectedRadius1 )
    {
        scalar r1(readScalar(dict.lookup("radius1")));
        const word c =
            gtk_entry_get_text
            (
                GTK_ENTRY(entry_coneRefinement_selectedRadius1)
            );
        dict.remove("radius1");
        r1 = help::textToScalar(c);
        dict.add("radius1", r1);
    }
    
    coneRefinement cr(name, dict);
    guiPtr->addObjectRefinement(cr);
}

static void updateConeRefinements(GtkWidget* coneRefinementFramePtr)
{
    GtkWidget* comboboxentry_coneRefinement_createdCones =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "comboboxentry_coneRefinement_createdCones"
        );
    
    const PtrList<entry> refs = guiPtr->objectRefinements();
    label nCones(0);
    GList* createdCones = 0;
    forAll(refs, refI)
    {
        const word& key = refs[refI].keyword();
        const dictionary& dict = refs[refI].dict();
        const word type(dict.lookup("type"));
        
        if( type == "cone" )
        {
            ++nCones;
            coneRefinement cone(key, dict);
            createdCones =
                g_list_append
                (
                    createdCones,
                    const_cast<char*>(cone.name().c_str())
                );
        }
    }
        
    gtk_combo_set_popdown_strings
    (
        GTK_COMBO
        (
            comboboxentry_coneRefinement_createdCones
        ), createdCones
    );
    
    if( nCones == 0 )
    {
        gtk_entry_set_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_coneRefinement_createdCones)->entry
            ),
            ""
        );
    }
        
    g_list_free(createdCones);
    
    showDataForConeRefinement(NULL, coneRefinementFramePtr);
}

static void addConeRefinement
(
    GtkWidget* button,
    GtkWidget* coneRefinementFramePtr
)
{
    GtkWidget* name =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneName"
        );
    GtkWidget* cellSize =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_cellSize"
        );
    GtkWidget* p0X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p0X"
        );
    GtkWidget* p0Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p0Y"
        );
    GtkWidget* p0Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p0Z"
        );
    GtkWidget* r0 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_radius0"
        );
    GtkWidget* p1X =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p1X"
        );
    GtkWidget* p1Y =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p1Y"
        );
    GtkWidget* p1Z =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_p1Z"
        );
    GtkWidget* r1 =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "entry_coneRefinement_radius1"
        );
    
    //- create all necessary settings
    const word coneName = gtk_entry_get_text(GTK_ENTRY(name));
    gtk_entry_set_text(GTK_ENTRY(name), "");
    const scalar coneCellSize =
        help::textToScalar(gtk_entry_get_text(GTK_ENTRY(cellSize)));
    gtk_entry_set_text(GTK_ENTRY(cellSize), "");
    point p0, p1;
    scalar radius0, radius1;
    p0.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0X)));
    gtk_entry_set_text(GTK_ENTRY(p0X), "");
    p0.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0Y)));
    gtk_entry_set_text(GTK_ENTRY(p0Y), "");
    p0.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p0Z)));
    gtk_entry_set_text(GTK_ENTRY(p0Z), "");
    radius0 = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(r0)));
    gtk_entry_set_text(GTK_ENTRY(r0), "");
    p1.x() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1X)));
    gtk_entry_set_text(GTK_ENTRY(p1X), "");
    p1.y() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1Y)));
    gtk_entry_set_text(GTK_ENTRY(p1Y), "");
    p1.z() = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(p1Z)));
    gtk_entry_set_text(GTK_ENTRY(p1Z), "");
    radius1 = help::textToScalar(gtk_entry_get_text(GTK_ENTRY(r1)));
    gtk_entry_set_text(GTK_ENTRY(r1), "");
        
    coneRefinement cone
    (
        coneName,
        coneCellSize,
        p0,
        radius0,
        p1,
        radius1
    );
    
    guiPtr->addObjectRefinement(cone);
    
    updateConeRefinements(coneRefinementFramePtr);
}

static void removeConeRefinement
(
    GtkWidget* button,
    GtkWidget* coneRefinementFramePtr
)
{
    GtkWidget* createdCones =
        (GtkWidget*)g_object_get_data
        (
            G_OBJECT(coneRefinementFramePtr),
            "comboboxentry_coneRefinement_createdCones"
        );
    
    const word coneName =
        gtk_entry_get_text
        (
            GTK_ENTRY
            (
                GTK_COMBO(createdCones)->entry
            )
        );

    guiPtr->removeObjectRefinement(coneName);

    updateConeRefinements(coneRefinementFramePtr);
}

}

void meshGenGTK::createConeRefinementWindowPage()
{
    guiPtr = &meshGui_;
    
    coneRefinementFramePtr_ = gtk_frame_new(NULL);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a table
    GtkWidget* table1 = gtk_table_new(5, 15, FALSE);
    gtk_widget_show(table1);
    gtk_container_add
    (
        GTK_CONTAINER(coneRefinementFramePtr_),
        table1
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for cone name
    GtkWidget* label_coneName = gtk_label_new("Cone name");
    gtk_widget_show(label_coneName);
    gtk_table_attach(GTK_TABLE(table1), label_coneName, 0, 1, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneName), 0, 0.5);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create an entry for cone name
    GtkWidget* entry_coneName = gtk_entry_new();
    gtk_widget_show(entry_coneName);
    
    gtk_table_attach(GTK_TABLE(table1), entry_coneName, 1, 4, 0, 1,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for desired cell size
    GtkWidget* label_coneRefinement_cellSize = gtk_label_new("Cell size [m]");
    gtk_widget_show(label_coneRefinement_cellSize);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_cellSize,
                    0, 1, 1, 2,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_cellSize), 0, 0.5);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for desired cell size
    GtkWidget* entry_coneRefinement_cellSize = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_cellSize);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_cellSize,
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
    GtkWidget* label_coneRefinement_p0 = gtk_label_new("Starting point");
    gtk_widget_show(label_coneRefinement_p0);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_p0,
                    0, 1, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_p0), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for p0
    GtkWidget* entry_coneRefinement_p0X = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p0X);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p0X,
                    1, 2, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_p0Y = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p0Y);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p0Y,
                    2, 3, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_p0Z = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p0Z);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p0Z,
                    3, 4, 3, 4,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for p1
    GtkWidget* label_coneRefinement_p1 = gtk_label_new("End point");
    gtk_widget_show(label_coneRefinement_p1);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_p1,
                    0, 1, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_p1), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for p1
    GtkWidget* entry_coneRefinement_p1X = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p1X);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p1X,
                    1, 2, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_p1Y = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p1Y);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p1Y,
                    2, 3, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_p1Z = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_p1Z);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_p1Z,
                    3, 4, 4, 5,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for radius0
    GtkWidget* label_coneRefinement_radius0 =
        gtk_label_new("Starting radius [m]");
    gtk_widget_show(label_coneRefinement_radius0);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_radius0,
                    0, 1, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_radius0), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for radius0
    GtkWidget* entry_coneRefinement_radius0 = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_radius0);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_radius0,
                    1, 4, 5, 6,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for radius1
    GtkWidget* label_coneRefinement_radius1 = gtk_label_new("End radius [m]");
    gtk_widget_show(label_coneRefinement_radius1);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_radius1,
                    0, 1, 6, 7,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_radius1), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for radius1
    GtkWidget* entry_coneRefinement_radius1 = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_radius1);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_radius1,
                    1, 4, 6, 7,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- set the button for adding cones
    GtkWidget* button1 = gtk_button_new();
    gtk_widget_show(button1);
    gtk_table_attach
    (
        GTK_TABLE(table1),
        button1,
        4, 5, 0, 7,
        (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 5
    );
    
    g_signal_connect
    (
        G_OBJECT(button1),
        "clicked",
        G_CALLBACK(addConeRefinement),
        coneRefinementFramePtr_
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
    gtk_table_attach(GTK_TABLE(table1), hseparator0, 0, 5, 7, 8,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(GTK_EXPAND | GTK_SHRINK | GTK_FILL),
                    0, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for existing cone objects
    GtkWidget* label_coneRefinement_existing = gtk_label_new("Existing cones");
    gtk_widget_show(label_coneRefinement_existing);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_existing,
                    0, 1, 8, 9,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- combobox for existing cone objects
    GtkWidget* comboboxentry_coneRefinement_createdCones =
        gtk_combo_new();
    gtk_widget_show(comboboxentry_coneRefinement_createdCones);
    gtk_table_attach
    (
        GTK_TABLE(table1),comboboxentry_coneRefinement_createdCones,
        1, 2, 8, 9,
        (GtkAttachOptions)(GTK_FILL),
        (GtkAttachOptions)(0), 0, 0
    );
                    
    gtk_editable_set_editable
    (
        GTK_EDITABLE
        (
            GTK_COMBO(comboboxentry_coneRefinement_createdCones)->entry
        ),
        0
    );
    
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for selected cell size
    GtkWidget* label_coneRefinement_selectedCellSize =
        gtk_label_new("Cell size [m]");
    gtk_widget_show(label_coneRefinement_selectedCellSize);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_selectedCellSize,
                    0, 1, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_coneRefinement_selectedCellSize), 0, 0.5
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create a label for selected cell size
    GtkWidget* entry_coneRefinement_selectedCellSize = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedCellSize);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedCellSize,
                    1, 4, 9, 10,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- create labels for axes
    label_axes = gtk_label_new("X");
    gtk_widget_show(label_axes);
    gtk_table_attach(GTK_TABLE(table1), label_axes,
                    1, 2, 10, 11,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
    
    label_axes = gtk_label_new("Y");
    gtk_widget_show(label_axes);
    gtk_table_attach(GTK_TABLE(table1), label_axes,
                    2, 3, 10, 11,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
    
    label_axes = gtk_label_new("Z");
    gtk_widget_show(label_axes);
    gtk_table_attach(GTK_TABLE(table1), label_axes,
                    3, 4, 10, 11,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_axes), 0, 0.5);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for p0
    GtkWidget* label_coneRefinement_selectedP0 =
        gtk_label_new("Starting point");
    gtk_widget_show(label_coneRefinement_selectedP0);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_selectedP0,
                    0, 1, 11, 12,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment(GTK_MISC(label_coneRefinement_selectedP0), 0, 0.5);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for selected p0
    GtkWidget* entry_coneRefinement_selectedP0X = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP0X);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP0X,
                    1, 2, 11, 12,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_selectedP0Y = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP0Y);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP0Y,
                    2, 3, 11, 12,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_selectedP0Z = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP0Z);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP0Z,
                    3, 4, 11, 12,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for selected p1
    GtkWidget* label_coneRefinement_selectedP1 =
        gtk_label_new("End point");
    gtk_widget_show(label_coneRefinement_selectedP1);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_selectedP1,
                    0, 1, 12, 13,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_coneRefinement_selectedP1), 0, 0.5
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entries for p1
    GtkWidget* entry_coneRefinement_selectedP1X = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP1X);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP1X,
                    1, 2, 12, 13,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_selectedP1Y = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP1Y);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP1Y,
                    2, 3, 12, 13,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    
    GtkWidget* entry_coneRefinement_selectedP1Z = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedP1Z);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedP1Z,
                    3, 4, 12, 13,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for selected radius0
    GtkWidget* label_coneRefinement_selectedRadius0 =
        gtk_label_new("Starting radius [m]");
    gtk_widget_show(label_coneRefinement_selectedRadius0);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_selectedRadius0,
                    0, 1, 13, 14,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_coneRefinement_selectedRadius0), 0, 0.5
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for selected radius0
    GtkWidget* entry_coneRefinement_selectedRadius0 = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedRadius0);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedRadius0,
                    1, 4, 13, 14,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
                    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- label for radius1
    GtkWidget* label_coneRefinement_selectedRadius1 =
        gtk_label_new("End radius [m]");
    gtk_widget_show(label_coneRefinement_selectedRadius1);
    gtk_table_attach(GTK_TABLE(table1), label_coneRefinement_selectedRadius1,
                    0, 1, 14, 15,
                    (GtkAttachOptions)(GTK_FILL),
                    (GtkAttachOptions)(0), 0, 0);
    gtk_misc_set_alignment
    (
        GTK_MISC(label_coneRefinement_selectedRadius1), 0, 0.5
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- entry for radius1
    GtkWidget* entry_coneRefinement_selectedRadius1 = gtk_entry_new();
    gtk_widget_show(entry_coneRefinement_selectedRadius1);
    gtk_table_attach(GTK_TABLE(table1), entry_coneRefinement_selectedRadius1,
                    1, 4, 14, 15,
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
        4, 5, 8, 15,
        (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 5
    );
    
    g_signal_connect
    (
        G_OBJECT(button2),
        "clicked",
        G_CALLBACK(removeConeRefinement),
        coneRefinementFramePtr_
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
        coneRefinementFramePtr_,
        coneRefinementFramePtr_,
        "coneRefinementFramePtr_");
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        table1,
        "table1"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneName,
        "entry_coneName"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_cellSize,
        "entry_coneRefinement_cellSize"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p0X,
        "entry_coneRefinement_p0X"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p0Y,
        "entry_coneRefinement_p0Y"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_radius0,
        "entry_coneRefinement_radius0"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_radius1,
        "entry_coneRefinement_radius1"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p0Z,
        "entry_coneRefinement_p0Z"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p1X,
        "entry_coneRefinement_p1X"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p1Y,
        "entry_coneRefinement_p1Y"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_p1Z,
        "entry_coneRefinement_p1Z"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        comboboxentry_coneRefinement_createdCones,
        "comboboxentry_coneRefinement_createdCones"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedCellSize,
        "entry_coneRefinement_selectedCellSize"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP0X,
        "entry_coneRefinement_selectedP0X"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP0Y,
        "entry_coneRefinement_selectedP0Y"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP0Z,
        "entry_coneRefinement_selectedP0Z"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedRadius0,
        "entry_coneRefinement_selectedRadius0"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP1X,
        "entry_coneRefinement_selectedP1X"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP1Y,
        "entry_coneRefinement_selectedP1Y"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedP1Z,
        "entry_coneRefinement_selectedP1Z"
    );
    GLADE_HOOKUP_OBJECT
    (
        coneRefinementFramePtr_,
        entry_coneRefinement_selectedRadius1,
        "entry_coneRefinement_selectedRadius1"
    );
    
    updateConeRefinements(coneRefinementFramePtr_);
    
    g_signal_connect
    (
        G_OBJECT
        (
            GTK_ENTRY
            (
                GTK_COMBO(comboboxentry_coneRefinement_createdCones)->entry
            )
        ),
        "changed",
        G_CALLBACK(showDataForConeRefinement),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedCellSize),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP0X),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP0Y),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP0Z),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP1X),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP1Y),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    g_signal_connect
    (
        G_OBJECT(entry_coneRefinement_selectedP1Z),
        "changed",
        G_CALLBACK(resetValuesForSelectedCone),
        coneRefinementFramePtr_
    );
    
    gtk_widget_show(coneRefinementFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
