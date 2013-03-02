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

void meshGenGTK::createObjectRefinementMainWindowPage()
{
    //- create pages within this page
    createBoxRefinementWindowPage();
    createLineRefinementWindowPage();
    createConeRefinementWindowPage();
    createSphereRefinementWindowPage();
    
    objectRefinementMainFramePtr_ = gtk_frame_new(NULL);
    
    GtkWidget* notebook1 = gtk_notebook_new ();

    gtk_container_add(GTK_CONTAINER(objectRefinementMainFramePtr_), notebook1);
    
    //- add page for box refinement
    GtkWidget* label_boxRefinement = gtk_label_new("Box refinement");
    gtk_widget_show(label_boxRefinement);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        boxRefinementFramePtr_,
        label_boxRefinement
    );
    
    //- add page for line refinement
    GtkWidget* label_lineRefinement = gtk_label_new("Line refinement");
    gtk_widget_show(label_lineRefinement);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        lineRefinementFramePtr_,
        label_lineRefinement
    );
    
    //- add page for cone refinement
    GtkWidget* label_coneRefinement = gtk_label_new("Cone refinement");
    gtk_widget_show(label_coneRefinement);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        coneRefinementFramePtr_,
        label_coneRefinement
    );
    
    //- add page for sphere refinement
    GtkWidget* label_sphereRefinement = gtk_label_new("Sphere refinement");
    gtk_widget_show(label_sphereRefinement);
    
    gtk_notebook_append_page
    (
        GTK_NOTEBOOK(notebook1),
        sphereRefinementFramePtr_,
        label_sphereRefinement
    );
    
    //- show the notebook
    gtk_widget_show(notebook1);
    
    /* Store pointers to all widgets, for use by lookup_widget(). */
    GLADE_HOOKUP_OBJECT_NO_REF
    (
        objectRefinementMainFramePtr_,
        objectRefinementMainFramePtr_,
        "objectRefinementMainFramePtr_"
    );
    GLADE_HOOKUP_OBJECT
    (
        objectRefinementMainFramePtr_,
        notebook1,
        "notebook1"
    );

    gtk_widget_show(objectRefinementMainFramePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
