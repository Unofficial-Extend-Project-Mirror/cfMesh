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

Class
    meshGenGTK

Description
    Graphical user inteface using GTK

SourceFiles
    meshGenGTK.C

\*---------------------------------------------------------------------------*/

#ifndef meshGenGTK_H
#define meshGenGTK_H

#include "meshGenGUI.H"

#include <gtk/gtk.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class triSurface;
    
static meshGenGUI* guiPtr(NULL);

/*---------------------------------------------------------------------------*\
                           Class meshGenGTK Declaration
\*---------------------------------------------------------------------------*/

class meshGenGTK
{
    // Private data
    
        //- mesh generation dictionary
        meshGenGUI meshGui_;
    
        //- pages
        GtkWidget* generalPageFramePtr_;
        GtkWidget* mainWindowPtr_;
        GtkWidget* localRefinementFramePtr_;
        GtkWidget* keepCellsIntersectingPatchesFramePtr_;
        GtkWidget* bndLayersFramePtr_;
        GtkWidget* localSettingsMainFramePtr_;

        GtkWidget* boxRefinementFramePtr_;
        GtkWidget* lineRefinementFramePtr_;
        GtkWidget* coneRefinementFramePtr_;
        GtkWidget* sphereRefinementFramePtr_;
        GtkWidget* objectRefinementMainFramePtr_;
    
        GtkWidget* renameBoundaryFramePtr_;

private:
    
    // Private member functions
        //- create pages
        void createGeneralPage();
        void createMainWindowPage();

        //- pages for local features
        void createLocalRefinementWindowPage();
        void createKeepCellsIntersectingPatchesWindowPage();
        void createBndLayersWindowPage();
        void createLocalSettingsMainWindowPage();

        //- pages for object refinement
        void createBoxRefinementWindowPage();
        void createLineRefinementWindowPage();
        void createConeRefinementWindowPage();
        void createSphereRefinementWindowPage();
        void createObjectRefinementMainWindowPage();

        //- pages for boundary renaming
        void createRenameBoundaryMainWindowPage();

        //- Disallow default bitwise copy construct
        meshGenGTK(const meshGenGTK&);

        //- Disallow default bitwise assignment
        void operator=(const meshGenGTK&);

public:

    // Constructors

        //- Construct from registry
        meshGenGTK(label argc, char* argv[], const objectRegistry&);

    // Destructor

        ~meshGenGTK();

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define GLADE_HOOKUP_OBJECT(component,widget,name) \
  g_object_set_data_full(G_OBJECT(component), name, \
    gtk_widget_ref(widget), (GDestroyNotify)gtk_widget_unref)

#define GLADE_HOOKUP_OBJECT_NO_REF(component,widget,name) \
  g_object_set_data(G_OBJECT(component), name, widget)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
