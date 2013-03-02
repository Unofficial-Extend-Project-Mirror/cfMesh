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
    meshGenQtGeneralSettings

Description
    General setting page

SourceFiles
    meshGenQtGeneralSettings.C

\*---------------------------------------------------------------------------*/

#ifndef meshGenQtGeneralSettings_H
#define meshGenQtGeneralSettings_H

#include "meshGenGUI.H"

#include <qvariant.h>
#include <qwidget.h>
#include <qapplication.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qcheckbox.h>
#include <qaction.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations

/*---------------------------------------------------------------------------*\
                           Class meshGenQtGeneralSettings Declaration
\*---------------------------------------------------------------------------*/

class meshGenQtGeneralSettings
:
    public QWidget
{
        Q_OBJECT
    
    // Private data
    
        //- mesh generation dictionary
        meshGenGUI& meshGui_;
    
        //- main window data
        QPushButton* openButtonPtr_;
        QLineEdit* surfaceFileEditPtr_;
        QLineEdit* maxCellSizeEditPtr_;
    
        QCheckBox* bndCellSizeCheckPtr_;
        QLineEdit* bndCellSizeEditPtr_;
    
        QCheckBox* autoRefineCheckPtr_;
        QLineEdit* autoRefineCellSizeEditPtr_;
        
        QCheckBox* keepBndCellsCheckPtr_;
        QCheckBox* gluedMeshCheckPtr_;
    
private slots:
    
    // Private callback functions

        //- open surface file
        void openSurfaceFile();
        
        //- set max cell size
        void setMaxCellSize(const QString&);

        //- set activity of bndCellSize
        void activateBndCellSize();

        //- set bndCellSize
        void setBndCellSize(const QString&);

        //- set activity of auto refinement
        void activateMinCellSize();

        //- set min size in auto refinement
        void setMinCellSize(const QString&);

        //- set activity of keepCellsIntersectingBoundary
        void activateKeepCellsIntersectingBoundary();
        
        //- set activity of glued mesh check
        void activateCheckForGluedMesh();
        
private:
    
    // Private member functions
    
        //- create window
        void createWindow();

        //- Disallow default bitwise copy construct
        meshGenQtGeneralSettings(const meshGenQtGeneralSettings&);

        //- Disallow default bitwise assignment
        void operator=(const meshGenQtGeneralSettings&);

public:

    // Constructors

        //- Construct parent and meshGenGUI
        meshGenQtGeneralSettings(::QWidget* parent, meshGenGUI& meshGui);

    // Destructor

        virtual ~meshGenQtGeneralSettings();

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
