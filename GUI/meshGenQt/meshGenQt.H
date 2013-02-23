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
    meshGenQt

Description
    Graphical user inteface using Qt

SourceFiles
    meshGenQt.C

\*---------------------------------------------------------------------------*/

#ifndef meshGenQt_H
#define meshGenQt_H

#include "meshGenGUI.H"

#include <qvariant.h>
#include "qwidget.h"
#include <qapplication.h>
#include <qmainwindow.h>
#include <qpushbutton.h>
#include <qtabwidget.h>
#include <qaction.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class triSurf;

class meshGenQtGeneralSettings;

/*---------------------------------------------------------------------------*\
                           Class meshGenQt Declaration
\*---------------------------------------------------------------------------*/

class meshGenQt
:
	public QMainWindow
{
		Q_OBJECT
	
	// Private data
	
		//- mesh generation dictionary
		meshGenGUI meshGui_;
	
		//- main window data
		QPushButton* saveButtonPtr_;
		QPushButton* quitButtonPtr_;
		QTabWidget* tabWidgetPtr_;
		meshGenQtGeneralSettings* tabGeneralSettingsPtr_;
		QWidget* tabRefinementPtr_;
		QWidget* tabPatchRenamePtr_;
	
private slots:
	
	// Private callback functions

		//- save meshDict file
		void saveMeshDictCallback()
		{
			meshGui_.writeDict();
		}

private:
	
	// Private member functions
	
		//- create main window
		void createMainWindow();

        //- Disallow default bitwise copy construct
        meshGenQt(const meshGenQt&);

        //- Disallow default bitwise assignment
        void operator=(const meshGenQt&);

public:

    // Constructors

        //- Construct from registry
        meshGenQt(label argc, char* argv[], const objectRegistry&);

    // Destructor

        ~meshGenQt();

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
