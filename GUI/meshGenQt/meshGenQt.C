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

#include "meshGenQt.H"
#include "objectRegistry.H"
#include "Time.H"

#include "meshGenQtGeneralSettings.H"

#include <qlayout.h>
#include <qvbox.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
meshGenQt::meshGenQt
(
    label argc, char* argv[], const objectRegistry& reg
)
:
	QMainWindow(),
	meshGui_(reg),
	saveButtonPtr_(NULL),
	quitButtonPtr_(NULL),
	tabWidgetPtr_(NULL),
	tabGeneralSettingsPtr_(NULL),
	tabRefinementPtr_(NULL),
	tabPatchRenamePtr_(NULL)
{
	createMainWindow();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshGenQt::~meshGenQt()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenQt::createMainWindow()
{
	QWidget* widget = new QWidget(this);
	widget->setGeometry(QRect(0, 0, 600, 480));
	
	//- create and add save button
	saveButtonPtr_ = new QPushButton(QString("Save meshDict"), widget);
	saveButtonPtr_->setGeometry(QRect(11, 11, 150, 24));
	connect
	(
		saveButtonPtr_,
		SIGNAL(clicked()),
		this,
		SLOT(saveMeshDictCallback())
	);
	
	//- create and add quit button
	quitButtonPtr_ = new QPushButton(QString("Quit"), widget);
	quitButtonPtr_->setGeometry(QRect(160, 11, 180, 24));
	connect(quitButtonPtr_, SIGNAL(clicked()), qApp, SLOT(quit()));

	//- create the table
    tabWidgetPtr_ = new QTabWidget(widget);
	tabWidgetPtr_->setGeometry(QRect(9, 41, 580, 406));

	tabGeneralSettingsPtr_ = new meshGenQtGeneralSettings(widget, meshGui_);
	tabWidgetPtr_->addTab(tabGeneralSettingsPtr_, "General settings");
	
	QWidget* tab1 = new QWidget();
	tabWidgetPtr_->addTab(tab1, "Refinement");
	
	//- set the front page 
	this->setCentralWidget(widget);
	
	//- resize to the desired size
	QSize size(600, 480);
    size = size.expandedTo(this->minimumSizeHint());
    this->resize(size);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
