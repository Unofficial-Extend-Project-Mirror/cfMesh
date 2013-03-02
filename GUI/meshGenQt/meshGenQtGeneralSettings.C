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

#include "meshGenQtGeneralSettings.H"

#include "helperFunctions.H"

#include <qlabel.h>
#include <qlayout.h>
#include <qfiledialog.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
meshGenQtGeneralSettings::meshGenQtGeneralSettings
(
    ::QWidget* parent, meshGenGUI& meshGui
)
:
    QWidget(parent),
    meshGui_(meshGui),
    openButtonPtr_(NULL),
    surfaceFileEditPtr_(NULL),
    maxCellSizeEditPtr_(NULL),
    bndCellSizeCheckPtr_(NULL),
    bndCellSizeEditPtr_(NULL),
    autoRefineCheckPtr_(NULL),
    autoRefineCellSizeEditPtr_(NULL),
    keepBndCellsCheckPtr_(NULL),
    gluedMeshCheckPtr_(NULL)
{
    createWindow();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshGenQtGeneralSettings::~meshGenQtGeneralSettings()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenQtGeneralSettings::createWindow()
{
    QGridLayout* gridLayout = new QGridLayout(this);
    gridLayout->setSpacing(6);
    gridLayout->setMargin(11);
    
    //- set surface file name
    QLabel* surfaceLabel = new QLabel(QString("Surface file name"), this);

    gridLayout->addWidget(surfaceLabel, 0, 0);

    surfaceFileEditPtr_ = new QLineEdit(this);
    surfaceFileEditPtr_->setText(meshGui_.surfaceFileName());
    gridLayout->addWidget(surfaceFileEditPtr_, 0, 1);

    openButtonPtr_ = new QPushButton(QString("Open"), this);
    connect(openButtonPtr_, SIGNAL(clicked()), this, SLOT(openSurfaceFile()));

    gridLayout->addWidget(openButtonPtr_, 0, 3);

    //- set max cell size
    QLabel* maxCellSize = new QLabel(QString("Max cell size [m]"), this);
    gridLayout->addWidget(maxCellSize, 1, 0);

    maxCellSizeEditPtr_ = new QLineEdit(this);
    maxCellSizeEditPtr_->setText(help::scalarToText(meshGui_.maxCellSize()));
    gridLayout->addWidget(maxCellSizeEditPtr_, 1, 1);
    
    connect
    (
        maxCellSizeEditPtr_,
        SIGNAL(textChanged(const QString&)),
        this,
        SLOT(setMaxCellSize(const QString&))
    );

    QFrame* line1 = new QFrame(this);
    //line1->setObjectName(QString::fromUtf8("line1"));
    line1->setFrameShape(QFrame::HLine);
    line1->setFrameShadow(QFrame::Sunken);

    gridLayout->addMultiCellWidget(line1, 2, 2, 0, 3);

    //- bnd cell size settings
    bndCellSizeCheckPtr_ =
        new QCheckBox(QString("Use different boundary cell size"), this);

    gridLayout->addMultiCellWidget(bndCellSizeCheckPtr_, 3, 3, 0, 2);

    bndCellSizeEditPtr_ = new QLineEdit(this);
    
    if( meshGui_.boundaryCellSizeEntryExist() )
    {
        bndCellSizeCheckPtr_->setChecked(true);
        bndCellSizeEditPtr_->setText
        (
            help::scalarToText(meshGui_.boundaryCellSize())
        );
    }
    else
    {
        bndCellSizeEditPtr_->setReadOnly(true);
    }
    
    connect
    (
        bndCellSizeCheckPtr_,
        SIGNAL(clicked()),
        this,
        SLOT(activateBndCellSize())
    );
    
    connect
    (
        bndCellSizeEditPtr_,
        SIGNAL(textChanged(const QString&)),
        this,
        SLOT(setBndCellSize(const QString&))
    );
    
    gridLayout->addWidget(bndCellSizeEditPtr_, 4, 1);

    QLabel* bndCellSize = new QLabel(QString("Bnd cell size [m]"), this);

    gridLayout->addWidget(bndCellSize, 4, 0);

    QFrame* line2 = new QFrame(this);
    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);

    gridLayout->addMultiCellWidget(line2, 5, 5, 0, 3);

    //- auto refinement settings
    autoRefineCheckPtr_ = new QCheckBox("Use automatic refinement", this);
    
    connect
    (
        autoRefineCheckPtr_,
        SIGNAL(clicked()),
        this,
        SLOT(activateMinCellSize())
    );

    gridLayout->addMultiCellWidget(autoRefineCheckPtr_, 6, 6, 0, 2);

    QLabel* minCellSize = new QLabel(QString("Min cell size [m]"), this);
    gridLayout->addWidget(minCellSize, 7, 0);
    
    autoRefineCellSizeEditPtr_ = new QLineEdit(this);
    gridLayout->addWidget(autoRefineCellSizeEditPtr_, 7, 1);
    
    if( meshGui_.minCellSizeEntryExist() )
    {
        autoRefineCheckPtr_->setChecked(true);
        autoRefineCellSizeEditPtr_->setText
        (
            help::scalarToText(meshGui_.minCellSize())
        );
    }
    else
    {
        autoRefineCellSizeEditPtr_->setReadOnly(true);
    }
    
    connect
    (
        autoRefineCellSizeEditPtr_,
        SIGNAL(textChanged(const QString&)),
        this,
        SLOT(setMinCellSize(const QString&))
    );

    QFrame* line3 = new QFrame(this);
    line3->setFrameShape(QFrame::HLine);
    line3->setFrameShadow(QFrame::Sunken);

    gridLayout->addMultiCellWidget(line3, 8, 8, 0, 3);

    //- keep boundary intersected cells
    keepBndCellsCheckPtr_ =
        new QCheckBox(QString("Keep cell intersecting boundary"), this);
    
    if( meshGui_.keepCellsIntersectingBoundaryEntryExist() )
    {
        keepBndCellsCheckPtr_->setChecked(true);
    }
    else
    {
        keepBndCellsCheckPtr_->setChecked(false);
    }
    
    connect
    (
        keepBndCellsCheckPtr_,
        SIGNAL(clicked()),
        this,
        SLOT(activateKeepCellsIntersectingBoundary())
    );
    
    gridLayout->addMultiCellWidget(keepBndCellsCheckPtr_, 9, 9, 0, 2);

    //- remove glued mesh
    gluedMeshCheckPtr_ =
        new QCheckBox(QString("Remove glued mesh parts"), this);
    
    if( meshGui_.checkForGluedMeshEntryExist() )
    {
        gluedMeshCheckPtr_->setChecked(true);
    }
    else
    {
        gluedMeshCheckPtr_->setChecked(false);
    }
    
    connect
    (
        gluedMeshCheckPtr_,
        SIGNAL(clicked()),
        this,
        SLOT(activateCheckForGluedMesh())
    );
    
    gridLayout->addMultiCellWidget(gluedMeshCheckPtr_, 10, 10, 0, 2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private slots

void meshGenQtGeneralSettings::openSurfaceFile()
{
    QString fn =
        QFileDialog::getOpenFileName
        (
            QString::null,
            "Surface files (*.ftr *.stl)",
            this
        );

    if( !fn.isEmpty() )
    {
        surfaceFileEditPtr_->setText(fn);
        fileName fName(fn.ascii());
        meshGui_.setSurfaceFileName(fName);
    }
}

void meshGenQtGeneralSettings::setMaxCellSize(const QString& s)
{
    const word sn(s.ascii());
    
    const scalar size = help::textToScalar(sn);
    meshGui_.setMaxCellSize(size);
}

void meshGenQtGeneralSettings::activateBndCellSize()
{
    if( bndCellSizeCheckPtr_->isChecked() )
    {
        bndCellSizeEditPtr_->setReadOnly(false);
    }
    else
    {
        bndCellSizeEditPtr_->setText(QString(""));
        bndCellSizeEditPtr_->setReadOnly(true);
        meshGui_.removeBoundaryCellSize();
    }
}

void meshGenQtGeneralSettings::setBndCellSize(const QString& s)
{
    if( bndCellSizeCheckPtr_->isChecked() )
    {
        const word sn(s.ascii());
    
        const scalar size = help::textToScalar(sn);
        meshGui_.setBoundaryCellSize(size);
    }
}

void meshGenQtGeneralSettings::activateMinCellSize()
{
    if( autoRefineCheckPtr_->isChecked() )
    {
        autoRefineCellSizeEditPtr_->setReadOnly(false);
    }
    else
    {
        autoRefineCellSizeEditPtr_->setText(QString(""));
        autoRefineCellSizeEditPtr_->setReadOnly(true);
        meshGui_.removeMinCellSize();
    }
}

void meshGenQtGeneralSettings::setMinCellSize(const QString& s)
{
    if( bndCellSizeCheckPtr_->isChecked() )
    {
        const word sn(s.ascii());
    
        const scalar size = help::textToScalar(sn);
        meshGui_.setMinCellSize(size);
    }
}

void meshGenQtGeneralSettings::activateKeepCellsIntersectingBoundary()
{
    if( keepBndCellsCheckPtr_->isChecked() )
    {
        meshGui_.setKeepCellsIntersectingBoundary();
    }
    else
    {
        meshGui_.removeKeepCellsIntersectingBoundary();
        
        if( gluedMeshCheckPtr_->isChecked() )
        {
            gluedMeshCheckPtr_->setChecked(false);
            meshGui_.removeCheckForGluedMesh();
        }
    }
}
        
void meshGenQtGeneralSettings::activateCheckForGluedMesh()
{
    if( keepBndCellsCheckPtr_->isChecked() )
    {
        meshGui_.setCheckForGluedMesh();
        
        if( gluedMeshCheckPtr_->isChecked() )
        {
            meshGui_.setCheckForGluedMesh();
        }
        else
        {
            meshGui_.removeCheckForGluedMesh();
        }
    }
    else
    {
        meshGui_.removeCheckForGluedMesh();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
