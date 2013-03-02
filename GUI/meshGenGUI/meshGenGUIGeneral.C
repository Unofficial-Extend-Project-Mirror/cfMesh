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

#include "meshGenGUI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshGenGUI::setSurfaceFileName(const fileName& fName)
{
    if( meshDict_.found("surfaceFile") )
        meshDict_.remove("surfaceFile");
    
    meshDict_.add("surfaceFile", fName);
}

fileName meshGenGUI::surfaceFileName() const
{
    if( !meshDict_.found("surfaceFile") )
    {
        fileName f("");
        const_cast<IOdictionary&>(meshDict_).add("surfaceFile", f);
    }
    
    return meshDict_.lookup("surfaceFile");
}

void meshGenGUI::setMaxCellSize(const scalar s)
{
    if( meshDict_.found("maxCellSize") )
        meshDict_.remove("maxCellSize");
    
    meshDict_.add("maxCellSize", s);
}

scalar meshGenGUI::maxCellSize() const
{
    if( !meshDict_.found("maxCellSize") )
        const_cast<IOdictionary&>(meshDict_).add("maxCellSize", 1.0);
    
    return readScalar(meshDict_.lookup("maxCellSize"));
}
        
bool meshGenGUI::boundaryCellSizeEntryExist() const
{
    word w("boundaryCellSize");
    return meshDict_.found(w);
}

void meshGenGUI::removeBoundaryCellSize()
{
    meshDict_.remove("boundaryCellSize");
}

void meshGenGUI::setBoundaryCellSize(const scalar s)
{
    if( boundaryCellSizeEntryExist() )
        meshDict_.remove("boundaryCellSize");
    
    meshDict_.add("boundaryCellSize", s);
}

scalar meshGenGUI::boundaryCellSize() const
{
    return readScalar(meshDict_.lookup("boundaryCellSize"));
}
        
bool meshGenGUI::minCellSizeEntryExist() const
{
    return meshDict_.found("minCellSize");
}

void meshGenGUI::removeMinCellSize()
{
    if( minCellSizeEntryExist() )
        meshDict_.remove("minCellSize");
}

void meshGenGUI::setMinCellSize(const scalar s)
{
    removeMinCellSize();
    
    meshDict_.add("minCellSize", s);
}

scalar meshGenGUI::minCellSize() const
{
    return readScalar(meshDict_.lookup("minCellSize"));
}
        
bool meshGenGUI::keepCellsIntersectingBoundaryEntryExist() const
{
    return meshDict_.found("keepCellsIntersectingBoundary");
}

void meshGenGUI::setKeepCellsIntersectingBoundary()
{
    removeKeepCellsIntersectingBoundary();
    
    meshDict_.add("keepCellsIntersectingBoundary", 1);
}

void meshGenGUI::removeKeepCellsIntersectingBoundary()
{
    if( keepCellsIntersectingBoundaryEntryExist() )
        meshDict_.remove("keepCellsIntersectingBoundary");
}

bool meshGenGUI::keepCellsIntersectingBoundary() const
{
    if( !keepCellsIntersectingBoundaryEntryExist() )
        return false;
    
    return meshDict_.lookup("keepCellsIntersectingBoundary");
}
        
bool meshGenGUI::checkForGluedMeshEntryExist() const
{
    return meshDict_.found("checkForGluedMesh");
}

void meshGenGUI::setCheckForGluedMesh()
{
    removeCheckForGluedMesh();
    
    meshDict_.add("checkForGluedMesh", 1);
}

void meshGenGUI::removeCheckForGluedMesh()
{
    if( checkForGluedMeshEntryExist() )
        meshDict_.remove("checkForGluedMesh");
}

bool meshGenGUI::checkForGluedMesh() const
{
    if( !checkForGluedMeshEntryExist() )
        return false;
    
    return readBool(meshDict_.lookup("checkForGluedMesh"));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
