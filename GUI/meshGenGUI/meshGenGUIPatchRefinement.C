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

bool meshGenGUI::patchCellSizeEntryExist() const
{
    return meshDict_.found("patchCellSize");
}

void meshGenGUI::addPatchCellSize(const patchRefinement& pr)
{
    if( patchCellSizeEntryExist() )
    {
        List<patchRefinement> prList(meshDict_.lookup("patchCellSize"));
        
        const label size = prList.size();
        prList.setSize(size + 1);
        prList[size] = pr;
        
        meshDict_.remove("patchCellSize");
        meshDict_.add("patchCellSize", prList);
    }
    else
    {
        List<patchRefinement> prList(1);
        prList[0] = pr;
        
        meshDict_.add("patchCellSize", prList);
    }
}

void meshGenGUI::removePatchCellSize(const word& patchName)
{
    if( patchCellSizeEntryExist() )
    {
        List<patchRefinement> prList(meshDict_.lookup("patchCellSize"));
        
        label pos(-1);
        forAll(prList, elI)
            if( prList[elI].patchName() == patchName )
            {
                pos = elI;
                break;
            }
            
        prList[pos] = prList[prList.size() - 1];
        prList.setSize(prList.size() - 1);
        
        meshDict_.remove("patchCellSize");
        meshDict_.add("patchCellSize", prList);
    }
}

List<patchRefinement> meshGenGUI::patchCellSize() const
{
    if( !patchCellSizeEntryExist() )
        return List<patchRefinement>();
    
    List<patchRefinement> pr(meshDict_.lookup("patchCellSize"));
    return pr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
