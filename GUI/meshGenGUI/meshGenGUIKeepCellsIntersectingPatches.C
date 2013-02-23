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

bool meshGenGUI::keepCellsIntersectingPatchesEntryExist() const
{
	return meshDict_.found("keepCellsIntersectingPatches");
}

void meshGenGUI::addKeepCellsIntersectingPatches(const word& pName)
{
	if( keepCellsIntersectingPatchesEntryExist() )
	{
		wordList pList(meshDict_.lookup("keepCellsIntersectingPatches"));
		
		const label size = pList.size();
		pList.setSize(size + 1);
		pList[size] = pName;
		
		meshDict_.remove("keepCellsIntersectingPatches");
		meshDict_.add("keepCellsIntersectingPatches", pList);
	}
	else
	{
		wordList pList(1);
		pList[0] = pName;
		
		meshDict_.add("keepCellsIntersectingPatches", pList);
	}
}

void meshGenGUI::removeKeepCellsIntersectingPatches(const word& pName)
{
	if( keepCellsIntersectingPatchesEntryExist() )
	{
		wordList pList(meshDict_.lookup("keepCellsIntersectingPatches"));
		
		label pos(-1);
		forAll(pList, elI)
			if( pList[elI] == pName )
			{
				pos = elI;
				break;
			}
			
		pList[pos] = pList[pList.size() - 1];
		pList.setSize(pList.size() - 1);
		
		meshDict_.remove("keepCellsIntersectingPatches");
		
		if( pList.size() != 0 )
			meshDict_.add("keepCellsIntersectingPatches", pList);
	}
}

wordList meshGenGUI::keepCellsIntersectingPatches() const
{
	if( !keepCellsIntersectingPatchesEntryExist() )
		return wordList();
	
	wordList pNames(meshDict_.lookup("keepCellsIntersectingPatches"));
	return pNames;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
