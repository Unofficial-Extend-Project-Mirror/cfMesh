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

#include "checkMeshDict.H"
#include "patchRefinementList.H"
#include "PtrList.H"
#include "objectRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkMeshDict::checkPatchCellSize()
{
    if( meshDict_.found("patchCellSize") )
    {
        patchRefinementList prl(meshDict_.lookup("patchCellSize"));
    }
}
    
void checkMeshDict::checkKeepCellsIntersectingPatches()
{
    if( meshDict_.found("keepCellsIntersectingPatches") )
    {
        wordList kcip(meshDict_.lookup("keepCellsIntersectingPatches"));
    }
}

void checkMeshDict::checkRemoveCellsIntersectingPatches()
{
    if( meshDict_.found("removeCellsIntersectingPatches") )
    {
        wordList kcip(meshDict_.lookup("removeCellsIntersectingPatches"));
    }
}
    
void checkMeshDict::checkObjectRefinements()
{
    if( meshDict_.found("objectRefinements") )
    {
        PtrList<objectRefinement> refObjects;
        Istream& is = meshDict_.lookup("objectRefinements");

        PtrList<entry> objectEntries(is);
        refObjects.setSize(objectEntries.size());
    
        forAll(refObjects, objectI)
        {
            refObjects.set
            (
                objectI,
                objectRefinement::New
                (
                    objectEntries[objectI].keyword(),
                    objectEntries[objectI].dict()
                )
            );
        }
    }
}

void checkMeshDict::checkBoundaryLayers()
{
    if( meshDict_.found("boundaryLayers") )
    {
        wordList bl(meshDict_.lookup("boundaryLayers"));
    }
}

void checkMeshDict::checkRenameBoundary()
{
    if( meshDict_.found("renameBoundary") )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        if( dict.found("newPatchNames") )
        {
            PtrList<entry> patchesToRename(dict.lookup("newPatchNames"));
        }
    }
}

void checkMeshDict::checkEntries()
{
    checkPatchCellSize();
    
    checkKeepCellsIntersectingPatches();
    
    checkRemoveCellsIntersectingPatches();
    
    checkObjectRefinements();
    
    checkBoundaryLayers();
    
    checkRenameBoundary();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

checkMeshDict::checkMeshDict
(
    const IOdictionary& meshDict
)
:
    meshDict_(meshDict)
{
    checkEntries();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

checkMeshDict::~checkMeshDict()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
