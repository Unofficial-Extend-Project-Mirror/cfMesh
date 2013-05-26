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
        if( meshDict_.isDict("patchCellSize") )
        {
            const dictionary& dict = meshDict_.subDict("patchCellSize");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            patchRefinementList prl(meshDict_.lookup("patchCellSize"));
            prl.size();
        }
    }
}

void checkMeshDict::checkSubsetCellSize()
{
    if( meshDict_.found("subsetCellSize") )
    {
        if( meshDict_.isDict("subsetCellSize") )
        {
            const dictionary& dict = meshDict_.subDict("subsetCellSize");

            const wordList subsetNames = dict.toc();
            subsetNames.size();
        }
        else
        {
            patchRefinementList prl(meshDict_.lookup("patchCellSize"));
        }
    }
}

void checkMeshDict::checkKeepCellsIntersectingPatches()
{
    if( meshDict_.found("keepCellsIntersectingPatches") )
    {
        if( meshDict_.isDict("keepCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("keepCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            wordList kcip(meshDict_.lookup("keepCellsIntersectingPatches"));
        }
    }
}

void checkMeshDict::checkRemoveCellsIntersectingPatches()
{
    if( meshDict_.found("removeCellsIntersectingPatches") )
    {
        if( meshDict_.isDict("removeCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("removeCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            wordList kcip(meshDict_.lookup("removeCellsIntersectingPatches"));
        }
    }
}

void checkMeshDict::checkObjectRefinements()
{
    if( meshDict_.found("objectRefinements") )
    {
        PtrList<objectRefinement> refObjects;

        if( meshDict_.isDict("objectRefinements") )
        {
            const dictionary& dict = meshDict_.subDict("objectRefinements");
            const wordList objectNames = dict.toc();

            refObjects.setSize(objectNames.size());

            forAll(refObjects, objectI)
            {
                const entry& objectEntry =
                    dict.lookupEntry(objectNames[objectI], false, false);

                refObjects.set
                (
                    objectI,
                    objectRefinement::New
                    (
                        objectEntry.keyword(),
                        objectEntry.dict()
                    )
                );
            }
        }
        else
        {
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
}

void checkMeshDict::checkBoundaryLayers()
{
    if( meshDict_.found("boundaryLayers") )
    {
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& dict = meshDict_.subDict("boundaryLayers");

            const wordList layerNames = dict.toc();
            layerNames.size();
        }
        else
        {
            wordList bl(meshDict_.lookup("boundaryLayers"));
        }
    }
}

void checkMeshDict::checkRenameBoundary()
{
    if( meshDict_.found("renameBoundary") )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        if( dict.found("newPatchNames") )
        {
            if( dict.isDict("newPatchNames") )
            {
                const dictionary& patchDicts = dict.subDict("newPatchNames");

                patchDicts.toc();
            }
            else
            {
                const PtrList<entry> patchesToRename
                (
                    dict.lookup("newPatchNames")
                );

                patchesToRename.size();
            }
        }
    }
}

void checkMeshDict::checkEntries()
{
    checkPatchCellSize();

    checkSubsetCellSize();

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
