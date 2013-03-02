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

#include "renameBoundaryPatches.H"
#include "demandDrivenData.H"
#include "IOdictionary.H"

#include <map>

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
void renameBoundaryPatches::calculateNewBoundary()
{
    Info << "Renaming boundary patches" << endl;
    
    const dictionary& dict = meshDict_.subDict("renameBoundary");
    
    std::map<word, label> patchToLabel;
    forAll(mesh_.boundaries(), patchI)
    {
        patchToLabel.insert
        (
            std::pair<word, label>
            (
                mesh_.boundaries()[patchI].patchName(),
                patchI
            )
        );
    }
    
    labelList patchToNew(mesh_.boundaries().size(), -1);
    
    wordList newPatchNames(patchToNew.size());
    wordList newPatchTypes(patchToNew.size());
    std::map<word, label> newNameToPos;
    label newPatchI(0);
    
    //- read new patch names and types
    if( dict.found("newPatchNames") )
    {
        PtrList<entry> patchesToRename(dict.lookup("newPatchNames"));
        
        forAll(patchesToRename, patchI)
        {
            const word patchName = patchesToRename[patchI].keyword();
            
            if( patchToLabel.find(patchName) == patchToLabel.end() )
            {
                Info<< "Patch " << patchName << " does not exist!!" << endl;
                continue;
            }
            
            if( !patchesToRename[patchI].isDict() )
            {
                Warning << "Cannot rename patch " << patchName << endl;
                Warning << "This is due to incorrect settings! Exitting."
                    << endl;
                return;
            }
            
            const dictionary& pDict = patchesToRename[patchI].dict();
            
            word newName(patchName);
            if( pDict.found("newName") )
                newName = word(pDict.lookup("newName"));
            
            if( newNameToPos.find(newName) != newNameToPos.end() )
            {
                //- patch with the same name already exists
                patchToNew[patchToLabel[patchName]] = newNameToPos[newName];
                continue;
            }
            
            //- add a new patch
            newNameToPos.insert(std::pair<word, label>(newName, newPatchI));
            newPatchNames[newPatchI] = newName;
            if( pDict.found("type") )
            {
                const word newType(pDict.lookup("type"));
                newPatchTypes[newPatchI] = newType;
            }
            else
            {
                newPatchTypes[newPatchI] = "patch";
            }
            
            patchToNew[patchToLabel[patchName]] = newPatchI;
            ++newPatchI;
        }
    }
    
    word defaultName("");
    if( dict.found("defaultName") )
        defaultName = word(dict.lookup("defaultName"));
    word defaultType("patch");
    if( dict.found("defaultType") )
        defaultType = word(dict.lookup("defaultType"));
    
    if( dict.found("defaultName") )
    {
        newNameToPos.insert(std::pair<word, label>(defaultName, newPatchI));
        newPatchNames[newPatchI] = defaultName;
        newPatchTypes[newPatchI] = defaultType;
        ++newPatchI;
    }
    else
    {
        forAll(patchToNew, patchI)
        {
            if( patchToNew[patchI] != -1 )
                continue;
            
            patchToNew[patchI] = newPatchI;
            newPatchNames[newPatchI] = mesh_.boundaries()[patchI].patchName();
            newPatchTypes[newPatchI] = mesh_.boundaries()[patchI].patchType();
            ++newPatchI;
        }
    }
    
    if( newPatchI == 0 )
        return;
    
    newPatchNames.setSize(newPatchI);
    newPatchTypes.setSize(newPatchI);
    
    //- start creating new boundary
    VRWGraph newBoundaryFaces;
    labelListPMG newBoundaryOwners;
    labelListPMG newBoundaryPatches;
    
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    forAll(boundaries, patchI)
    {
        const writePatch& wp = boundaries[patchI];
        const label start = wp.patchStart();
        const label end = start + wp.patchSize();
        
        if( patchToNew[patchI] == -1 )
        {
            //- this patch is moved to the default patch
            for(label faceI=start;faceI<end;++faceI)
            {
                newBoundaryFaces.appendList(faces[faceI]);
                newBoundaryPatches.append(newPatchI-1);
                newBoundaryOwners.append(owner[faceI]);
            }
        }
        else
        {
            //- this patch is renamed
            for(label faceI=start;faceI<end;++faceI)
            {
                newBoundaryFaces.appendList(faces[faceI]);
                newBoundaryPatches.append(patchToNew[patchI]);
                newBoundaryOwners.append(owner[faceI]);
            }
        }
    }
    
    //- execute the modifier
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.replaceBoundary
    (
        newPatchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );
    forAll(meshModifier.boundariesAccess(), patchI)
        meshModifier.boundariesAccess()[patchI].patchType() =
            newPatchTypes[patchI];
    
    Info << "Finished renaming boundary patches" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

renameBoundaryPatches::renameBoundaryPatches
(
    polyMeshGen& mesh,
    const IOdictionary& meshDict
)
:
    mesh_(mesh),
    meshDict_(meshDict)
{
    if( meshDict.found("renameBoundary") )
        calculateNewBoundary();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

renameBoundaryPatches::~renameBoundaryPatches()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
