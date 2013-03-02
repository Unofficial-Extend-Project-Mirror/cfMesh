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
    Creates surface patches from surface subsets

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurf.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurf* makePatchFromSubset
(
    triSurf& origSurf,
    const DynList<word>& subsetNames
)
{
    //- create new list of patches
    geometricSurfacePatchList newPatches
    (
        origSurf.patches().size() + subsetNames.size()
    );
    
    //- set names of the new patches
    forAll(origSurf.patches(), patchI)
        newPatches[patchI].name() = origSurf.patches()[patchI].name();
    
    forAll(subsetNames, subsetI)
        newPatches[origSurf.patches().size()+subsetI].name() =
            subsetNames[subsetI];
    
    //- create new triangles
    List<labelledTri> newTriangles(origSurf.localFaces());
    
    //- set patches for all triangles
    forAll(subsetNames, subsetI)
    {
        const labelListPMG& subsetFaces =
            origSurf.facesInSubset(subsetNames[subsetI]);
        
        const label regionI = origSurf.patches().size() + subsetI;
    
        forAll(subsetFaces, fI)
        {
            newTriangles[subsetFaces[fI]].region() = regionI;
        }
    }
    
    //- remove patches with no elements
    labelList nTrianglesInPatch(newPatches.size(), 0);
    forAll(newTriangles, triI)
        ++nTrianglesInPatch[newTriangles[triI].region()];
    
    Map<label> newPatchLabel;
    label counter(0);
    forAll(nTrianglesInPatch, patchI)
    {
        if( nTrianglesInPatch[patchI] )
            newPatchLabel.insert(patchI, counter++);
    }
    
    geometricSurfacePatchList copyPatches(counter);
    counter = 0;
    forAll(newPatches, patchI)
    {
        if( newPatchLabel.found(patchI) )
        {
            copyPatches[newPatchLabel[patchI]].name() =
                newPatches[patchI].name();
        }
    }
    
    newPatches = copyPatches;
    
    //- renumber the patches in the list of triangles
    forAll(newTriangles, triI)
        newTriangles[triI].region() =
            newPatchLabel[newTriangles[triI].region()];
    
    //- copy subsets for the new surface
    DynList<word> existingSubsets;
    origSurf.existingFaceSubsets(existingSubsets);
    
    std::map<word, labelListPMG> newSubsets;
    forAll(existingSubsets, sI)
        newSubsets.insert
        (
            std::pair<word, labelListPMG>
            (
                existingSubsets[sI],
                origSurf.facesInSubset(existingSubsets[sI])
            )
        );
    
    //- delete subsets converted to patches
    forAll(subsetNames, subsetI)
        newSubsets.erase(subsetNames[subsetI]);
    
    triSurf* newSurfPtr =
        new triSurf(newTriangles, newPatches, origSurf.points(), newSubsets);
    
    return newSurfPtr;
}
    

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("subset file name");
    argList::validArgs.append("subset");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName subsetFileName(args.args()[2]);
    word subsetName(args.args()[3]);
    
    triSurf* origSurfPtr = new triSurf(inFileName);
    origSurfPtr->readFaceSubsets(subsetFileName);
    
    DynList<word> subsetNames;
    if( !origSurfPtr->doesFaceSubsetExist(subsetName) )
    {
        Warning << "Subset " << subsetName
        << " checking subsets containing this string!" << endl;
        DynList<word> existingSubsets;
        origSurfPtr->existingFaceSubsets(existingSubsets);
        
        forAll(existingSubsets, subsetI)
        {
            if(
                existingSubsets[subsetI].substr(0, subsetName.size()) ==
                subsetName
            )
            {
                subsetNames.append(existingSubsets[subsetI]);
            }
        }
        
        Info << "Converting " << subsetNames.size() << " subsets" << endl;
    }
    else
    {
        subsetNames.append(subsetName);
    }
    
    triSurf* newSurfPtr = makePatchFromSubset(*origSurfPtr, subsetNames);
    deleteDemandDrivenData(origSurfPtr);
        
    newSurfPtr->write(inFileName, false);
    newSurfPtr->writeFaceSubsets(subsetFileName);
    deleteDemandDrivenData(newSurfPtr);
    
    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
