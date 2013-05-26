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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

#include <omp.h>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::addProcessorFaces
(
    const VRWGraph& procFaces,
    const labelListPMG& facePatches
)
{
    Info << "Adding processor faces" << endl;

    PtrList<writeProcessorPatch>& procBoundaries = mesh_.procBoundaries_;

    labelList nAddedFaces(procBoundaries.size(), 0);
    forAll(facePatches, fI)
        ++nAddedFaces[facePatches[fI]];

    labelList newPatchStart(procBoundaries.size());
    newPatchStart[0] = procBoundaries[0].patchStart();
    for(label i=1;i<procBoundaries.size();++i)
        newPatchStart[i] =
            newPatchStart[i-1] +
            procBoundaries[i-1].patchSize() + nAddedFaces[i-1];

    //- set new size of the faceListPMG
    faceListPMG& faces = mesh_.faces_;
    const label nFaces = faces.size();
    faces.setSize(nFaces+procFaces.size());

    label endProcFaces(0);
    forAllReverse(procBoundaries, patchI)
    {
        const writeProcessorPatch& wp = procBoundaries[patchI];
        endProcFaces = Foam::max(endProcFaces, wp.patchStart()+wp.patchSize());
    }

    //- move faces to their new positions
    labelListPMG newFaceLabel(nFaces, -1);

    if( endProcFaces != nFaces )
    {
        for(label faceI=nFaces-1;faceI>=endProcFaces;--faceI)
        {
            newFaceLabel[faceI] = faceI+facePatches.size();
            faces[faceI+facePatches.size()].transfer(faces[faceI]);
        }
    }

    labelList faceIndex(procBoundaries.size());
    for(label patchI=procBoundaries.size()-1;patchI>=0;--patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        const label shift = newPatchStart[patchI] - start;

        if( shift != 0 )
        {
            for(label faceI=end-1;faceI>=start;--faceI)
            {
                faces[faceI+shift].transfer(faces[faceI]);
                newFaceLabel[faceI] = faceI+shift;
            }
        }

        //- set new start for the given patch
        procBoundaries[patchI].patchStart() = newPatchStart[patchI];
        faceIndex[patchI] =
            newPatchStart[patchI] + procBoundaries[patchI].patchSize();
        procBoundaries[patchI].patchSize() += nAddedFaces[patchI];
    }

    //- add new faces into patches
    forAll(procFaces, fI)
    {
        face f(procFaces.sizeOfRow(fI));
        forAll(f, pI)
            f[pI] = procFaces(fI, pI);

        faces[faceIndex[facePatches[fI]]++].transfer(f);
    }

    //- renumber cells
    cellListPMG& cells = mesh_.cells_;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            if( newFaceLabel[c[fI]] != -1 )
                c[fI] = newFaceLabel[c[fI]];
    }

    this->clearOut();
    mesh_.clearOut();
    mesh_.updateFaceSubsets(newFaceLabel);

    Info << "Finished adding processor faces" << endl;
}

label polyMeshGenModifier::addProcessorPatch(const label otherProcLabel)
{
    const label nProcPatches = mesh_.procBoundaries().size();

    PtrList<writeProcessorPatch>& procBoundaries =
        this->procBoundariesAccess();

    procBoundaries.setSize(nProcPatches + 1);

    std::ostringstream ss;
    ss << Pstream::myProcNo();
    std::ostringstream ssNei;
    ssNei << otherProcLabel;
    const word name("processor"+ss.str()+"to"+ssNei.str());

    procBoundaries.set
    (
        nProcPatches,
        new writeProcessorPatch
        (
            name,
            "processor",
            0,
            0,
            Pstream::myProcNo(),
            otherProcLabel
        )
    );

    return nProcPatches;
}

bool polyMeshGenModifier::removeEmptyProcessorPatches()
{
    PtrList<writeProcessorPatch>& procBoundaries =
        this->procBoundariesAccess();

    label nValidPatches(0);
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].patchSize() != 0 )
            ++nValidPatches;
    }

    if( nValidPatches == procBoundaries.size() )
        return false;

    PtrList<writeProcessorPatch> newProcBoundaries(nValidPatches);

    nValidPatches = 0;
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].patchSize() != 0 )
        {
            newProcBoundaries.set
            (
                nValidPatches++,
                new writeProcessorPatch(procBoundaries[patchI])
            );
        }
    }

    procBoundaries.transfer(newProcBoundaries);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
