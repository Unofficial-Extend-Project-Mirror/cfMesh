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

#include "boundaryFacesGenerator.H"
#include "createFacesFromChain.H"

//#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private member functions * * * * * * * * * * * //

void boundaryFacesGenerator::createFacesForChain
(
    const labelList& chainVertices,
    const DynList<label>& patches,
    List< DynList<face> >& facesForChain
)
{
    createFacesFromChain cffc(chainVertices, pointRegions_);

    //- create faces which do not need an additional corner first
    cffc.createFacesWithoutACorner();

    if( cffc.unresolvedPoints().size() == 0 )
    {
        //- store created faces and return
        const DynList<face>& createdFaces = cffc.createdFaces();
        const DynList<label>& faceRegion = cffc.faceRegion();

        forAll(createdFaces, fI)
        {
            const label pos = patches.containsAtPosition(faceRegion[fI]);
            const face& f = createdFaces[fI];
            # ifdef DEBUGCutter
            Info << "Storing face " << f << " at position " << pos << endl;
            # endif
            facesForChain[pos].append(f);
        }

        return;
    }
    else if( cffc.unresolvedPoints().size() < 3 )
    {
        FatalErrorIn
        (
            "void boundaryFacesGenerator::createFacesForChain"
            "("
            "const labelList& chainVertices,"
            "const DynList<label>& patches,"
            "List< DynList<face> >& facesForChain"
            ")"
        ) << "I am not sure if this should ever happen!" << exit(FatalError);
    }

    const labelList& unresolvedPoints = cffc.unresolvedPoints();
    DynList<label> chainPatches;
    forAll(unresolvedPoints, upI)
    {
        const labelList pp = patchesForPoint(unresolvedPoints[upI]);

        forAll(pp, ppI)
            chainPatches.appendIfNotIn(pp[ppI]);
    }

    if( chainPatches.size() < 3 )
    {
        //- expected three or more patches
        //- this will be resolved by storing the face into a default patch
        facesForChain[patches.size()].append(face(chainVertices));
    }
    else
    {
        //- find the cornerLabel
        DynList<label> cornersCandidates;

        forAll(cornersPatches_, cornerI)
        {
            const DynList<label>& cPatches = cornersPatches_[cornerI];

            bool allFound(true);
            forAll(cPatches, cpI)
            {
                bool found = chainPatches.contains(cPatches[cpI]);

                if( !found )
                {
                    allFound = false;
                    break;
                }
            }

            if( allFound )
                cornersCandidates.append(surfaceCorners_[cornerI]);
        }

        if( cornersCandidates.size() != 0 )
        {
            label cornerLabel(-1);

            if( cornersCandidates.size() == 1 )
            {
                cornerLabel = cornersCandidates[0];
            }
            else
            {
                //- find nearest corner
                point c(vector::zero);
                forAll(chainVertices, cvI)
                    c += mesh_.points()[chainVertices[cvI]];
                c /= chainVertices.size();

                scalar dist(VGREAT);
                forAll(cornersCandidates, cornerI)
                    if(
                        mag
                        (
                            surface_.points()[cornersCandidates[cornerI]] -
                            c
                        ) < dist
                    )
                    {
                        dist =
                        mag
                        (
                            surface_.points()[cornersCandidates[cornerI]] -
                            c
                        );
                        cornerLabel = cornersCandidates[cornerI];
                    }
            }

            //- create faces including the corner
            //- create faces with a corner vertex
            cffc.createFacesWithACorner(nPoints_);

            //- add new mesh vertex
            polyMeshGenModifier modifier(mesh_);
            modifier.pointsAccess().newElmt(nPoints_++) =
                surface_.points()[cornerLabel];

            const DynList<face>& createdFaces = cffc.createdFaces();
            const DynList<label>& faceRegion = cffc.faceRegion();
            forAll(createdFaces, fI)
            {
                const label pos = patches.containsAtPosition(faceRegion[fI]);
                const face& f = createdFaces[fI];
                # ifdef DEBUGCutter
                Info << "Storing face " << f << " at position " << pos << endl;
                # endif
                facesForChain[pos].append(f);
            }
        }
        else
        {
            //- this is a tricky combination
            //- store the face into the default patch
            facesForChain[patches.size()].append(face(chainVertices));
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
