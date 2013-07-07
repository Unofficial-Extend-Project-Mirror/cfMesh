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

#include "error.H"
#include "meshSurfaceEngine.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctionsPar.H"
#include "DynList.H"
#include "labelPair.H"
#include "HashSet.H"

#include <omp.h>

#define DEBUGExchangeMap

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void distributeBoundaryFaces
(
    const meshSurfaceEngine& surfaceEngine,
    const meshOctree& octree,
    labelListPMG& facePatch
)
{
    const faceList::subList& bFaces = surfaceEngine.boundaryFaces();
    const pointFieldPMG& points = surfaceEngine.points();

    //- set the size of the facePatch list
    facePatch.setSize(bFaces.size());

    //- set size of patchNames, newBoundaryFaces_ and newBoundaryOwners_
    const triSurf& surface = octree.surface();
    const label nPatches = surface.patches().size();

    //- find the region for face by finding the patch nearest
    //- to the face centre
    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 100 ) schedule(dynamic, 40)
    # endif
    forAll(bFaces, bfI)
    {
        const point c = bFaces[bfI].centre(points);

        label fPatch;
        point p;
        scalar distSq;

        octree.findNearestSurfacePoint(p, distSq, fPatch, c);

        if( (fPatch > -1) && (fPatch < nPatches) )
        {
            facePatch[bfI] = fPatch;
        }
        else
        {
            FatalErrorIn
            (
                "void meshSurfaceEdgeExtractorNonTopo::"
                "distributeBoundaryFaces()"
            ) << "Cannot distribute a face " << bFaces[bfI] << " into any "
                << "surface patch!. Exiting.." << exit(FatalError);
        }
    }
}

void findCornerCandidates
(
    const meshSurfaceEngine& surfaceEngine,
    const meshOctree& octree,
    LongList<labelPair>& pointMap
)
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
