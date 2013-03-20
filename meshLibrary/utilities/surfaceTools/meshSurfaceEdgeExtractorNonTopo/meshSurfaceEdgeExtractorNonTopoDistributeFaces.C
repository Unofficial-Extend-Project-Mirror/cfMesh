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

#include "demandDrivenData.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "helperFunctions.H"

#include <omp.h>

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEdgeExtractorNonTopo::distributeBoundaryFaces()
{
    meshSurfaceEngine mse(mesh_);

    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwner = mse.faceOwners();
    const pointFieldPMG& points = mse.points();

    //- set size of patchNames, newBoundaryFaces_ and newBoundaryOwners_
    const triSurf& surface = meshOctree_.surface();
    const label nPatches = surface.patches().size();

    wordList patchNames(nPatches);
    VRWGraph newBoundaryFaces;
    labelListPMG newBoundaryOwners(bFaces.size());
    labelListPMG newBoundaryPatches(bFaces.size());

    //- set patchNames
    forAll(surface.patches(), patchI)
            patchNames[patchI] = surface.patches()[patchI].name();

    //- append boundary faces
    forAll(bFaces, bfI)
    {
        newBoundaryFaces.appendList(bFaces[bfI]);
        newBoundaryOwners[bfI] = faceOwner[bfI];
    }

    //- find the region for face by finding the patch nearest
    //- to the face centre
    # pragma omp parallel for if( bFaces.size() > 100 ) schedule(guided)
    forAll(bFaces, bfI)
    {
        const point c = bFaces[bfI].centre(points);

        label facePatch;
        point p;
        scalar distSq;

        meshOctree_.findNearestSurfacePoint(p, distSq, facePatch, c);

        if( (facePatch > -1) && (facePatch < nPatches) )
        {
            newBoundaryPatches[bfI] = facePatch;
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

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );
}

void meshSurfaceEdgeExtractorNonTopo::remapBoundaryPoints()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceMapper mapper(mse, meshOctree_);

    mapper.mapVerticesOntoSurfacePatches();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
