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

#include "triSurfaceDetectFeatureEdges.H"
#include "helperFunctions.H"
#include "triSurfaceDetectPlanarRegions.H"
#include "demandDrivenData.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::detectFeatureEdgesAngleCriterion()
{
    const scalar tol = Foam::cos(angleTolerance_*M_PI/180.0);

    const vectorField& normals = surf_.facetNormals();

    const VRWGraph& edgeFaces = surf_.edgeFacets();

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(edgeFaces, edgeI)
    {
        const constRow eFaces = edgeFaces[edgeI];

        if( edgeFaces.sizeOfRow(edgeI) != 2 )
        {
            featureEdges_[edgeI] |= 8;
            continue;
        }

        const scalar cosAngle =
            (normals[eFaces[0]] & normals[eFaces[1]]) /
            (mag(normals[eFaces[0]]) * mag(normals[eFaces[1]]) + VSMALL);

        if( cosAngle < tol )
            featureEdges_[edgeI] |= 1;
    }
}

void triSurfaceDetectFeatureEdges::detectFeatureEdgesPointAngleCriterion()
{

}

void triSurfaceDetectFeatureEdges::detectOuterBoundariesOfPlanarRegions()
{
    triSurfaceDetectPlanarRegions dpr(surf_, 0.1);

    VRWGraph planarRegions;
    dpr.detectedRegions(planarRegions);

    labelListPMG facetInPlanarRegion(surf_.size(), -1);

    forAll(planarRegions, regionI)
        forAllRow(planarRegions, regionI, rfI)
            facetInPlanarRegion[planarRegions(regionI, rfI)] = regionI;

    const VRWGraph& edgeFaces = surf_.edgeFacets();

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(edgeFaces, edgeI)
    {
        const constRow eFaces = edgeFaces[edgeI];

        if( eFaces.size() != 2 )
            continue;

        if( facetInPlanarRegion[eFaces[0]] != facetInPlanarRegion[eFaces[1]] )
            featureEdges_[edgeI] |= 4;
    }
}

void triSurfaceDetectFeatureEdges::createPatches()
{
    nPatches_ = 0;
    facetInPatch_.setSize(surf_.size());
    facetInPatch_ = -1;

    const VRWGraph& faceEdges = surf_.facetEdges();
    const VRWGraph& edgeFaces = surf_.edgeFacets();

    forAll(facetInPatch_, triI)
    {
        if( facetInPatch_[triI] != -1 )
            continue;

        labelListPMG front;
        front.append(triI);
        facetInPatch_[triI] = nPatches_;

        while( front.size() )
        {
            const label fLabel = front.removeLastElement();

            const constRow fEdges = faceEdges[fLabel];

            forAll(fEdges, feI)
            {
                const label edgeI = fEdges[feI];

                //- check if th edges is marked as a feature edge
                if( featureEdges_[edgeI] )
                    continue;

                const constRow eFaces = edgeFaces[edgeI];

                //- stop at non-manifold edges
                if( eFaces.size() != 2 )
                    continue;

                label neiTri = eFaces[0];
                if( neiTri == fLabel )
                    neiTri = eFaces[1];

                //- do not overwrite existing patch information
                if( surf_[fLabel].region() != surf_[neiTri].region() )
                    continue;
                if( facetInPatch_[neiTri] != -1 )
                    continue;

                facetInPatch_[neiTri] = nPatches_;
                front.append(neiTri);
            }
        }

        ++nPatches_;
    }

    Info << "Created " << nPatches_ << " surface patches" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
