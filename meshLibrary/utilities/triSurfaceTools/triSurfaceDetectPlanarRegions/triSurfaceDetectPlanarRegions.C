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

#include "triSurfaceDetectPlanarRegions.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectPlanarRegions::checkPlanarRegions()
{
    planarRegion_.setSize(surf_.size());
    planarRegion_ = -1;

    nRegions_ = 0;

    const scalar cosTol = Foam::cos(tol_);

    const vectorField& normals = surf_.faceNormals();
    const labelListList& faceEdges = surf_.faceEdges();
    const labelListList& edgeFaces = surf_.edgeFaces();

    forAll(planarRegion_, triI)
    {
        if( planarRegion_[triI] != -1 )
            continue;

        vector n = normals[triI];
        if( magSqr(n) > VSMALL )
            n /= mag(n);

        planarRegion_[triI] = nRegions_;
        labelListPMG front;
        front.append(triI);
        label nFacetsInRegion(1);

        while( front.size() )
        {
            const label fLabel = front.removeLastElement();

            const labelList& fEdges = faceEdges[fLabel];
            forAll(fEdges, feI)
            {
                const label edgeI = fEdges[feI];

                const labelList& eFaces = edgeFaces[edgeI];

                if( eFaces.size() != 2 )
                    continue;

                forAll(eFaces, efI)
                {
                    const label fNei = eFaces[efI];

                    if( planarRegion_[fNei] != -1 )
                        continue;

                    vector neiNormal = normals[fNei];
                    if( magSqr(neiNormal) > VSMALL )
                        neiNormal /= mag(neiNormal);

                    if( mag(neiNormal & n) > cosTol )
                    {
                        planarRegion_[fNei] = nRegions_;
                        front.append(fNei);
                        ++nFacetsInRegion;
                    }
                }
            }
        }

        if( nFacetsInRegion > 1 )
        {
            ++nRegions_;
        }
        else
        {
            planarRegion_[triI] = -1;
        }
    }

    Info << "Detected " << nRegions_ << " planar regions" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectPlanarRegions::triSurfaceDetectPlanarRegions
(
    const triSurf& surface,
    const scalar angleDeviation
)
:
    surf_(surface),
    planarRegion_(),
    nRegions_(0),
    tol_(angleDeviation)
{
    if( Pstream::parRun() )
        FatalError << "Material detection does not run in parallel"
            << exit(FatalError);

    checkPlanarRegions();
}

triSurfaceDetectPlanarRegions::~triSurfaceDetectPlanarRegions()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label triSurfaceDetectPlanarRegions::numPlanarRegions() const
{
    return nRegions_;
}

void triSurfaceDetectPlanarRegions::detectedRegions(VRWGraph& graph) const
{
    labelList nFacetsInRegion(nRegions_, 0);
    forAll(planarRegion_, triI)
    {
        if( planarRegion_[triI] < 0 )
            continue;

        ++nFacetsInRegion[planarRegion_[triI]];
    }

    graph.setSizeAndRowSize(nFacetsInRegion);
    nFacetsInRegion = 0;

    forAll(planarRegion_, triI)
    {
        if( planarRegion_[triI] < 0 )
            continue;

        const label regionI = planarRegion_[triI];
        graph(regionI, nFacetsInRegion[regionI]) = triI;
    }
}

void triSurfaceDetectPlanarRegions::detectedRegions(const word prefix) const
{
    List<labelListPMG> facetsInRegion(nRegions_);
    forAll(planarRegion_, triI)
    {
        if( planarRegion_[triI] == -1 )
            continue;

        facetsInRegion[planarRegion_[triI]].append(triI);
    }

    for(label regionI=0;regionI<nRegions_;++regionI)
    {
        const word subsetName = prefix+help::scalarToText(regionI);
        const_cast<triSurf&>(surf_).addFacetsToSubset
        (
            subsetName,
            facetsInRegion[regionI]
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
