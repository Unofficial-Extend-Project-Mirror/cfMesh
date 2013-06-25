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
#include "demandDrivenData.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectFeatureEdges::triSurfaceDetectFeatureEdges
(
    triSurf& surface,
    const scalar angleDeviation
)
:
    surf_(surface),
    featureEdges_(surf_.edges().size(), direction(0)),
    angleTolerance_(angleDeviation)
{
    if( Pstream::parRun() )
        FatalError << "Feature edges detection does not run in parallel"
            << exit(FatalError);

    detectFeatureEdgesAngleCriterion();

    detectFeatureEdgesPointAngleCriterion();

    //detectOuterBoundariesOfPlanarRegions();
}

triSurfaceDetectFeatureEdges::~triSurfaceDetectFeatureEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::detectFeatureEdges()
{
    const edgeListPMG& edges = surf_.edges();
    triSurfModifier surfMod(surf_);
    edgeListPMG& featureEdges = surfMod.featureEdgesAccess();
    featureEdges.clear();

    forAll(featureEdges_, eI)
    {
        if( featureEdges_[eI] )
            featureEdges.append(edges[eI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
