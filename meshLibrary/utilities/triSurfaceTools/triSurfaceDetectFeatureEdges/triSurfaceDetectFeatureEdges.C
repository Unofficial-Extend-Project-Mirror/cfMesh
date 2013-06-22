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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectFeatureEdges::triSurfaceDetectFeatureEdges
(
    const triSurf& surface,
    const scalar angleDeviation
)
:
    surf_(surface),
    featureEdges_(surf_.edges().size(), direction(0)),
    angleTolerance_(angleDeviation),
    facetInPatch_(),
    nPatches_(),
    newPatchNames_(),
    newPatchTypes_()
{
    if( Pstream::parRun() )
        FatalError << "Material detection does not run in parallel"
            << exit(FatalError);

    detectFeatureEdgesAngleCriterion();

    detectFeatureEdgesPointAngleCriterion();

    detectOuterBoundariesOfPlanarRegions();

    createPatches();
}

triSurfaceDetectFeatureEdges::~triSurfaceDetectFeatureEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::detectedSurfaceRegions
(
    VRWGraph& graph
) const
{
    graph.setSize(nPatches_);

    labelListPMG nFacetsInPatch(nPatches_, 0);

    forAll(facetInPatch_, triI)
        ++nFacetsInPatch[facetInPatch_[triI]];

    graph.setSizeAndRowSize(nFacetsInPatch);

    nFacetsInPatch = 0;
    forAll(facetInPatch_, triI)
    {
        const label patchI = facetInPatch_[triI];

        graph(patchI, nFacetsInPatch[patchI]) = triI;
        ++nFacetsInPatch[patchI];
    }
}

const triSurf* triSurfaceDetectFeatureEdges::surfaceWithPatches
(
    const word prefix,
    const bool forceOverwrite
) const
{
    //- collect patch information
    VRWGraph facetsInPatch;
    detectedSurfaceRegions(facetsInPatch);

    //- create new list of boundary patches
    LongList<labelledTri> newTriangles(facetInPatch_.size());
    label counter(0);
    geometricSurfacePatchList newPatches(nPatches_);

    if( forceOverwrite )
    {
        forAll(newPatches, patchI)
        {
            newPatches[patchI].name() = prefix+help::scalarToText(patchI);
            newPatches[patchI].geometricType() = "patch";
            newPatches[patchI].index() = patchI;
        }
    }
    else
    {
        forAll(facetsInPatch, patchI)
        {
            forAllRow(facetsInPatch, patchI, fpI)
            {
                const label origPatchI =
                    surf_[facetsInPatch(patchI, fpI)].region();
                newPatches[patchI].name() =
                    surf_.patches()[origPatchI].name() +
                    help::scalarToText(patchI);
                newPatches[patchI].geometricType() =
                    surf_.patches()[origPatchI].geometricType();
                newPatches[patchI].index() = patchI;
            }
        }
    }

    //- create triangles for the new surface
    labelListPMG newFacetLabel(newTriangles.size(), -1);

    forAll(facetsInPatch, patchI)
        forAllRow(facetsInPatch, patchI, tI)
        {
            newFacetLabel[facetsInPatch(patchI, tI)] = counter;
            labelledTri tria = surf_[facetsInPatch(patchI, tI)];
            tria.region() = patchI;
            newTriangles[counter++] = tria;
        }

    //- create and return a new surface mesh
    triSurf* newSurfPtr =
        new triSurf
        (
            newTriangles,
            newPatches,
            edgeListPMG(),
            surf_.points()
        );

    //- transfer facet subsets
    DynList<label> subsetIDs;
    surf_.facetSubsetIndices(subsetIDs);
    forAll(subsetIDs, subsetI)
    {
        const word sName = surf_.facetSubsetName(subsetIDs[subsetI]);

        const label newID = newSurfPtr->addFacetSubset(sName);

        labelListPMG facetsInSubset;
        surf_.facetsInSubset(subsetIDs[subsetI], facetsInSubset);

        forAll(facetsInSubset, i)
        {
            const label fI = newFacetLabel[facetsInSubset[i]];

            newSurfPtr->addFacetToSubset(newID, fI);
        }
    }

    //- transfer point subsets
    surf_.pointSubsetIndices(subsetIDs);
    forAll(subsetIDs, subsetI)
    {
        const word sName = surf_.pointSubsetName(subsetIDs[subsetI]);

        const label newID = newSurfPtr->addPointSubset(sName);

        labelListPMG pointsInSubset;
        surf_.pointsInSubset(subsetIDs[subsetI], pointsInSubset);

        forAll(pointsInSubset, i)
            newSurfPtr->addPointToSubset(newID, pointsInSubset[i]);
    }

    return newSurfPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
