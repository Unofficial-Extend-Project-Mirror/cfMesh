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

#include "triSurfFacets.H"
#include "pointIOField.H"
#include "IOobjectList.H"
#include "pointSet.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfFacets::triSurfFacets()
:
    triangles_(),
    patches_(),
    facetSubsets_()
{}

triSurfFacets::triSurfFacets(const LongList<labelledTri>& triangles)
:
    triangles_(triangles),
    patches_(1),
    facetSubsets_()
{
    forAll(triangles_, triI)
        triangles_[triI].region() = 0;

    patches_[0].name() = "patch";
}

triSurfFacets::triSurfFacets
(
    const LongList<labelledTri>& triangles,
    const geometricSurfacePatchList& patches
)
:
    triangles_(triangles),
    patches_(patches),
    facetSubsets_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
triSurfFacets::~triSurfFacets()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label triSurfFacets::addFacetSubset(const word& subsetName)
{
    label id = facetSubsetIndex(subsetName);
    if( id >= 0 )
    {
        Warning << "Point subset " << subsetName << " already exists!" << endl;
        return id;
    }

    id = 0;
    forAllConstIter(Map<meshSubset>, facetSubsets_, it)
        id = Foam::max(id, it.key()+1);

    facetSubsets_.insert
    (
        id,
        meshSubset(subsetName, meshSubset::FACESUBSET)
    );

    return id;
}

void triSurfFacets::removeFacetSubset(const label subsetID)
{
    if( facetSubsets_.find(subsetID) == facetSubsets_.end() )
        return;

    facetSubsets_.erase(subsetID);
}

word triSurfFacets::facetSubsetName(const label subsetID) const
{
    Map<meshSubset>::const_iterator it = facetSubsets_.find(subsetID);
    if( it == facetSubsets_.end() )
    {
        Warning << "Subset " << subsetID << " is not a facet subset" << endl;
        return word();
    }

    return it().name();
}

label triSurfFacets::facetSubsetIndex(const word& subsetName) const
{
    forAllConstIter(Map<meshSubset>, facetSubsets_, it)
    {
        if( it().name() == subsetName )
            return it.key();
    }

    return -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
