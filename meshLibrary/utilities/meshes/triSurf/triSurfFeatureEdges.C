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

#include "triSurfFeatureEdges.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfFeatureEdges::triSurfFeatureEdges()
:
    featureEdges_(),
    featureEdgeSubsets_()
{}

triSurfFeatureEdges::triSurfFeatureEdges(const edgeLongList& featureEdges)
:
    featureEdges_(featureEdges),
    featureEdgeSubsets_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
triSurfFeatureEdges::~triSurfFeatureEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label triSurfFeatureEdges::addEdgeSubset(const word& subsetName)
{
    label id = edgeSubsetIndex(subsetName);
    if( id >= 0 )
    {
        Warning << "Edge subset " << subsetName << " already exists!" << endl;
        return id;
    }

    id = 0;
    forAllConstIter(Map<meshSubset>, featureEdgeSubsets_, it)
        id = Foam::max(id, it.key()+1);

    featureEdgeSubsets_.insert
    (
        id,
        meshSubset(subsetName, meshSubset::FEATUREEDGESUBSET)
    );

    return id;
}

void triSurfFeatureEdges::removeEdgeSubset(const label subsetID)
{
    if( featureEdgeSubsets_.find(subsetID) == featureEdgeSubsets_.end() )
        return;

    featureEdgeSubsets_.erase(subsetID);
}

word triSurfFeatureEdges::edgeSubsetName(const label subsetID) const
{
    Map<meshSubset>::const_iterator it = featureEdgeSubsets_.find(subsetID);
    if( it == featureEdgeSubsets_.end() )
    {
        Warning << "Subset " << subsetID << " is not an edge subset" << endl;
        return word();
    }

    return it().name();
}

label triSurfFeatureEdges::edgeSubsetIndex(const word& subsetName) const
{
    forAllConstIter(Map<meshSubset>, featureEdgeSubsets_, it)
    {
        if( it().name() == subsetName )
            return it.key();
    }

    return -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
