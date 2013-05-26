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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfAddressing.H"

#include <omp.h>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const VRWGraph& triSurfAddressing::pointFacets() const
{
    if( !pointFacetsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& triSurfAddressing::pointFacets() const"
            ) << "Cannot calculate pointFacets" << abort(FatalError);
        # endif

        calculatePointFacets();
    }

    return *pointFacetsPtr_;
}

inline const edgeListPMG& triSurfAddressing::edges() const
{
    if( !edgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const edgeListPMG& triSurfAddressing::edges() const"
            ) << "Cannot calculate edges" << abort(FatalError);
        # endif

        calculateEdges();
    }

    return *edgesPtr_;
}

inline const VRWGraph& triSurfAddressing::facetEdges() const
{
    if( !facetEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& triSurfAddressing::facetEdges() const"
            ) << "Cannot calculate facetEdges" << abort(FatalError);
        # endif

        calculateFacetEdges();
    }

    return *facetEdgesPtr_;
}

inline const VRWGraph& triSurfAddressing::edgeFacets() const
{
    if( !edgeFacetsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& triSurfAddressing::edgeFacets() const"
            ) << "Cannot calculate edgeFacets" << abort(FatalError);
        # endif

        calculateEdgeFacets();
    }

    return *edgeFacetsPtr_;
}

inline const VRWGraph& triSurfAddressing::pointEdges() const
{
    if( !pointEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& triSurfAddressing::pointEdges() const"
            ) << "Cannot calculate pointEdges" << abort(FatalError);
        # endif

        calculatePointEdges();
    }

    return *pointEdgesPtr_;
}

inline const VRWGraph& triSurfAddressing::facetFacets() const
{
    if( !facetFacetsEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& triSurfAddressing::facetFacets() const"
            ) << "Cannot calculate facetFacets" << abort(FatalError);
        # endif

        calculateFacetFacetsEdges();
    }

    return *facetFacetsEdgesPtr_;
}

inline const vectorField& triSurfAddressing::pointNormals() const
{
    if( !pointNormalsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const vectorField&"
                " triSurfAddressing::pointNormals() const"
            ) << "Cannot calculate pointNormals" << abort(FatalError);
        # endif

        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}

inline const vectorField& triSurfAddressing::facetNormals() const
{
    if( !facetNormalsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const vectorField&"
                " triSurfAddressing::facetNormals() const"
            ) << "Cannot calculate facetNormals" << abort(FatalError);
        # endif

        calculateFacetNormals();
    }

    return *facetNormalsPtr_;
}

inline const vectorField& triSurfAddressing::facetCentres() const
{
    if( !facetCentresPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const vectorField&"
                " triSurfAddressing::facetCentres() const"
            ) << "Cannot calculate facetCentres" << abort(FatalError);
        # endif

        calculateFacetCentres();
    }

    return *facetCentresPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
