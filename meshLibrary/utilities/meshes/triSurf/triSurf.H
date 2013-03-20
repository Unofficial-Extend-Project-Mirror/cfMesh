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

Class
    triSurf

Description
    A class for triangulated surface used in the meshing process.
    It is derived from points and facets with some additional subsets.
    Subsets are vert useful for local mesh refinement.

SourceFiles
    triSurf.C

\*---------------------------------------------------------------------------*/

#ifndef triSurf_H
#define triSurf_H

#include "triSurfPoints.H"
#include "triSurfFacets.H"
#include "triSurfAddressing.H"
#include <map>
#include "DynList.H"
#include "labelListPMG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Froward declarations
class triSurfModifier;

/*---------------------------------------------------------------------------*\
                        Class triSurf Declaration
\*---------------------------------------------------------------------------*/

class triSurf
:
    public triSurfPoints,
    public triSurfFacets,
    public triSurfAddressing
{
    // Private data

        //- groups of facets
        //- which can be visited over manifold edges
        //mutable label nFaceGroups_;
        //mutable labelListPMG* faceGroupPtr_;

    // Private member functions
        //- calculate groups of facets
        //void calculateFaceGroups() const;

        void readFromSTL(const fileName&);
        void writeToSTL(const fileName&) const;

        void readFromSTLB(const fileName&);
        void writeToSTLB(const fileName&) const;

        void readFromFTR(const fileName&);
        void writeToFTR(const fileName&) const;

        inline LongList<labelledTri>& accessToFacets();
        inline geometricSurfacePatchList& accessToPatches();

        //- Disallow default bitwise assignment
        void operator=(const triSurf&);

public:

    // Friend classes
        friend class triSurfModifer;

    // Constructors

        //- Default construct
        triSurf();

        //- Construct from parts
        triSurf
        (
            const LongList<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            const pointField& points
        );

        //- Read from file
        triSurf(const fileName& fName);

    // Destructor

        ~triSurf();


    // Member Functions
        //- read and write the surface
        void readSurface(const fileName&);
        void writeSurface(const fileName&) const;

/*        //- returns the names of the existing subsets
        inline void existingFaceSubsets(DynList<word>&) const;

        //- return facets in subset
        inline bool doesFaceSubsetExist(const word&) const;
        inline const labelListPMG& facesInSubset(const word&) const;

        //- creates a new subset if it does not yet exist
        //- adds elements which are not yet in the subset
        inline void addFacetsToSubset(const word&, const labelListPMG&);

        //- remove subset of facets
        inline void removeFaceSubset(const word&);

        //- number of groups of facets that are connected over manifold edges
        inline label numberOfFaceGroups() const;

        //- groups assigned to the surface facets
        inline const labelListPMG& faceGroups() const;
        */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
