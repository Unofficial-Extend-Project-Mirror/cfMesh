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
    triSurfaceDetectFeatureEdges

Description
    Divides the surface mesh into regions bounded by feature edges

SourceFiles
    triSurfaceDetectFeatureEdges.C
    triSurfaceDetectFeatureEdgesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDetectFeatureEdges_H
#define triSurfaceDetectFeatureEdges_H

#include "triSurf.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceDetectFeatureEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDetectFeatureEdges
{
    // Private data
        //- reference to triSurf
        const triSurf& surf_;

        //- detected feature edges
        List<direction> featureEdges_;

        //- angle tolerance
        scalar angleTolerance_;

        //- surface patches
        labelList facetInPatch_;

        //- number of patches
        label nPatches_;

        //- patch names
        wordList newPatchNames_;

        //- patch types
        wordList newPatchTypes_;

    // Private member functions
        //- detect feature edges by check angle between the normals
        //- of the two neighbouring facets
        void detectFeatureEdgesAngleCriterion();

        //- detect feature edges by checking the angle between
        //- the point normals in the vicinity of an edge
        void detectFeatureEdgesPointAngleCriterion();

        //- detect feature edges at the boundaries of planar regions
        void detectOuterBoundariesOfPlanarRegions();

        //- create patches bounded by a set of feature edges
        void createPatches();

        //- Disallow default bitwise copy construct
        triSurfaceDetectFeatureEdges(const triSurfaceDetectFeatureEdges&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceDetectFeatureEdges&);

public:

    // Constructors

        //- Construct from triSurface and tolerance
        triSurfaceDetectFeatureEdges
        (
            const triSurf& surface,
            const scalar angleDeviation = 45.0
        );

    // Destructor

        ~triSurfaceDetectFeatureEdges();


    // Member Functions

        //- store regions in a graph
        void detectedSurfaceRegions(VRWGraph&) const;

        //- store regions into subsets with a given prefix
        const triSurf* surfaceWithPatches
        (
            const word prefix = "patch_",
            const bool forceOverwrite = false
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
