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
    triSurfaceDetectPlanarRegions

Description
    Provides information which parts of the surface mesh are planar

SourceFiles
    triSurfaceDetectPlanarRegions.C
    triSurfaceDetectPlanarRegionsFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDetectPlanarRegions_H
#define triSurfaceDetectPlanarRegions_H

#include "triSurf.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceDetectPlanarRegions Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDetectPlanarRegions
{
    // Private data
        //- reference to triSurf
        const triSurf& surf_;

        //- planar region to which a surface facet belongs to
        labelList planarRegion_;

        //- number of planar regions
        label nRegions_;

        //- tolerance
        scalar tol_;

    // Private member functions
        //- check if the partitions are well behaved
        void checkPlanarRegions();

        //- Disallow default bitwise copy construct
        triSurfaceDetectPlanarRegions(const triSurfaceDetectPlanarRegions&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceDetectPlanarRegions&);

public:

    // Constructors

        //- Construct from triSurface and tolerance
        triSurfaceDetectPlanarRegions
        (
            const triSurf& surface,
            const scalar angleDeviation = 0.5 * M_PI / 180.0
        );

    // Destructor

        ~triSurfaceDetectPlanarRegions();


    // Member Functions

        //- return the number of planar regions
        label numPlanarRegions() const;

        //- store regions in a graph
        void detectedRegions(VRWGraph&) const;

        //- store regions into subsets with a given prefix
        void detectedRegions(const word prefix = "planarRegion_") const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
