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
    meshSurfaceDetectPlanarRegions

Description
    Detects planar regions in the surface of the volume mesh

SourceFiles
    meshSurfaceDetectPlanarRegions.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceDetectPlanarRegions_H
#define meshSurfaceDetectPlanarRegions_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceDetectPlanarRegions Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceDetectPlanarRegions
{
    // Private data
        //- const reference to mesh surface
        const meshSurfaceEngine& meshSurface_;

        //- angle deviation tolerance
        scalar angleTol_;

        //- planar region index fo each boundary face
        labelList faceInPlanarRegion_;

        //- number of detected planar regions
        label nPlanarRegions_;

    // Private member functions
        //- calculate new boundary and replace the existing one
        void findPlanarRegions();

        //- Disallow default bitwise copy construct
        meshSurfaceDetectPlanarRegions(const meshSurfaceDetectPlanarRegions&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceDetectPlanarRegions&);

public:

    // Constructors

        //- Construct from mesh and IOdictionary
        meshSurfaceDetectPlanarRegions
        (
            const meshSurfaceEngine& meshSurface,
            const scalar angleTol = (M_PI / 180.0)
        );

    // Destructor

        ~meshSurfaceDetectPlanarRegions();

    // Member Functions
        //- return the number of detected planar regions
        inline label nPlanarRegions() const
        {
            return nPlanarRegions_;
        }

        //- returns a list containing indices of planar regions
        inline const labelList& planarRegions() const
        {
            return faceInPlanarRegion_;
        }

        //- return a region index for a given face
        inline label planarRegion(const label bfI) const
        {
            return faceInPlanarRegion_[bfI];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
