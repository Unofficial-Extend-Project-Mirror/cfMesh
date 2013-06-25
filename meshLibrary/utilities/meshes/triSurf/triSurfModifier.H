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
    triSurfModifier

Description
    A class allowing non-const access to the member of the surface mesh

SourceFiles
    triSurfModifier.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfModifier_H
#define triSurfModifier_H

#include "triSurf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class triSurfModifier Declaration
\*---------------------------------------------------------------------------*/

class triSurfModifier
{
    // Private data
        //- reference to the surface mesh
        triSurf& surface_;

    // Disallow bitwise assignment
        //- Default construct
        triSurfModifier();

        //- Disallow default bitwise assignment
        void operator=(const triSurfModifier&);

public:

    // Constructors

        //- Construct from surface mesh
        triSurfModifier(triSurf& surface);

    // Destructor

        ~triSurfModifier();


    // Member Functions
        //- non-const access to points
        inline pointField& pointsAccess();

        //- access to facets
        inline LongList<labelledTri>& facetsAccess();

        //- non-const access to feature edges
        inline edgeListPMG& featureEdgesAccess();

        //- access to patches
        inline geometricSurfacePatchList& patchesAccess();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
