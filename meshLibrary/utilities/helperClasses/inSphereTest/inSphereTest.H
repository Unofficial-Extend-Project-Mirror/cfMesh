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
    inSphereTest

Description
    Calculates if a point is inside a tetrahedron sphere or not. Uses adaptive-
	precision arithmetic by J. R. Schewchuk

SourceFiles
    inSphereTest.C

\*---------------------------------------------------------------------------*/

#ifndef inSphereTest_H
#define inSphereTest_H

#include "triSurface.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inSphereTest Declaration
\*---------------------------------------------------------------------------*/

class inSphereTest
{
    // Private members
        const point& pa_;
        const point& pb_;
        const point& pc_;
        const point& pd_;

    // Private member functions

public:

    // Constructors
        inline inSphereTest
        (
            const point&,
            const point&,
            const point&,
            const point&
        );

    // Destructor
        ~inSphereTest();


    // Member functions

    scalar inSphereFast(const point& pe) const;
    scalar inSphere(const point& pe) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
