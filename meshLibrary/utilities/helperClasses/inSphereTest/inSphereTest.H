/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

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
