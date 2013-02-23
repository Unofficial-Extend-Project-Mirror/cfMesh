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

#include "simplexSmoother.H"
#include "tetrahedron.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simplexSmoother::simplexSmoother(partTetMeshSimplex& simplex)
:
    points_(simplex.pts()),
    tets_(simplex.tets()),
    pointI_(tets_[0][3]),
    bb_()
{
    point min(VGREAT, VGREAT, VGREAT), max(-VGREAT, -VGREAT, -VGREAT);
    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        min = Foam::min(min, tet.a());
        max = Foam::max(max, tet.a());

        min = Foam::min(min, tet.b());
        max = Foam::max(max, tet.b());

        min = Foam::min(min, tet.c());
        max = Foam::max(max, tet.c());
    }

    bb_.max() = max;
    bb_.min() = min;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simplexSmoother::~simplexSmoother()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
