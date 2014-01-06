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

#include "demandDrivenData.H"
#include "volumeOptimizer.H"
#include "tetrahedron.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volumeOptimizer::volumeOptimizer(partTetMeshSimplex& simplex)
:
    simplexSmoother(simplex)
{
}

volumeOptimizer::~volumeOptimizer()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

void volumeOptimizer::optimizeNodePosition(const scalar tol)
{
    label iter(0);

    point& p = points_[pointI_];

    if( !bb_.contains(p) )
        p = 0.5 * (bb_.max() + bb_.min());

    const scalar scale = 1.0 / bb_.mag();
    forAll(points_, pI)
        points_[pI] *= scale;

    # ifdef DEBUGSmooth
    Info << nl << "Smoothing point " << pointI_
         << " with coordinates " << p << endl;
    scalar Vmina(VGREAT);
    forAll(tets_, tetI)
    Vmina = Foam::min(Vmina, tets_[tetI].mag(points_));
    Info << "Vmin before " << Vmina << endl;
    # endif

    vector gradF;
    vector disp(vector::zero);
    tensor gradGradF;
    point pOrig;

    scalar funcBefore, funcAfter(evaluateFunc());

    bool finished;
    do
    {
        finished = false;
        pOrig = p;
        funcBefore = funcAfter;

        evaluateGradientsExact(gradF, gradGradF);

        const scalar determinant = Foam::det(gradGradF);
        if( determinant > SMALL )
        {
            disp = (inv(gradGradF, determinant) & gradF);

            p -= disp;

            funcAfter = evaluateFunc();

            # ifdef DEBUGSmooth
            Info << nl << "gradF " << gradF << endl;
            Info << "gradGradF " << gradGradF << endl;
            Info << "det(gradGradF) " << determinant << endl;
            Info << "disp " << disp << endl;
            Info << "Func before " << funcBefore << endl;
            Info << "Func after " << funcAfter << endl;
            # endif

            scalar relax(0.8);
            label nLoops(0);

            while( funcAfter > funcBefore )
            {
                p = pOrig - relax * disp;
                relax *= 0.5;
                funcAfter = evaluateFunc();

                if( funcAfter < funcBefore )
                    continue;

                if( ++nLoops == 5 )
                {
                    //- it seems that this direction is wrong, stop the loop
                    p = pOrig;
                    disp = vector::zero;
                    finished = true;
                    funcAfter = funcBefore;
                }
            }

            if( mag(funcBefore - funcAfter) / funcBefore < tol )
                finished = true;
        }
        else
        {
            //- move in random direction
            //- this is usually needed to move the point off the zero volume
            disp = vector::zero;
            forAll(tets_, tetI)
            {
                const partTet& tet = tets_[tetI];
                const scalar Vtri = tet.mag(points_);

                if( Vtri < SMALL )
                {
                    triangle<point, point> tri
                    (
                        points_[tet.a()],
                        points_[tet.b()],
                        points_[tet.c()]
                    );

                    vector n = tri.normal();
                    const scalar d = mag(n);

                    if( d > VSMALL )
                        disp += 0.01 * (n / d);
                }
            }

            p += disp;
            funcAfter = evaluateFunc();
        }
    } while( (++iter < 10) && !finished );

    # ifdef DEBUGSmooth
    scalar Vmin(VGREAT);
    forAll(tets_, tetI)
        Vmin = Foam::min(Vmin, tets_[tetI].mag(points_));

    Info << nl << "New coordinates for point "
        << pointI_ << " are " << p << endl;
    Info << "Num iterations " << iter << " gradient " << gradF << endl;
    Info << "Vmin " << Vmin << endl;
    # endif

    //- scale back to the original size
    forAll(points_, pI)
        points_[pI] /= scale;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
