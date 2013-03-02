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
#include "surfaceOptimizer.H"
#include "matrix2D.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
surfaceOptimizer::surfaceOptimizer
(
    DynList<point>& pts,
    const DynList<triFace>& trias
)
:
    pts_(pts),
    trias_(trias)
{
}

surfaceOptimizer::~surfaceOptimizer()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point surfaceOptimizer::optimizePoint(const scalar tol)
{
    scalar avgEdge(0.0);
    forAll(trias_, triI)
        avgEdge +=
            mag(pts_[trias_[triI][1]] - pts_[trias_[triI][0]]);
    avgEdge /= trias_.size();

    avgEdge *= tol;
    
    tensor gradGradLt(tensor::zero);
    gradGradLt.xx() = 4.0;
    gradGradLt.yy() = 4.0;
    
    point newPoint(vector::zero);
    
    vector disp;
    disp.z() = 0.0;
    direction nIterations(0);
    do
    {
        //- find the minimum area
        scalar Amin(SMALL);
        forAll(trias_, triI)
        {
            const point& p0 = pts_[trias_[triI][0]];
            const point& p1 = pts_[trias_[triI][1]];
            const point& p2 = pts_[trias_[triI][2]];
            
            const scalar Atri =
                0.5 *
                (
                    (p1.y() - p2.y()) * p0.x() +
                    (p2.x() - p1.x()) * p0.y() +
                    (p1.x() * p2.y() - p2.x() * p1.y())
                );
            
            if( Atri < Amin )
                Amin = Atri;
        }
        
        //- K is greater than zero in case the stabilisation is needed
        scalar K = 0.0;
        if( Amin < SMALL )
        {
            K = max(10.0 * sqrt(SMALL * (SMALL - Amin)), K);
            if( sqrt(sqr(Amin) + 4.0 * sqr(K)) < (1.0 + SMALL) * mag(Amin) )
                K = max(sqrt(SMALL) * mag(Amin), K);
        }
        //const scalar K = sqrt(SMALL * (SMALL - Amin));
        
        # ifdef DEBUGSmooth
        Info << "Point " << bpI << " K = " << K << endl;
        # endif
        
        //- start assembling the system
        vector gradF(vector::zero);
        tensor gradGradF(tensor::zero);
        forAll(trias_, triI)
        {
            const point& p0 = pts_[trias_[triI][0]];
            const point& p1 = pts_[trias_[triI][1]];
            const point& p2 = pts_[trias_[triI][2]];
            
            if( magSqr(p1 - p2) < VSMALL ) continue;
            
            const scalar LSqrTri
            (
                magSqr(p0 - p1) +
                magSqr(p1 - p2) +
                magSqr(p2 - p0)
            );
            
            const scalar Atri =
                0.5 *
                (
                    (p1.y() - p2.y()) * p0.x() +
                    (p2.x() - p1.x()) * p0.y() +
                    (p1.x() * p2.y() - p2.x() * p1.y())
                );
            
            const scalar stab = sqrt(sqr(Atri) + K);
            const scalar Astab = 0.5 * (Atri + stab);
            const vector B((p1.y() - p2.y()), (p2.x() - p1.x()), 0.0);
            const vector gradAstab = 0.5 * (B + (Atri * B) / stab);
            const tensor gradGradAstab =
                0.5 *
                (
                    (B * B) / stab -
                    (B * B) * sqr(Atri) / pow(stab, 3)
                );
                
            const vector gradLt(4 * p0 - 2.0 * p1 - 2.0 * p2);
            
            //- calculate the gradient
            gradF += (gradLt * Astab - LSqrTri * gradAstab) / sqr(Astab);
        
            //- calculate the second gradient
            gradGradF +=
                gradGradLt / Astab -
                2.0 * symm(gradLt * gradAstab) / sqr(Astab) -
                LSqrTri * gradGradAstab / sqr(Astab) +
                2.0 * LSqrTri * (gradAstab * gradAstab) / pow(Astab, 3);
        }
        
        if( mag(gradGradF.xx()) < VSMALL ) gradGradF.xx() = VSMALL;
        if( mag(gradGradF.yy()) < VSMALL ) gradGradF.yy() = VSMALL;
        
        matrix2D mat;
        mat[0][0] = gradGradF.xx();
        mat[0][1] = gradGradF.xy();
        mat[1][0] = gradGradF.yx();
        mat[1][1] = gradGradF.yy();
        FixedList<scalar, 2> source;
        source[0] = gradF.x();
        source[1] = gradF.y();
        
        const scalar det = mat.determinant();
        
        if( mag(det) < VSMALL )
        {
            disp = vector::zero;
        }
        else
        {
            disp.x() = mat.solveFirst(source);
            disp.y() = mat.solveSecond(source);
            
            if( mag(disp) > 0.7 * avgEdge )
            {
                vector dir = disp / mag(disp);
            
                disp = dir * 0.7 * avgEdge;
            }
        }
        
        # ifdef DEBUGSmooth
        Info << "Second gradient " << gradGradF << endl;
        Info << "Gradient " << gradF << endl;
        Info << "Displacement " << disp << endl;
        # endif
        
        newPoint -= disp;
        pts_[trias_[0][0]] = newPoint;

        #ifdef DEBUGSmooth
        Info << "New coordinates " << newPoint << endl;
        # endif

    } while( (++nIterations < 100) && (mag(disp) > avgEdge) );
    
    //- do not move the vertex if the process has not converged
    if( nIterations >= 100 )
        newPoint = vector::zero;
    
    return newPoint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
