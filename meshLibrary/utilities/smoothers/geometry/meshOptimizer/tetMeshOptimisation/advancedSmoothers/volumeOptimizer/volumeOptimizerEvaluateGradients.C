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

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar volumeOptimizer::evaluateFunc() const
{
    const scalar K = evaluateStabilisationFactor();
    
	scalar func(0.0);
	
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
		
		const scalar LSqrTri
		(
			magSqr(tet.d() - tet.a()) +
			magSqr(tet.d() - tet.b()) +
			magSqr(tet.d() - tet.c())
		);
		
		const scalar Vtri = tet.mag();
        const scalar Vstab = 0.5 * (Vtri + mag(Vtri)) + K;
		
		func += LSqrTri / pow(Vstab, 2./3.);
	}
	
	return func;
}

scalar volumeOptimizer::evaluateStabilisationFactor() const
{
	scalar K = 0.0;
	
	scalar Vmin(VGREAT), LSqMax(0.0);
	
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
		
		const scalar Vtri = tet.mag();
        
        Vmin = Foam::min(Vmin, Vtri);
        
        const scalar LSqrTri
		(
			magSqr(tet.d() - tet.a()) +
			magSqr(tet.d() - tet.b()) +
			magSqr(tet.d() - tet.c())
		);
        
        LSqMax = Foam::max(LSqMax, LSqrTri);
	}
	
	if( Vmin < SMALL * LSqMax )
        K = SMALL * LSqMax;
	
	return K;
}

void volumeOptimizer::evaluateGradientsExact
(
	vector& gradF,
	tensor& gradGradF
) const
{
	gradF = vector::zero;
	gradGradF = tensor::zero;
    
    const scalar K = evaluateStabilisationFactor();
	
	tensor gradGradLsq(tensor::zero);
	gradGradLsq.xx() = 6.0;
	gradGradLsq.yy() = 6.0;
	gradGradLsq.zz() = 6.0;
	
	const point& p = points_[pointI_];
	
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
		
		const vector B
		(
			(1.0/6.0) *
			(
				(tet.b() - tet.a()) ^
				(tet.c() - tet.a())
			)
		);
		
		//- add this trihedron to the metric
		const scalar LSqrTri
		(
			magSqr(tet.d() - tet.a()) +
			magSqr(tet.d() - tet.b()) +
			magSqr(tet.d() - tet.c())
		);
	
		const scalar Vtri = tet.mag();
	
		//- evaluate gradients
		const scalar Vstab = 0.5 * (Vtri + mag(Vtri)) + K;
		
		if( Vstab < VSMALL )
		{
			Info << "Tet " << tet << endl;
			Info << "B " << B << endl;
			Info << "Vtri " << Vtri << endl;
			IOstream::defaultPrecision(20);
			Info << "Vstab " << Vstab << endl;

			FatalErrorIn
			(
				"void nodeDisplacementVolumeOptimizer()"
			) << "I cannot continue " << exit(FatalError);
		}
		
		const vector gradLsq = 2. * (3. * p - tet.a() - tet.b() - tet.c());
		const vector gradVstab = 0.5 * (B + Foam::sign(Vtri) * B);

		//- calculate the gradients
		const scalar Vs = pow(Vstab, 2./3.) + VSMALL;
        const scalar Vs53 = Vs * Vstab + VSMALL;

		gradF += (gradLsq / Vs) - (2./3. * LSqrTri * gradVstab / Vs53);
	
		//- calculate the second gradient
		gradGradF +=
			10./9. * LSqrTri * (gradVstab * gradVstab) / (Vstab * Vs53 + VSMALL)
			- 2./3. * twoSymm(gradLsq * gradVstab) / Vs53
			+ gradGradLsq / Vs;
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
