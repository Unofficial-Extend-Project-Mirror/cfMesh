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
    matrix3D

Description
    Implementation od 3 x 3 matrix

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
inline void matrix3D::calculateDeterminant() const
{
	if( calculatedDet_ )
		return;
	
	const FixedList<FixedList<scalar, 3>, 3>& mat = *this;
	det_ =
		mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
		mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
		mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
	calculatedDet_ = true;
}
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
inline matrix3D::matrix3D()
:
	det_(),
	calculatedDet_(false)
{
}

inline matrix3D::matrix3D(const matrix3D& m)
:
	FixedList<FixedList<scalar, 3>, 3>(m),
	det_(),
	calculatedDet_(false)
{
}

inline matrix3D::~matrix3D()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar matrix3D::determinant() const
{
	calculateDeterminant();
	
	return det_;
}

inline matrix3D matrix3D::inverse() const
{
	calculateDeterminant();
	const FixedList<FixedList<scalar, 3>, 3>& mat = *this;
	matrix3D imat;
	imat[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / det_;
	imat[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) / det_;
	imat[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / det_; 
	imat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) / det_;
	imat[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) / det_;
	imat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) / det_;
	imat[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / det_;
	imat[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) / det_;
	imat[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / det_;
	
	return imat;
}

inline FixedList<scalar, 3> matrix3D::solve
(
	const FixedList<scalar, 3>& source
) const
{
	FixedList<scalar, 3> result;
	result[0] = solveFirst(source);
	result[1] = solveSecond(source);
	result[2] = solveThird(source);
	
	return result;
}

inline scalar matrix3D::solveFirst(const FixedList<scalar, 3>& source) const
{
	calculateDeterminant();
	
	const FixedList<FixedList<scalar, 3>, 3>& mat = *this;
	
	return 
	(
		(mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * source[0] +
		(mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * source[1] +
		(mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * source[2]
	) / det_;
}
		
inline scalar matrix3D::solveSecond(const FixedList<scalar, 3>& source) const
{
	calculateDeterminant();
	
	const FixedList<FixedList<scalar, 3>, 3>& mat = *this;
	
	return
	(
		(mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * source[0] +
		(mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * source[1] +
		(mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * source[2]
	) / det_;
}
		
inline scalar matrix3D::solveThird(const FixedList<scalar, 3>& source) const
{
	calculateDeterminant();
	
	const FixedList<FixedList<scalar, 3>, 3>& mat = *this;
	
	return
	(
		(mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * source[0] +
		(mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * source[1] +
		(mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * source[2]
	) / det_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
