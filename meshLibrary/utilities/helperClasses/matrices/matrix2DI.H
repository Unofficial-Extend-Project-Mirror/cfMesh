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
    matrix2D

Description
    Implementation od 2 x 2 matrix

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
inline void matrix2D::calculateDeterminant() const
{
	if( calculatedDet_ )
		return;
	
	const FixedList<FixedList<scalar, 2>, 2>& mat = *this;
	det_ = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		
	calculatedDet_ = true;
}
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
inline matrix2D::matrix2D()
:
	det_(),
	calculatedDet_(false)
{
}

inline matrix2D::matrix2D(const matrix2D& m)
:
	FixedList<FixedList<scalar, 2>, 2>(m),
	det_(),
	calculatedDet_(false)
{
}

inline matrix2D::~matrix2D()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar matrix2D::determinant() const
{
	calculateDeterminant();
	
	return det_;
}

inline matrix2D matrix2D::inverse() const
{
	calculateDeterminant();
	const FixedList<FixedList<scalar, 2>, 2>& mat = *this;
	matrix2D imat;
	imat[0][0] = mat[1][1] / det_;
	imat[0][1] = -1.0 * mat[0][1] / det_;
	imat[1][0] = -1.0 * mat[1][0] / det_; 
	imat[1][1] = mat[0][0] / det_;
	
	return matrix2D(imat);
}

inline FixedList<scalar, 2> matrix2D::solve
(
	const FixedList<scalar, 2>& source
) const
{
	FixedList<scalar, 2> result;
	result[0] = solveFirst(source);
	result[1] = solveSecond(source);
	
	return result;
}

inline scalar matrix2D::solveFirst(const FixedList<scalar, 2>& source) const
{
	calculateDeterminant();
	
	const FixedList<FixedList<scalar, 2>, 2>& mat = *this;
	
	return 
	(
		mat[1][1] * source[0] -
		mat[0][1] * source[1]
	) / det_;
}
		
inline scalar matrix2D::solveSecond(const FixedList<scalar, 2>& source) const
{
	calculateDeterminant();
	
	const FixedList<FixedList<scalar, 2>, 2>& mat = *this;
	
	return
	(
		-1.0 * mat[1][0] * source[0] +
		mat[0][0] * source[1]
	) / det_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
