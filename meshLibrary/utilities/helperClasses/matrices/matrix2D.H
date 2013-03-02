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
    Implementation of 2 x 2 matrix

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef matrix2D_H
#define matrix2D_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                           Class matrix2D Declaration
\*---------------------------------------------------------------------------*/

class matrix2D
: public FixedList<FixedList<scalar, 2>, 2>
{
    // Private members
        mutable scalar det_;
        mutable bool calculatedDet_;
    
    // Private member functions
        //- calculate matrix determinant
        inline void calculateDeterminant() const;
    
public:

    // Constructors
        //- Null constructor
        explicit inline matrix2D();

        //- Copy constructor
        inline matrix2D(const matrix2D&);

    // Destructor
        inline ~matrix2D();

    // Member functions
        //- return matrix determinant
        inline scalar determinant() const;

        //- return inverse matrix
        inline matrix2D inverse() const;

        //- find the solution of the system with the given rhs
        inline FixedList<scalar, 2> solve
        (
            const FixedList<scalar, 2>& source
        ) const;

        //- return the first component of the solution vector
        inline scalar solveFirst(const FixedList<scalar, 2>& source) const;
        
        //- return the second component of the solution vector
        inline scalar solveSecond(const FixedList<scalar, 2>& source) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "matrix2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
