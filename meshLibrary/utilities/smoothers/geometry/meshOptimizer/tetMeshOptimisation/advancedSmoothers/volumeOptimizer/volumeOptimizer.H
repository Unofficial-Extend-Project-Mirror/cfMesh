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
    volumeOptimizer

Description
    Mesh smoothing without any topological changes. The vertex is put into
    location which minimises the objective function.

SourceFiles
    volumeOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef volumeOptimizer_H
#define volumeOptimizer_H

#include "simplexSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                    Class volumeOptimizer Declaration
\*---------------------------------------------------------------------------*/

//- class for volume optimizer
class volumeOptimizer
:
    public simplexSmoother
{
    // Private member functions
        //- evaluate functional
        scalar evaluateFunc() const;

        //- find appropriate value of K
        scalar evaluateStabilisationFactor() const;

        //- evaluate gradients of the functional
        void evaluateGradientsExact(vector&, tensor&) const;

        //- optimize position using a divide and conquer algorithm
        scalar optimiseDivideAndConquer(const scalar tol);

        //- optimise using the steepest descent
        scalar optimiseSteepestDescent(const scalar tol);

    public:

    // Constructor
        //- construct from partTetMeshSimplex
        volumeOptimizer(partTetMeshSimplex& simplex);

    // Destructor
        ~volumeOptimizer();

    // Member functions
        //- find the best position for the node
        void optimizeNodePosition(const scalar tol = 0.001);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
