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
    knuppMetric

Description
    Mesh smoothing using Knupp's metric. The metric is non-zero for inverted
	tets only

SourceFiles
    knuppMetric.C

\*---------------------------------------------------------------------------*/

#ifndef knuppMetric_H
#define knuppMetric_H

#include "simplexSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                           Class knuppMetric Declaration
\*---------------------------------------------------------------------------*/

//- class for mesh untangler
class knuppMetric
:
    public simplexSmoother
{
    // Private data
	//- free vertex
	point& p_;
    
        //- normals of triangles forming the outer hull
        DynList<vector, 64> normals_;
    
        //- centres of triangles forming the outer hull
        DynList<point, 64> centres_;
    
        //- control parameter. It helps the point get inside the feasible region
        scalar beta_;
    
    // Private member functions
        //- evaluate the value of the metric
        scalar evaluateMetric() const;
        scalar evaluateMetricNoBeta() const;
    
        //- evaluate metric gradients
        void evaluateGradients(vector& grad, tensor& gradGrad) const;

    public:

	knuppMetric(partTetMeshSimplex& simplex);

	// Destructor
	~knuppMetric();

	// Member functions
	//- find the best position such that all tets making a simplex
	//- have a positive volume
	void optimizeNodePosition(const scalar tol = 0.001);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
