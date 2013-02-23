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
    surfaceOptimizer

Description
    Performs optimisation of a central vertex in a simplex

SourceFiles
    surfaceOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceOptimizer_H
#define surfaceOptimizer_H

#include "point.H"
#include "triFace.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceOptimizer Declaration
\*---------------------------------------------------------------------------*/

class surfaceOptimizer
{
    // Private data
		//- reference to the simplex points
        DynList<point>& pts_;
    
        //- reference to the triangles forming a simplex
        const DynList<triFace>& trias_;

        //- Disallow default bitwise copy construct
        surfaceOptimizer(const surfaceOptimizer&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceOptimizer&);

public:

    // Constructors

        //- Construct from transformed points and triangles forming a simplex
        surfaceOptimizer
        (
			DynList<point>& pts,
            const DynList<triFace>& trias
        );


    // Destructor

        ~surfaceOptimizer();

    // Member Functions
		
		//- optimizes position of a central point in the simplex
		point optimizePoint(const scalar tol = 0.1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
