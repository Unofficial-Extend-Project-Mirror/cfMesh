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
    delaunayElement

Description
    Class for Delaunay tetrahedra

SourceFiles
    

\*---------------------------------------------------------------------------*/

#ifndef tessellationElement_H
#define tessellationElement_H

#include "LongList.H"
#include "point.H"
#include "triFace.H"
#include "partTet.H"
#include "direction.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class elmt Declaration
\*---------------------------------------------------------------------------*/

class tessellationElement
: public partTet
{
    // Private data
	
		//- neighbours of the given element
		label neighbours_[4];
	
public:

    // Constructors

        //- Construct null
		tessellationElement();

		//- Construct from nodes
        tessellationElement
		(
			const label a,
			const label b,
			const label c,
			const label d
		);

    // Destructor

        ~tessellationElement();

	// Enumerations and public members

		//- types of influence a node can have on an tessellationElement
        enum typesOfInfluence
        {
            NONE = 0,
            GOOD = 1,
            BOUND = 2,
            RMVISIT = 4,
            BVISIT = 8
        };

        direction influence_;

    // Member functions
		//- swap node with another one
		inline void setNode(const direction, const label);

		//- Return index of a neighbour in the given direction
        inline label neighbour(const direction) const;
		
		//- set neighbour in a given direction
		inline void setNeighbour(const direction, const label);
		
		//- return face in the given direction
		inline triFace face(const direction) const;

        //- Return the influence the point r has on this tessellationElement.
        //  Result is the difference between the sqr of the radius of
        //  the circumcircle and the sqr distance between the point and the
        //  centre of the circumcircle.  The result is positive if the point
        //  is inside the circle (i.e. influences it).
        scalar influencedBy
		(
			const LongList<point>& points,
			const point& r
		) const;

		// Ostream Operator
        friend Ostream& operator<<(Ostream&, const tessellationElement&);
		
		inline void operator=(const tessellationElement&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tessellationElementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
