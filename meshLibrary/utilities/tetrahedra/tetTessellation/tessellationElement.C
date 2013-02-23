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

#include "error.H"

#include "tessellationElement.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tessellationElement::tessellationElement()
:
	partTet(),
	influence_(NONE)
{
}

tessellationElement::tessellationElement
(
    const label a,
	const label b,
	const label c,
	const label d
)
:
	partTet(a, b, c, d),
    influence_(NONE)
{
    for(label dir=0;dir<4;++dir)
        neighbours_[dir] = -1;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tessellationElement::~tessellationElement()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar tessellationElement::influencedBy
(
	const LongList<point>& points,
	const point& r
) const
{
	const point crcm = crcmCentre(points);
	
    const scalar d = (magSqr(crcm - points[a()]) - magSqr(crcm - r));

#   ifdef DEBUGtessalation
    if (d > -VSM2 && d < VSM2)
    {
        Info<< "    influenced : new point lies on voronoi circle\n";
    }
#   endif

    return d;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const tessellationElement& elmt)
{
	partTet t(elmt);
    os << t;

	os << " neighbours are ";
    for(label ineigh=0;ineigh<4;++ineigh)
    {
        os << " " << elmt.neighbours_[ineigh];
    }
	os << endl;
	
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
