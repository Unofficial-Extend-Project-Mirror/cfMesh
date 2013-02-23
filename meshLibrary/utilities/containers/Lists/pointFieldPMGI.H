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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constructors

inline pointFieldPMG::pointFieldPMG(const IOobject& io)
:
	regIOobject(io),
	pointField(readStream(typeName)),
	nElmts_(pointField::size())
{
	close();
}

inline pointFieldPMG::pointFieldPMG(const IOobject& io, const label s)
:
	regIOobject(io),
	pointField(s),
	nElmts_(s)
{
}

inline pointFieldPMG::pointFieldPMG(const IOobject& io, const pointField& pts)
:
	regIOobject(io),
	pointField(pts),
	nElmts_(pts.size())
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor

inline pointFieldPMG::~pointFieldPMG()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
			
inline label pointFieldPMG::size() const
{
	return nElmts_;
}

inline void pointFieldPMG::setSize(const label nElmts)
{
	if( nElmts >= pointField::size() )
	{
		Info << "Resizing points!" << endl;
		pointField::setSize(label(1.5*nElmts)+1);
	}

	nElmts_ = nElmts;
}

inline void pointFieldPMG::clear()
{
	nElmts_ = 0;
}

inline void pointFieldPMG::append(const point& p)
{
	const label i = nElmts_;
	setSize(i+1);
	this->operator[](i) = p;
}

inline point& pointFieldPMG::newElmt(const label pI)
{
	setSize(pI+1);
	return this->operator[](pI);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool pointFieldPMG::writeData(Ostream& os) const
{
	return (os << *this).good();
}

inline void pointFieldPMG::operator=(const pointField& pts)
{
	setSize(pts.size());
	forAll(pts, pI)
		this->operator[](pI) = pts[pI];
}
	
inline Ostream& operator<<(Ostream& os, const pointFieldPMG& pts)
{
	SubList<point> p(pts, pts.nElmts_, 0);
	
	os << p;
	return os;
}

inline Istream& operator>>(Istream& is, pointFieldPMG& pts)
{
	pointField& points = static_cast<pointField&>(pts);
	is >> points;
	pts.nElmts_  = points.size();
	
	return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
