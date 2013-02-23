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

#include "patchRefinement.H"
#include "triSurface.H"
#include "Istream.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchRefinement::patchRefinement()
:
    patchName_(),
    cellSize_(0.0)
{
}

patchRefinement::patchRefinement(const word& n, const scalar& s)
:
    patchName_(n),
    cellSize_(s)
{
}

patchRefinement::patchRefinement(Istream& is)
:
    patchName_(is),
    cellSize_(readScalar(is))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchRefinement::~patchRefinement()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word patchRefinement::patchName() const
{
    return patchName_;
}

label patchRefinement::patchInSurface(const triSurface& ts) const
{
    forAll(ts.patches(), patchI)
        if( ts.patches()[patchI].name() == patchName_ )
            return patchI;

    FatalErrorIn
    (
        "label patchRefinement::patchInSurface(const triSurface& ts) const"
    ) << "Patch " << patchName_
        << " does not exist in surface " << ts.patches()
        << exit(FatalError);

    return -1;
}

scalar patchRefinement::cellSize() const
{
    return cellSize_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void patchRefinement::operator=(const patchRefinement& pr)
{
	patchName_ = pr.patchName_;
	cellSize_ = pr.cellSize_;
}

Istream& operator>>(Istream& is, patchRefinement& pr)
{
    pr.patchName_ = word(is);
    pr.cellSize_ = readScalar(is);
    return is;
}

Ostream& operator<<(Ostream& os, const patchRefinement& pr)
{
    os << pr.patchName_ << " " << pr.cellSize_ << nl;
    return os;
}

bool operator==(const patchRefinement& p1, const patchRefinement& p2)
{
    if( p1.patchName() == p2.patchName() )
        return true;

    return false;
}

bool operator!=(const patchRefinement& p1, const patchRefinement& p2)
{
    return !(p1 == p2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
