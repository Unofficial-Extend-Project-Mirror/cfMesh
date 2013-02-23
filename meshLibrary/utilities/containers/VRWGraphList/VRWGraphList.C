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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "VRWGraphList.H"
#include "token.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::VRWGraphList& DL
)
{
	os << DL.size() << nl << token::BEGIN_LIST;
	
	for(register label i=0;i<DL.size();++i)
	{
		os << nl << DL[i];
	}
	
	os << nl << token::END_LIST;
    return os;
}

/*
template<class T, Foam::label width>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::VRWGraphList<T, width>& DL
)
{
	label size;
	T e;
	is >> size;
	DL.setSize(size);
	for(IndexType i=0;i<size;++i)
	{
		is >> e;
		DL[i] = e;
	}

    return is;
}
*/

// ************************************************************************* //
