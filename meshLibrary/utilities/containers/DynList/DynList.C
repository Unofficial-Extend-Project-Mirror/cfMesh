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

#include "DynList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Construct from Istream
template<class T, Foam::label staticSize>
Foam::DynList<T, staticSize>::DynList(Istream& is)
:
    UList<T>(),
    nextFree_(0)
{
    FatalErrorIn
    (
        "template<class T, Foam::label staticSize>"
        "\nFoam::DynList<T, staticSize>::DynList(Istream& is)"
    ) << "Not implemented" << exit(FatalError);
    
    List<T> helper(is);
    
    nextFree_ = helper.size();
    UList<T>::swap(helper);
}


template<class T, Foam::label staticSize>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::DynList<T, staticSize>& DL
)
{
    UList<T> helper(const_cast<T*>(DL.begin()), DL.nextFree_);
    os << helper;

    return os;
}


template<class T, Foam::label staticSize>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::DynList<T, staticSize>& DL
)
{
    FatalErrorIn
    (
        "template<class T, Foam::label staticSize>"
        "\nFoam::Istream& Foam::operator>>"
        "(Foam::Istream& is, Foam::DynList<T, staticSize>& DL)"
    ) << "Not implemented" << exit(FatalError);
    
    is >> static_cast<List<T>&>(DL);
    DL.nextFree_ = DL.List<T>::size();

    return is;
}


// ************************************************************************* //
