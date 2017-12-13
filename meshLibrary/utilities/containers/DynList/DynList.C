/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DynList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, Foam::label staticSize>
Foam::Module::DynList<T, staticSize>::DynList(Istream&)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    NotImplemented;
}


template<class T, Foam::label staticSize>
Foam::Ostream& Foam::Module::operator<<
(
    Foam::Ostream& os,
    const Foam::Module::DynList<T, staticSize>& DL
)
{
    UList<T> helper(DL.dataPtr_, DL.nextFree_);
    os << helper;

    return os;
}


template<class T, Foam::label staticSize>
Foam::Istream& Foam::Module::operator>>
(
    Foam::Istream& is,
    Foam::Module::DynList<T, staticSize>& DL
)
{
    NotImplemented;

    UList<T> helper(DL.dataPtr_, DL.nextFree_);
    //is >> static_cast<List<T>&>(DL);
    is >> helper;
    DL.nextFree_ = helper.size();

    return is;
}


// ************************************************************************* //
