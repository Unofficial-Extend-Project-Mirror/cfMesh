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

Description

\*---------------------------------------------------------------------------*/

#include "meshSubset.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline meshSubset::meshSubset()
:
    name_(),
    type_(UNKNOWN),
    data_()
{}

inline meshSubset::meshSubset
(
    const word& name,
    const meshSubset::subsetType_& t
)
:
    name_(name),
    type_(t),
    data_()
{}

template<class ListType>
inline meshSubset::meshSubset
(
    const word& name,
    const meshSubset::subsetType_& type,
    const ListType& elements
)
:
    name_(name),
    type_(type),
    data_()
{
    forAll(elements, i)
        data_.insert(elements[i]);
}

inline meshSubset::meshSubset(Istream& is)
:
    name_(),
    type_(UNKNOWN),
    data_()
{
    is >> *this;
}

inline meshSubset::~meshSubset()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& meshSubset::name() const
{
    return name_;
}

inline label meshSubset::type() const
{
    return type_;
}

template<class ListType>
inline void meshSubset::containedElements(ListType& l) const
{
    l.setSize(data_.size());

    label counter(0);
    forAllConstIter(HashSet<label>, data_, it)
        l[counter++] = it.key();
}

inline void meshSubset::addElement(const label elmt)
{
    data_.insert(elmt);
}

inline void meshSubset::removeElement(const label elmt)
{
    data_.erase(elmt);
}

template<class ListType>
inline void meshSubset::updateSubset(const ListType& newLabels)
{
    HashSet<label> newData;

    forAllConstIter(HashSet<label>, data_, it)
    {
        if( newLabels[it.key()] < 0 )
            continue;

        newData.insert(newLabels[it.key()]);
    }

    data_.clear();
    data_ = newData;
}

inline void meshSubset::updateSubset(const VRWGraph& newLabels)
{
    HashSet<label> newData;

    forAllConstIter(HashSet<label>, data_, it)
    {
        forAllRow(newLabels, it.key(), i)
            newData.insert(newLabels(it.key(), i));
    }

    data_.clear();
    data_ = newData;
}

inline bool meshSubset::contains(const label elmt) const
{
    return data_.found(elmt);
}

inline void meshSubset::operator=(const meshSubset& ms)
{
    name_ = ms.name_;
    type_ = ms.type_;
    data_ = ms.data_;
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const meshSubset& sel)
{
    os.check("inline Ostream& operator<<(Ostream&, const meshSubset&)");

    os << sel.name_ << nl << sel.type_ << nl << sel.data_;

    return os;
}

inline Istream& operator>>(Istream& is, meshSubset& sel)
{
    is.check("friend Istream& operator>>(Istream&, meshSubset&)");

    is >> sel.name_ >> sel.type_ >> sel.data_;

    is.check("friend Istream& operator>>(Istream&, meshSubset&)");

    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
