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

#include "writePatch.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(writePatch, 0);
addToRunTimeSelectionTable(writePatchBase, writePatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

writePatch::writePatch
(
    const word& n,
    const word& t,
    const label nF,
    const label sF
)
:
    writePatchBase(n, t, nF, sF)
{}
    
writePatch::writePatch(const word& name, const dictionary& dict)
:
    writePatchBase(name, dict)
{
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary writePatch::dict() const
{
    dictionary dict;

	dict.add("type", type_);

    dict.add("nFaces", nFaces_);
    dict.add("startFace", startFace_);

    return dict;
}

void writePatch::write(Ostream& os) const
{
    this->operator<<(os);
}

void writePatch::writeDict(Ostream& os) const
{
    
}

Ostream& writePatch::operator<<(Ostream& os) const
{
    os  << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type_ << token::END_STATEMENT << nl
        << "    nFaces " << nFaces_ << token::END_STATEMENT << nl
        << "    startFace " << startFace_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;

    return os;
}

Istream& writePatch::operator>>(Istream& is)
{
	token t;
	is >> name_ >> t;
	is >> t >> type_ >> t;
	is >> t >> nFaces_ >> t;
	is >> t >> startFace_ >> t;
	is >> t;

    return is;
}

void writePatch::operator=(const writePatch& wp)
{
    name_ = wp.name_;
    type_ = wp.type_;
    nFaces_ = wp.nFaces_;
    startFace_ = wp.startFace_;
}

bool writePatch::operator!=(const writePatch& wp) const
{
    if( name_ != wp.name_ )
    {
        return true;
    }
    else if( type_ != wp.name_ )
    {
        return true;
    }
    else if( (nFaces_ != wp.nFaces_) || (startFace_ != wp.startFace_) )
    {
        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
