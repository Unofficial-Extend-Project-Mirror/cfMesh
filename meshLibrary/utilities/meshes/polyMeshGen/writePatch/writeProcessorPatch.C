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

#include "writeProcessorPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(writeProcessorPatch, 0);
addToRunTimeSelectionTable(writePatchBase, writeProcessorPatch, dictionary);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

writeProcessorPatch::writeProcessorPatch
(
    const word& name,
    const word& type,
    const label nFaces,
    const label startFace,
    const label myProcNo,
    const label neighbProcNo
)
:
    writePatchBase(name, type, nFaces, startFace),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo)
{}
    
writeProcessorPatch::writeProcessorPatch
(
    const word& name,
    const dictionary& dict
)
:
    writePatchBase(name, dict),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo")))
{
}   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
dictionary writeProcessorPatch::dict() const
{
    dictionary dict;

    dict.add("type", type_);

    dict.add("nFaces", nFaces_);
    dict.add("startFace", startFace_);
    dict.add("myProcNo", myProcNo_);
    dict.add("neighbProcNo", neighbProcNo_);

    return dict;
}

void writeProcessorPatch::write(Ostream& os) const
{
    this->operator<<(os);
}

void writeProcessorPatch::writeDict(Ostream& os) const
{
    
}

Ostream& writeProcessorPatch::operator<<(Ostream& os) const
{
    os  << patchName() << nl << token::BEGIN_BLOCK << nl
        << "    type         " << patchType() << token::END_STATEMENT << nl
        << "    nFaces       " << patchSize() << token::END_STATEMENT << nl
        << "    startFace    " << patchStart() << token::END_STATEMENT << nl
        << "    myProcNo     " << myProcNo_ << token::END_STATEMENT << nl
        << "    neighbProcNo " << neighbProcNo_
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;

    return os;
}

Istream& writeProcessorPatch::operator>>(Istream& is)
{
    token t;
    is >> name_ >> t;
    is >> t >> type_ >> t;
    is >> t >> nFaces_ >> t;
    is >> t >> startFace_ >> t;
    is >> t >> myProcNo_ >> t;
    is >> t >> neighbProcNo_ >> t;
    is >> t;

    return is;
}

void writeProcessorPatch::operator=(const writeProcessorPatch& wp)
{
    name_ = wp.name_;
    type_ = wp.type_;
    nFaces_ = wp.nFaces_;
    startFace_ = wp.startFace_;
    myProcNo_ = wp.myProcNo_;
    neighbProcNo_ = wp.neighbProcNo_;
}

bool writeProcessorPatch::operator!=(const writeProcessorPatch& wp) const
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
    else if(
        (myProcNo_ != wp.myProcNo_) || (neighbProcNo_ != wp.neighbProcNo_)
    )
    {
        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
