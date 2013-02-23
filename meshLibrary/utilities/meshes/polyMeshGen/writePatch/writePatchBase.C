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

#include "writePatchBase.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "IOPtrList.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebugWithName
(
	IOPtrList<writePatchBase>,
	"polyBoundaryMesh",
	0
);
    
defineTypeNameAndDebug(writePatchBase, 0);
defineRunTimeSelectionTable(writePatchBase, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<writePatchBase> writePatchBase::New
(
    const word& name,
    const dictionary& dict
)
{
    word type(dict.lookup("type"));
    //- check the type of processor. Allowed types are processor and patch
    //- Other patch types are treated as ordinary patches
    if( type != "processor" )
        type = "patch";
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if( cstrIter == dictionaryConstructorTablePtr_->end() )
    {
        FatalIOErrorIn
        (
            "writePatchBase::New(const word&, const dictionary&)",
            dict
        )   << "Unknown writePatchBase type " << type << nl << nl
            << "Valid writePatchBase types are :" << nl
            << "[default: " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }
    
    return autoPtr<writePatchBase>(cstrIter()(name, dict));
}

autoPtr<writePatchBase> writePatchBase::New
(
    Istream& is
)
{
    word name(is);
    dictionary dict(is);

    return writePatchBase::New(name, dict);
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

writePatchBase::writePatchBase
(
    const word& n,
    const word& t,
    const label nF,
    const label sF
)
:
    name_(n),
    type_(t),
    nFaces_(nF),
    startFace_(sF)
{}
    
writePatchBase::writePatchBase(const word& name, const dictionary& dict)
:
    name_(name),
    type_(),
    nFaces_(),
    startFace_()
{
    word type(dict.lookup("type"));
    type_ = type;
    nFaces_ = readLabel(dict.lookup("nFaces"));
    startFace_ = readLabel(dict.lookup("startFace"));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const writePatchBase& wpb)
{
    wpb.write(os);
    os.check("Ostream& operator<<(Ostream& f, const writePatchBase& wpb");
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
