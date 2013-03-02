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

#include "meshGenGUI.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
bool meshGenGUI::objectRefinementEntryExist() const
{
    return meshDict_.found("objectRefinements");
}

PtrList<entry> meshGenGUI::objectRefinements() const
{
    if( !objectRefinementEntryExist() )
    {
        return PtrList<entry>();
    }
    
    // Read polyPatchList
    Istream& is = meshDict_.lookup("objectRefinements");

    PtrList<entry> objectEntries(is);
    
    return objectEntries;
}

void meshGenGUI::addObjectRefinement(const objectRefinement& object)
{
    if( objectRefinementEntryExist() )
    {
        PtrList<entry> objectEntries(meshDict_.lookup("objectRefinements"));
        const label s = objectEntries.size();
        
        objectEntries.setSize(s+1);
        objectEntries.set
        (
            s,
            new dictionaryEntry(object.name(), object.dict())
        );
        
        meshDict_.remove("objectRefinements");
        meshDict_.add("objectRefinements", objectEntries);
    }
    else
    {
        PtrList<entry> objectEntries(1);
        objectEntries.set
        (
            0,
            new dictionaryEntry
            (
                object.name(),
                object.dict()
            )
        );
        
        meshDict_.add("objectRefinements", objectEntries);
    }
}

void meshGenGUI::removeObjectRefinement(const word& name)
{
    if( !objectRefinementEntryExist() )
        return;
    
    PtrList<entry> objectEntries(meshDict_.lookup("objectRefinements"));
    meshDict_.remove("objectRefinements");
    
    if( objectEntries.size() == 1 )
        return;
    
    PtrList<entry> refs(objectEntries.size()-1);
    label counter(0);
    forAll(objectEntries, i)
        if( objectEntries[i].keyword() != name )
        {
            refs.set
            (
                counter++,
                new dictionaryEntry
                (
                    objectEntries[i].keyword(),
                    objectEntries[i].dict()
                )
            );
        }
        
    meshDict_.add("objectRefinements", refs);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
