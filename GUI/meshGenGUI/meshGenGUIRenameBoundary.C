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

bool meshGenGUI::renameBoundaryEntryExist() const
{
    if( meshDict_.found("renameBoundary") )
        return true;
    
    return false;
}

bool meshGenGUI::defaultPatchNameEntryExist() const
{
    if( renameBoundaryEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        
        if( dict.found("defaultName") )
            return true;
    }
    
    return false;
}

word meshGenGUI::defaultPatchName() const
{
    if( defaultPatchNameEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        const word name(dict.lookup("defaultName"));
        return name;
    }
    
    return word();
}

void meshGenGUI::setDefaultPatchName(const word& name)
{
    if( renameBoundaryEntryExist() )
    {
        dictionary dict = meshDict_.subDict("renameBoundary");
        meshDict_.remove("renameBoundary");
        
        if( dict.found("defaultName") )
            dict.remove("defaultName");
        
        dict.add("defaultName", name);
        meshDict_.add("renameBoundary", dict);
    }
    else
    {
        dictionary dict;
        dict.add("defaultName", name);
        meshDict_.add("renameBoundary", dict);
    }
}

void meshGenGUI::removeDefaultPatchName()
{
    if( renameBoundaryEntryExist() )
    {
        dictionary dict = meshDict_.subDict("renameBoundary");
        meshDict_.remove("renameBoundary");
        
        dict.remove("defaultName");
        
        if( dict.toc().size() != 0 )
            meshDict_.add("renameBoundary", dict);
    }
}

bool meshGenGUI::defaultPatchTypeEntryExist() const
{
    if( renameBoundaryEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        
        if( dict.found("defaultType") )
            return true;
    }
    
    return false;
}

word meshGenGUI::defaultPatchType() const
{
    if( defaultPatchNameEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        const word type(dict.lookup("defaultType"));
        return type;
    }
    
    return word();
}

void meshGenGUI::setDefaultPatchType(const word& type)
{
    if( renameBoundaryEntryExist() )
    {
        dictionary dict = meshDict_.subDict("renameBoundary");
        meshDict_.remove("renameBoundary");
        
        if( dict.found("defaultType") )
            dict.remove("defaultType");
        
        dict.add("defaultType", type);
        meshDict_.add("renameBoundary", dict);
    }
    else
    {
        dictionary dict;
        dict.add("defaultType", type);
        meshDict_.add("renameBoundary", dict);
    }
}

void meshGenGUI::removeDefaultPatchType()
{
    if( renameBoundaryEntryExist() )
    {
        dictionary dict = meshDict_.subDict("renameBoundary");
        meshDict_.remove("renameBoundary");
        
        dict.remove("defaultType");
        
        if( dict.toc().size() != 0 )
            meshDict_.add("renameBoundary", dict);
    }
}

bool meshGenGUI::newPatchNamesEntryExist() const
{
    if( renameBoundaryEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        
        if( dict.found("newPatchNames") )
            return true;
    }
    
    return false;
}

PtrList<entry> meshGenGUI::newPatchNames() const
{
    if( newPatchNamesEntryExist() )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        Istream& is = dict.lookup("newPatchNames");

        PtrList<entry> patchEntries(is);
        return patchEntries;
    }
    
    return PtrList<entry>();
}

void meshGenGUI::addNewPatchName(const word& name, const dictionary& pDict)
{
    if( renameBoundaryEntryExist() )
    {
        dictionary& dict =
            const_cast<dictionary&>(meshDict_.subDict("renameBoundary"));
        
        if( dict.found("newPatchNames") )
        {
            Istream& is = dict.lookup("newPatchNames");
            PtrList<entry> entries(is);
            
            const label s = entries.size();
            entries.setSize(s+1);
            entries.set(s, new dictionaryEntry(name, pDict));
            
            dict.remove("newPatchNames");
            dict.add("newPatchNames", entries);
        }
        else
        {
            PtrList<entry> entries(1);
            entries.set(0, new dictionaryEntry(name, pDict));
            dict.add("newPatchNames", entries);
        }
    }
    else
    {
        PtrList<entry> entries(1);
        entries.set(0, new dictionaryEntry(name, pDict));
        dictionary dict;
        dict.add("newPatchNames", entries);
        meshDict_.add("renameBoundary", dict);
    }
}

void meshGenGUI::removePatchName(const word& name)
{
    if( newPatchNamesEntryExist() )
    {
        dictionary& dict =
            const_cast<dictionary&>(meshDict_.subDict("renameBoundary"));
        
        Istream& is = dict.lookup("newPatchNames");
        PtrList<entry> entries(is);
        
        PtrList<entry> newEntries(entries.size()-1);
        label counter(0);
        
        forAll(entries, i)
        {
            if( entries[i].keyword() == name )
                continue;
            
            newEntries.set
            (
                counter++,
                new dictionaryEntry(entries[i].keyword(), entries[i].dict())
            );
        }
        
        dict.remove("newPatchNames");
        
        if( newEntries.size() != 0 )
            dict.add("newPatchNames", newEntries);
        
        if( dict.toc().size() == 0 )
            meshDict_.remove("renameBoundary");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
