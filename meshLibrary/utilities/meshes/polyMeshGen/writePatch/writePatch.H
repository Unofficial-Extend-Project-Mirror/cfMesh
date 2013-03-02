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

Class
    writePatch

Description
    Writes information for a given boundary patch

SourceFiles
    writePatch.C

\*---------------------------------------------------------------------------*/

#ifndef writePatch_H
#define writePatch_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "writePatchBase.H"
#include "typeInfo.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class writePatch
:   public writePatchBase
{

public:
    
    //- Runtime type information
    TypeName("patch");

    // Constructors
    
        writePatch
        (
            const word&,
            const word&,
            const label,
            const label
        );

        writePatch(const word& name, const dictionary&);

    // Return clone of the object
        virtual autoPtr<writePatchBase> clone() const
        {
            return autoPtr<writePatchBase>
            (
                new writePatch
                (
                    patchName(),
                    patchType(),
                    patchSize(),
                    patchStart()
                )
            );
        }
    
        virtual autoPtr<writePatchBase> clone(const writePatch& wp) const
        {
            return autoPtr<writePatchBase>
            (
                new writePatch
                (
                    wp.patchName(),
                    wp.patchType(),
                    wp.patchSize(),
                    wp.patchStart()
                )
            );
        }
        
        //- Return as dictionary of entries
        dictionary dict() const;
    
    // Write
        //- write to Ostream
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;
    
    // Member operators
    
        Ostream& operator<<(Ostream&) const;
    
        Istream& operator>>(Istream&);
    
        void operator=(const writePatch&);
    
        bool operator!=(const writePatch&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
