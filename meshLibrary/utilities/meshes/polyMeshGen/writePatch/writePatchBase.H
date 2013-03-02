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
    writePatchBase

Description
    Writes information for a given boundary patch

SourceFiles
    writePatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef writePatchBase_H
#define writePatchBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "label.H"
#include "objectRegistry.H"
#include "typeInfo.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class writePatchBase
{
protected:
    
    // Protected data
        word name_;
    
        word type_;
    
        label nFaces_;
    
        label startFace_;

public:
    
    //- Runtime type information
    TypeName("writePatchBase");

    // Constructors
    
        writePatchBase
        (
            const word&,
            const word&,
            const label,
            const label
        );

        writePatchBase(const word& name, const dictionary&);

    // Destructor

        virtual ~writePatchBase()
        {}

    // Declare run-time constructor selection table
            
        declareRunTimeSelectionTable
        (
            autoPtr,
            writePatchBase,
            dictionary,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );

    // Selectors
        
        static autoPtr<writePatchBase> New
        (
            const word& name,
            const dictionary& dict
        );
        
        static autoPtr<writePatchBase> New
        (
            Istream&
        );
        
    // Construct and return a clone
    
        autoPtr<writePatchBase> clone() const
        {
            return New
            (
                name_,
                this->dict()
            );
        }
        
    //- Write
        virtual void write(Ostream&) const = 0;

        //- Write dictionary
        virtual void writeDict(Ostream&) const = 0;

    // Member functions

        inline const word& patchName() const
        {
            return name_;
        }
        
        inline word& patchName()
        {
            return name_;
        }
        
        inline const word& patchType() const
        {
            return type_;
        }
        
        inline word& patchType()
        {
            return type_;
        }
        
        inline label patchStart() const
        {
            return startFace_;
        }
        
        inline label& patchStart()
        {
            return startFace_;
        }
        
        inline label patchSize() const
        {
            return nFaces_;
        }
        
        inline label& patchSize()
        {
            return nFaces_;
        }
        
        //- Return as dictionary of entries
        virtual dictionary dict() const = 0;
        
    // Member operators
        friend Ostream& operator<<(Ostream&, const writePatchBase&);
        virtual Istream& operator>>(Istream&) = 0;
        //virtual void operator=(const writePatchBase&) = 0;
        //virtual bool operator!=(const writePatchBase&) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
