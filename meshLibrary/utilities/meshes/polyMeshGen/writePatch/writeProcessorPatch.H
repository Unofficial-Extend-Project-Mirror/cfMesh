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
    writeProcessorPatch

Description
    Writes information for a given processor patch

SourceFiles
    writeProcessorPatch.C

\*---------------------------------------------------------------------------*/

#ifndef writeProcessorPatch_H
#define writeProcessorPatch_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "writePatchBase.H"
#include "typeInfo.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class writeProcessorPatch
:   public writePatchBase 
{
    // private data
    
        label myProcNo_;
    
        label neighbProcNo_;

public:
    
    //- Runtime type information
    TypeName("processor");

    // Construct from components
        writeProcessorPatch
        (
            const word& name,
            const word& type,
            const label nFaces,
            const label startFace,
            const label myProcNo,
            const label neighbProcNo
        );

        writeProcessorPatch(const word& name, const dictionary&);

    // Construct and return a clone
        virtual autoPtr<writePatchBase> clone() const
        {
            return autoPtr<writePatchBase>
			(
				new writeProcessorPatch
				(
					patchName(),
                    patchType(),
                    patchSize(),
                    patchStart(),
                    myProcNo_,
                    neighbProcNo_
				)
			);
        }

        virtual autoPtr<writePatchBase> clone
		(
			const writeProcessorPatch& pp
		) const
        {
            return autoPtr<writePatchBase>
			(
				new writeProcessorPatch
				(
					pp.patchName(),
                    pp.patchType(),
                    pp.patchSize(),
                    pp.patchStart(),
                    pp.myProcNo_,
                    pp.neighbProcNo_
				)
			);
        }
        
    // Member functions
        //- return the owner processor
        inline label myProcNo() const
        {
            return myProcNo_;
        }
        
        //- return the neighbour processor
        inline label neiProcNo() const
        {
            return neighbProcNo_;
        }
        
        //- check if the processor is the owner of the interface
        inline bool owner() const
        {
            if( myProcNo_ < neighbProcNo_ )
                return true;
            
            return false;
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
    
        void operator=(const writeProcessorPatch&);
    
        bool operator!=(const writeProcessorPatch&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
