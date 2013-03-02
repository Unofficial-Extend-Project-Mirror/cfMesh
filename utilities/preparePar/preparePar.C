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

Application
    Prepares the case for a parallel mesh generation run

Description
    - creates processor* directories which contain data for processors

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"

#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
    
    objectRegistry registry(runTime);
    
    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            registry.time().system(),
            registry,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary decomposeParDict
    (
        IOobject
        (
            "decomposeParDict",
            registry.time().system(),
            registry,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const label nProcessors
    (
        readLabel(decomposeParDict.lookup("numberOfSubdomains"))
    );
    
    for(label procI=0;procI<nProcessors;++procI)
    {
        fileName file("processor");
        std::ostringstream ss;
        ss << procI;
        file += ss.str();
        Info << "Creating " << file << endl;
        
        //- create a directory for processor data
        mkDir(runTime.path()/file);
        
        //- copy the contents of the const directory into processor*
        cp(registry.path()/"constant", runTime.path()/file);
    }
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
