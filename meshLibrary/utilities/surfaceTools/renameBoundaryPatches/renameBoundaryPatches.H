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
    renameBoundaryPatches

Description
    Renames boundary patches based on user intput and sets correct patch types

SourceFiles
    renameBoundaryPatches.C

\*---------------------------------------------------------------------------*/

#ifndef renameBoundaryPatches_H
#define renameBoundaryPatches_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class IOdictionary;

/*---------------------------------------------------------------------------*\
                           Class renameBoundaryPatches Declaration
\*---------------------------------------------------------------------------*/

class renameBoundaryPatches
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;
    
        //- dictionary containing relevant information
        const IOdictionary& meshDict_;

    // Private member functions
        //- calculate new boundary and replace the existing one
        void calculateNewBoundary();
        
        //- Disallow default bitwise copy construct
        renameBoundaryPatches(const renameBoundaryPatches&);

        //- Disallow default bitwise assignment
        void operator=(const renameBoundaryPatches&);

public:

    // Constructors

        //- Construct from mesh and IOdictionary
        renameBoundaryPatches
        (
            polyMeshGen& mesh,
            const IOdictionary& meshDict
        );

    // Destructor

        ~renameBoundaryPatches();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
