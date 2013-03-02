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
    patchRefinement

Description
    A class used for octree refinement. Refinement of boxes intersected
    by a given patch to the prescribed size.

SourceFiles
    patchRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef patchRefinement_H
#define patchRefinement_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurface;
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class patchRefinement Declaration
\*---------------------------------------------------------------------------*/

class patchRefinement
{
    // Private members

        word patchName_;

        scalar cellSize_;

public:

    // Constructors
        //- Null contructor
        patchRefinement();

        //- Construct from word and desired cell size
        patchRefinement(const word&, const scalar&);

        //- contruct from stream
        patchRefinement(Istream&);

    // Destructor

        ~patchRefinement();


    // Member Functions
        word patchName() const;

        label patchInSurface(const triSurface&) const;

        scalar cellSize() const;
  
    // Operators

        void operator=(const patchRefinement&);
        friend Istream& operator>>(Istream&, patchRefinement&);
        friend Ostream& operator<<(Ostream&, const patchRefinement&);
        friend bool operator==(const patchRefinement&, const patchRefinement&);
        friend bool operator!=(const patchRefinement&, const patchRefinement&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
