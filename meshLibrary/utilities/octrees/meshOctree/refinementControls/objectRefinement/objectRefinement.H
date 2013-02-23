/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    objectRefinement

Description
    Base class for coordinate systems.
    All systems are defined by an origin and a coordinate rotation.

SourceFiles
    objectRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef objectRefinement_H
#define objectRefinement_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class boundBox;

/*---------------------------------------------------------------------------*\
                     Class objectRefinement Declaration
\*---------------------------------------------------------------------------*/

class objectRefinement
{
    // Private data

        //- Name of object refinement type
        word name_;

        //- cell size for this object
        scalar cellSize_;

public:

    // Runtime type information
        TypeName("objectRefinement");


    // Constructors

        //- Null construct
        objectRefinement();

        //- Construct from dictionary
        objectRefinement(const word& name, const dictionary&);

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectRefinement,
            dictionary,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<objectRefinement> New
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~objectRefinement();


    // Member Functions

        //- check if a boundBox intersects or is inside the object
        virtual bool intersectsObject(const boundBox&) const = 0;

      // Access

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- set the name
        void setName(const word& name)
        {
            name_ = name;
        }

        //- set cell size for this object
        void setCellSize(const scalar s)
        {
            cellSize_ = s;
        }

        //- return cell size for this object
        scalar cellSize() const
        {
            return cellSize_;
        }

        //- Return as dictionary of entries
        virtual dictionary dict(bool ignoreType = false) const = 0;

      // Write

        //- Write
        virtual void write(Ostream&) const = 0;

        //- Write dictionary
        virtual void writeDict(Ostream&, bool subDict = true) const = 0;

    // Member Operators

        //- assign from dictionary
        virtual void operator=(const dictionary&) = 0;

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const objectRefinement&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
