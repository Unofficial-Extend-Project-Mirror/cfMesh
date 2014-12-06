/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    planeTranslation

Description
    Applies translation to points on the positive side of the plane.

SourceFiles
    planeTranslation.C

\*---------------------------------------------------------------------------*/

#ifndef planeTranslation_H
#define planeTranslation_H

#include "coordinateModification.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class planeTranslation Declaration
\*---------------------------------------------------------------------------*/

class planeTranslation
:
    public coordinateModification
{
    // Private data
        //- origin of the plane
        point origin_;

        //- normal vector
        vector normal_;

        //- translation distance
        scalar translationDistance_;

        //- scaling factor
        scalar scalingFactor_;

public:

    //- Runtime type information
    TypeName("planeTranslation");


    // Constructors

        //- Null construct
        planeTranslation();

        //- Construct from name, origin, normal, translation distance and scale
        planeTranslation
        (
            const word& name,
            const point& original,
            const vector& normal,
            const scalar translationDistance,
            const scalar scalingFactor
        );

        //- Construct from dictionary
        planeTranslation(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const planeTranslation& pt
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new planeTranslation
                (
                    pt.name(),
                    pt.origin_,
                    pt.normal_,
                    pt.translationDistance_,
                    pt.scalingFactor_
                )
            );
        }

    // Member Functions

        //- return the origin of the plane
        virtual point origin() const;

        //- translate the object to the modified coordinates
        //- this is needed for backward transformation
        virtual void translateAndModifyObject(const vector&);

        //- calculate the displacement vector for plane translation
        virtual vector displacement(const point&) const;

        //- calculate the displacement vector for plane translation
        virtual vector backwardDisplacement(const point&) const;

    //- Return as dictionary of entries
        dictionary dict(bool ignoreType = false) const;

      // Write

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&, bool subDict = true) const;

    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);

    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
