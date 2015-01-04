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
    sphereScaling

Description
    Scales distances inside the sphere and translates points outside

SourceFiles
    sphereScaling.C

\*---------------------------------------------------------------------------*/

#ifndef sphereScaling_H
#define sphereScaling_H

#include "coordinateModification.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sphereScaling Declaration
\*---------------------------------------------------------------------------*/

class sphereScaling
:
    public coordinateModification
{
    // Private data
        //- centre of the sphere
        point centre_;

        //- radius of the sphere
        scalar radius_;

        //- scaling factor in the radial direction
        scalar radialScaling_;

public:

    //- Runtime type information
    TypeName("sphereScaling");


    // Constructors

        //- Null construct
        sphereScaling();

        //- Construct from name, centre, radius and scaling factor
        sphereScaling
        (
            const word& name,
            const point& centre,
            const scalar radius,
            const scalar radialScaling = 1.0
        );

        //- Construct from dictionary
        sphereScaling(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const sphereScaling& ss
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new sphereScaling
                (
                    ss.name(),
                    ss.centre_,
                    ss.radius_,
                    ss.radialScaling_
                )
            );
        }

    // Member Functions

        //- return the centre of the box
        virtual point origin() const;

        //- translate the object to the modified coordinates
        //- this is needed for backward transformation
        virtual void translateAndModifyObject(const vector&);

        //- calculate the displacement vector for box scaling
        virtual vector displacement(const point&) const;

        //- calculate the displacement vector for box scaling
        virtual vector backwardDisplacement(const point&) const;

        //- can this modification object be combined with other ones
        virtual bool combiningPossible() const;

        //- return that "bounding planes" of the scaling region for
        //- the given object
        virtual void boundingPlanes(PtrList<plane>&) const;

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
