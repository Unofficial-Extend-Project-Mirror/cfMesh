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
    coneScaling

Description
    Scales coordinates inside the cone and translates points
    outside of the cone

SourceFiles
    coneScaling.C

\*---------------------------------------------------------------------------*/

#ifndef coneScaling_H
#define coneScaling_H

#include "coordinateModification.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class coneScaling Declaration
\*---------------------------------------------------------------------------*/

class coneScaling
:
    public coordinateModification
{
    // Private data
        //- starting point and the radius
        point p0_;
        scalar r0_;

        //- end point and the radius
        point p1_;
        scalar r1_;

        //- scaling factor in the radial direction
        scalar radialScaling_;

        //- scaling factor in the axial direction
        scalar axialScaling_;

        //- axial vector
        vector axialVec_;

        //- squared length of the axial vector
        scalar axialVecLengthSq_;

    // Private member functions
        //- calculate axial vector
        void calculateAxialVector();

public:

    //- Runtime type information
    TypeName("cone");


    // Constructors

        //- Null construct
        coneScaling();

        //- Construct from components
        coneScaling
        (
            const word& name,
            const point& p0,
            const scalar& r0,
            const point& p1,
            const scalar& r1,
            const scalar radialScaling = 1.0,
            const scalar axialScaling = 1.0
        );

        //- Construct from dictionary
        coneScaling(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const coneScaling& cs
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new coneScaling
                (
                    cs.name(),
                    cs.p0_,
                    cs.r0_,
                    cs.p1_,
                    cs.r1_,
                    cs.radialScaling_,
                    cs.axialScaling_
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
