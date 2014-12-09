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

\*---------------------------------------------------------------------------*/

#include "coordinateModifier.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coordinateModifier::coordinateModifier(const dictionary& geomModDict)
:
    modificationDict_(geomModDict),
    modifiers_(),
    backwardModifiers_()
{
    const wordList modifiers = modificationDict_.toc();

    //- setup modification
    modifiers_.setSize(modifiers.size());
    backwardModifiers_.setSize(modifiers.size());
    forAll(modifiers, modI)
    {
        const word& mName = modifiers[modI];
        const dictionary& modDict = modificationDict_.subDict(mName);
        modifiers_.set(modI, coordinateModification::New(mName, modDict));

        backwardModifiers_.set
        (
            modI,
            coordinateModification::New(mName, modDict)
        );
    }

    //- setup backward modification
    forAll(backwardModifiers_, modI)
    {
        vector disp(vector::zero);
        const point pOrigin = backwardModifiers_[modI].origin();

        forAll(modifiers_, i)
            disp += modifiers_[i].displacement(pOrigin);

        backwardModifiers_[modI].translateAndModifyObject(disp);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coordinateModifier::~coordinateModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

point coordinateModifier::modifiedPoint(const point& p) const
{
    point pNew = p;

    forAll(modifiers_, modI)
        pNew += modifiers_[modI].displacement(p);

    return pNew;
}

point coordinateModifier::backwardModifiedPoint(const point& p) const
{
    point pNew = p;

    forAll(backwardModifiers_, modI)
        pNew += backwardModifiers_[modI].backwardDisplacement(p);

    return pNew;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
