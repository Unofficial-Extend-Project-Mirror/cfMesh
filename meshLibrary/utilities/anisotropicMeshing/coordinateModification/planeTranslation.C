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

#include "planeTranslation.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(planeTranslation, 0);
addToRunTimeSelectionTable
(
    coordinateModification,
    planeTranslation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

planeTranslation::planeTranslation()
:
    coordinateModification(),
    origin_(vector::zero),
    normal_(1, 1, 1),
    translationDistance_(0.0),
    scalingFactor_(1.0)
{}

planeTranslation::planeTranslation
(
    const word& name,
    const point& origin,
    const vector& normal,
    const scalar translationDistance,
    const scalar scalingFactor
)
:
    coordinateModification(),
    origin_(origin),
    normal_(normal/mag(normal)),
    translationDistance_(translationDistance),
    scalingFactor_(scalingFactor)
{
    if( scalingFactor_ < SMALL )
    {
        Warning << "Scaling factor for plane " << name << " is less than 0.0 "
                << endl;

        scalingFactor_= 1.0;
    }

    setName(name);
}

planeTranslation::planeTranslation
(
    const word& name,
    const dictionary& dict
)
:
    coordinateModification(name, dict)
{
    this->operator=(dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point planeTranslation::origin() const
{
    return origin_;
}

void planeTranslation::translateAndModifyObject(const vector& disp)
{
    origin_ += disp;

    translationDistance_ /= scalingFactor_;
}

vector planeTranslation::displacement(const point& p) const
{
    const scalar dist = (p - origin_) & normal_;

    const vector translationVec =
        normal_ * translationDistance_ * ((1.0/scalingFactor_) - 1.0);

    const scalar t = dist / translationDistance_;

    const scalar tBnd = Foam::max(0.0, Foam::min(1.0, t));

    return tBnd * translationVec;
}

vector planeTranslation::backwardDisplacement(const point& p) const
{
    const scalar dist = (p - origin_) & normal_;

    const vector translationVec =
        normal_ * translationDistance_ * (scalingFactor_ - 1.0);

    const scalar t = dist / translationDistance_;

    const scalar tBnd = Foam::max(0.0, Foam::min(1.0, t));

    return tBnd * translationVec;
}

bool planeTranslation::combiningPossible() const
{
    return true;
}

void planeTranslation::boundingPlanes(PtrList<plane>& pl) const
{
    if( Foam::mag(scalingFactor_ - 1.0) > VSMALL )
    {
        pl.setSize(2);

        pl.set(0, new plane(origin_, normal_));
        pl.set(1, new plane(origin_ + translationDistance_ * normal_, normal_));
    }
    else
    {
        pl.clear();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary planeTranslation::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("type", type());

    dict.add("origin", origin_);
    dict.add("normal", normal_);
    dict.add("translationDistance", translationDistance_);
    dict.add("scalingFactor", scalingFactor_);

    return dict;
}

void planeTranslation::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " origin: " << origin_
        << " normal: " << normal_
        << " translationDistance: " << translationDistance_
        << " scalingFactor: " << scalingFactor_;
}

void planeTranslation::writeDict(Ostream& os, bool subDict) const
{
    if( subDict )
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if( type() != typeName_() )
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("normal") << normal_ << token::END_STATEMENT << nl;
    os.writeKeyword("translationDistamce") << translationDistance_
                                           << token::END_STATEMENT << nl;
    os.writeKeyword("scalingFactor") << scalingFactor_
                                    << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void planeTranslation::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("origin") )
    {
        dict.lookup("origin") >> origin_;
    }
    else
    {
        FatalErrorIn
        (
            "void planeTranslation::operator=(const dictionary& d)"
        ) << "Entry origin is not specified!" << exit(FatalError);

        origin_ = vector::zero;
    }

    // specify normal
    if( dict.found("normal") )
    {
        dict.lookup("normal") >> normal_;
    }
    else
    {
        FatalErrorIn
        (
            "void planeTranslation::operator=(const dictionary& d)"
        ) << "Entry lengthX is not specified!" << exit(FatalError);

        normal_ = vector(1, 1, 1);
    }

    // specify translation distance
    if( dict.found("translationDistance") )
    {
        translationDistance_ = readScalar(dict.lookup("translationDistance"));
    }
    else
    {
        FatalErrorIn
        (
            "void planeTranslation::operator=(const dictionary& d)"
        ) << "Entry translationDistance is not specified!" << exit(FatalError);

        translationDistance_ = 0.0;
    }

    // specify scaling factor
    if( dict.found("scalingFactor") )
    {
        scalingFactor_ = readScalar(dict.lookup("scalingFactor"));
    }
    else
    {
        WarningIn
        (
            "void planeTranslation::operator=(const dictionary& d)"
        ) << "Entry scalingFactor is not specified!" << endl;

        scalingFactor_ = 1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& planeTranslation::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    write(os);
    return os;
}

Ostream& operator<<(Ostream& os, const planeTranslation& pt)
{
    return pt.operator<<(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
