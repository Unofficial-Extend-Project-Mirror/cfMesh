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

#include "sphereScaling.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sphereScaling, 0);
addToRunTimeSelectionTable(coordinateModification, sphereScaling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphereScaling::sphereScaling()
:
    coordinateModification(),
    centre_(vector::zero),
    radius_(0.0),
    radialScaling_(1.0)
{}

sphereScaling::sphereScaling
(
    const word& name,
    const point& centre,
    const scalar radius,
    const scalar radialScaling
)
:
    coordinateModification(),
    centre_(centre),
    radius_(radius),
    radialScaling_(radialScaling)
{
    setName(name);
}

sphereScaling::sphereScaling
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

point sphereScaling::origin() const
{
    return centre_;
}

void sphereScaling::translateAndModifyObject(const vector& disp)
{
    centre_ += disp;

    radius_ /= radialScaling_;

    if( radius_ < VSMALL )
    {
        FatalErrorIn
        (
            "void sphereScaling::translateAndModifyObject(const vector&)"
        ) << "Radius of the sphere is zero " << *this << exit(FatalError);
    }
}

vector sphereScaling::displacement(const point& p) const
{
    vector disp(vector::zero);

    const vector rVec = p - centre_;
    const scalar r = mag(rVec);

    if( r > VSMALL )
    {
        const scalar tRadial = r / radius_;
        const scalar tRadialBnd = Foam::max(0.0, Foam::min(1.0, tRadial));

        const scalar rScale = (1.0/radialScaling_) - 1.0;
        const vector dispRadial = radius_ * (rVec / r) * rScale;

        //- calculate displacements in the radial direction
        disp += tRadialBnd * dispRadial;
    }

    return disp;
}

vector sphereScaling::backwardDisplacement(const point& p) const
{
    vector disp(vector::zero);

    const vector rVec = p - centre_;
    const scalar r = mag(rVec);

    if( r > VSMALL )
    {
        const scalar tRadial = r / radius_;
        const scalar tRadialBnd = Foam::max(0.0, Foam::min(1.0, tRadial));

        const scalar rScale = radialScaling_ - 1.0;
        const vector dispRadial = radius_ * (rVec / r) * rScale;

        //- scale the distance in the radial direction
        disp += tRadialBnd * dispRadial;
    }

    return disp;
}

bool sphereScaling::combiningPossible() const
{
    if( Foam::mag(radialScaling_ - 1.0) > VSMALL )
        return false;

    return true;
}

void sphereScaling::boundingPlanes(PtrList<plane>& pl) const
{
    pl.clear();

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary sphereScaling::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("type", type());

    dict.add("centre", centre_);
    dict.add("radius", radius_);

    dict.add("radialScaling", radialScaling_);

    return dict;
}

void sphereScaling::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " centre: " << centre_
        << " radius:  " << radius_
        << " radialScaling:  " << radialScaling_
        << endl;
}

void sphereScaling::writeDict(Ostream& os, bool subDict) const
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

    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius") << radius_ << token::END_STATEMENT << nl;
    os.writeKeyword("radialScaling") << radialScaling_
                                     << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void sphereScaling::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("centre") )
    {
        centre_ = vector(dict.lookup("centre"));
    }
    else
    {
        FatalErrorIn
        (
            "void sphereScaling::operator=(const dictionary& d)"
        ) << "Entry centre is not specified!" << exit(FatalError);
        centre_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius") )
    {
        radius_ = readScalar(dict.lookup("radius"));
    }
    else
    {
        FatalErrorIn
        (
            "void sphereScaling::operator=(const dictionary& d)"
        ) << "Entry radius is not specified!" << exit(FatalError);

        radius_ = 0.0;
    }

    // specify radialScaling
    if( dict.found("radialScaling") )
    {
        radialScaling_ = readScalar(dict.lookup("radialScaling"));
    }
    else
    {
        WarningIn
        (
            "void sphereScaling::operator=(const dictionary& d)"
        ) << "Entry radialScaling is not specified!" << endl;

        radialScaling_ = 1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& sphereScaling::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    write(os);
    return os;
}

Ostream& operator<<(Ostream& os, const sphereScaling& ss)
{
    return ss.operator<<(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
