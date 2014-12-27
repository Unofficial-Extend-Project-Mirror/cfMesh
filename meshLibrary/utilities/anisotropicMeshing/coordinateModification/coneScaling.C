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

#include "coneScaling.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(coneScaling, 0);
addToRunTimeSelectionTable(coordinateModification, coneScaling, dictionary);

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void coneScaling::calculateAxialVector()
{
    axialVec_ = p1_ - p0_;

    axialVecLengthSq_ = magSqr(axialVec_);

    if( axialVecLengthSq_ < VSMALL )
    {
        FatalErrorIn
        (
            "void coneScaling::calculateAxialVector()"
        ) << "The length of the cone in the axial direction is zero"
          << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coneScaling::coneScaling()
:
    coordinateModification(),
    p0_(vector::zero),
    r0_(0.0),
    p1_(vector::zero),
    r1_(0.0),
    radialScaling_(1.0),
    axialScaling_(1.0),
    axialVec_(vector::zero),
    axialVecLengthSq_(0.0)
{}

coneScaling::coneScaling
(
    const word& name,
    const point& p0,
    const scalar& r0,
    const point& p1,
    const scalar& r1,
    const scalar radialScaling,
    const scalar axialScaling
)
:
    coordinateModification(),
    p0_(p0),
    r0_(r0),
    p1_(p1),
    r1_(r1),
    radialScaling_(radialScaling),
    axialScaling_(axialScaling),
    axialVec_(vector::zero),
    axialVecLengthSq_(0.0)
{
    setName(name);
    calculateAxialVector();
}

coneScaling::coneScaling
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

point coneScaling::origin() const
{
    return p0_;
}

void coneScaling::translateAndModifyObject(const vector& disp)
{
    p0_ += disp;
    p1_ += disp;

    r0_ /= radialScaling_;
    r1_ /= radialScaling_;

    calculateAxialVector();
}

vector coneScaling::displacement(const point& p) const
{
    vector disp(vector::zero);

    //- calculate point position in the axial direction
    const scalar t = ((p - p0_) & axialVec_) / (axialVecLengthSq_ + VSMALL);

    //- calculate point position in the radial direction
    const vector rVec = p - (p0_ + t * axialVec_);
    const scalar r = mag(rVec);

    //- calculate the radius of the cone at the axial position
    const scalar tBnd = Foam::max(0.0, Foam::min(t, 1.0));
    const scalar rCone = (1.0 - tBnd) * r0_ + tBnd * r1_;

    //- calculate displacement in the axial direction
    if( t > 1.0 )
    {
        //- translate in the axial direction
        disp += axialVec_ * ((1.0/axialScaling_) - 1.0);
    }
    else if( t > 0.0 )
    {
        //- scale the distance in the axial direction
        disp += t * axialVec_ * ((1.0/axialScaling_) - 1.0);
    }

    //- calculate displacement in the radial direction
    if( r > rCone )
    {
        //- translate in the radial direction
        disp += rVec * ((1.0/radialScaling_) - 1.0);
    }
    else
    {
        //- scale the distance in the radial direction
        disp += (r / rCone) * rVec * ((1.0/radialScaling_) - 1.0);
    }

    return disp;
}

vector coneScaling::backwardDisplacement(const point& p) const
{
    vector disp(vector::zero);

    //- calculate point position in the axial direction
    const scalar t = ((p - p0_) & axialVec_) / (axialVecLengthSq_ + VSMALL);

    //- calculate point position in the radial direction
    const vector rVec = p - (p0_ + t * axialVec_);
    const scalar r = mag(rVec);

    //- calculate the radius of the cone at the axial position
    const scalar tBnd = Foam::max(0.0, Foam::min(t, 1.0));
    const scalar rCone = (1.0 - tBnd) * r0_ + tBnd * r1_;

    //- calculate displacement in the axial direction
    const scalar aScale = 1.0 - axialScaling_;
    if( t > 1.0 )
    {
        //- translate in the axial direction
        disp -= axialVec_ * aScale;
    }
    else if( t > 0.0 )
    {
        //- scale the distance in the axial direction
        disp -= t * axialVec_ * aScale;
    }

    //- calculate displacement in the radial direction
    const scalar rScale = 1.0 - radialScaling_;
    if( r > rCone )
    {
        //- translate in the radial direction
        disp -= rVec * rScale;
    }
    else
    {
        //- scale the distance in the radial direction
        disp -= (r / rCone) * rVec * rScale;
    }

    return disp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary coneScaling::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("type", type());

    dict.add("p0", p0_);
    dict.add("radius0", r0_);
    dict.add("p1", p1_);
    dict.add("radius1", r1_);

    dict.add("radialScaling", radialScaling_);
    dict.add("axialScaling", axialScaling_);

    return dict;
}

void coneScaling::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " p0: " << p0_
        << " radius0: " << r0_
        << " p1: " << p1_
        << " radius1: " << r1_
        << " radialScaling:  " << radialScaling_
        << " axialScaling:  " << axialScaling_
        << endl;
}

void coneScaling::writeDict(Ostream& os, bool subDict) const
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

    os.writeKeyword("p0") << p0_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius0") << r0_ << token::END_STATEMENT << nl;
    os.writeKeyword("p1") << p1_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius1") << r1_ << token::END_STATEMENT << nl;

    os.writeKeyword("radialScaling") << radialScaling_
                                   << token::END_STATEMENT << nl;
    os.writeKeyword("axialScaling") << axialScaling_
                                  << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void coneScaling::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // specify p0
    if( dict.found("p0") )
    {
        p0_ = vector(dict.lookup("p0"));
    }
    else
    {
        FatalErrorIn
        (
            "void coneScaling::operator=(const dictionary& d)"
        ) << "Entry p0 is not specified!" << exit(FatalError);

        p0_ = vector::zero;
    }

    // specify r0
    if( dict.found("radius0") )
    {
        r0_ = readScalar(dict.lookup("radius0"));
    }
    else
    {
        FatalErrorIn
        (
            "void coneScaling::operator=(const dictionary& d)"
        ) << "Entry radius0 is not specified!" << exit(FatalError);

        r0_ = 0.0;
    }

    // specify p1
    if( dict.found("p1") )
    {
        p1_ = vector(dict.lookup("p1"));
    }
    else
    {
        FatalErrorIn
        (
            "void coneScaling::operator=(const dictionary& d)"
        ) << "Entry p1 is not specified!" << exit(FatalError);

        p1_ = vector::zero;
    }

    // specify radius1
    if( dict.found("radius1") )
    {
        r1_ = readScalar(dict.lookup("radius1"));
    }
    else
    {
        r1_ = 1.0;
    }

    // specify radialScaling
    if( dict.found("radialScaling") )
    {
        radialScaling_ = readScalar(dict.lookup("radialScaling"));
    }
    else
    {
        radialScaling_ = 1.0;
    }

    // specify axialScaling
    if( dict.found("axialScaling") )
    {
        axialScaling_ = readScalar(dict.lookup("axialScaling"));
    }
    else
    {
        axialScaling_ = 1.0;
    }

    calculateAxialVector();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& coneScaling::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    write(os);
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
