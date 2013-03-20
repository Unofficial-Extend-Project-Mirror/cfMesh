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

Description

\*---------------------------------------------------------------------------*/

#include "triSurf.H"
#include "boundBox.H"
#include "surfaceIntersectionsOctreeCube.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions
List< List<point> > surfaceIntersectionsOctreeCube::edges() const
{
    const point& min_ = cubeBox_.min();
    const point& max_ = cubeBox_.max();
    pointField cP(8);
    cP[0] = min_;
    cP[1] = point(max_.x(), min_.y(), min_.z());
    cP[2] = point(max_.x(), max_.y(), min_.z());
    cP[3] = point(min_.x(), max_.y(), min_.z());
    cP[4] = point(min_.x(), min_.y(), max_.z());
    cP[5] = point(max_.x(), min_.y(), max_.z());
    cP[6] = max_;
    cP[7] = point(min_.x(), max_.y(), max_.z());

    List< List<point> > edg(12);
    List<point> helper(2);
    //- edges in the x-direction
    helper[0] = cP[0];
    helper[1] = cP[1];
    edg[0] = helper;

    helper[0] = cP[3];
    helper[1] = cP[2];
    edg[1] = helper;

    helper[0] = cP[7];
    helper[1] = cP[6];
    edg[2] = helper;

    helper[0] = cP[4];
    helper[1] = cP[5];
    edg[3] = helper;

    //- edges in the y-direction
    helper[0] = cP[0];
    helper[1] = cP[3];
    edg[4] = helper;

    helper[0] = cP[4];
    helper[1] = cP[7];
    edg[5] = helper;

    helper[0] = cP[5];
    helper[1] = cP[6];
    edg[6] = helper;

    helper[0] = cP[1];
    helper[1] = cP[2];
    edg[7] = helper;

    //- edges in the z-direction
    helper[0] = cP[0];
    helper[1] = cP[4];
    edg[8] = helper;

    helper[0] = cP[1];
    helper[1] = cP[5];
    edg[9] = helper;

    helper[0] = cP[2];
    helper[1] = cP[6];
    edg[10] = helper;

    helper[0] = cP[3];
    helper[1] = cP[7];
    edg[11] = helper;

    //- return edges
    return edg;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
surfaceIntersectionsOctreeCube::surfaceIntersectionsOctreeCube
(
    const triSurf& surface,
    const boundBox& bb,
    const direction& l
)
:
    surface_(surface),
    level_(l),
    cubeBox_(bb),
    containedElements_(),
    subCubesPtr_(NULL),
    cubeType_(surfaceIntersectionsOctreeCube::OUTSIDE)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceIntersectionsOctreeCube::~surfaceIntersectionsOctreeCube()
{
    if( subCubesPtr_ )
    {
        const FixedList<surfaceIntersectionsOctreeCube*, 8>& subCubes = *subCubesPtr_;
        forAll(subCubes, cubeI)
        {
            delete subCubes[cubeI];
        }

        delete subCubesPtr_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FixedList<point, 8> surfaceIntersectionsOctreeCube::vertices() const
{
    FixedList<point, 8> vrt;

    const point& min_ = cubeBox_.min();
    const point& max_ = cubeBox_.max();

    vrt[0] = point(min_.x(), min_.y(), min_.z());
    vrt[1] = point(max_.x(), min_.y(), min_.z());
    vrt[2] = point(min_.x(), max_.y(), min_.z());
    vrt[3] = point(max_.x(), max_.y(), min_.z());
    vrt[4] = point(min_.x(), min_.y(), max_.z());
    vrt[5] = point(max_.x(), min_.y(), max_.z());
    vrt[6] = point(min_.x(), max_.y(), max_.z());
    vrt[7] = point(max_.x(), max_.y(), max_.z());

    return vrt;
}

bool surfaceIntersectionsOctreeCube::isVertexInside(const point& p) const
{
    const point max_ = cubeBox_.max() + point(SMALL, SMALL, SMALL);
    const point min_ = cubeBox_.min() - point(SMALL, SMALL, SMALL);
    if(
            p.x() > max_.x() ||
            p.y() > max_.y() ||
            p.z() > max_.z() ||
            p.x() < min_.x() ||
            p.y() < min_.y() ||
            p.z() < min_.z()
    )
        return false;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
