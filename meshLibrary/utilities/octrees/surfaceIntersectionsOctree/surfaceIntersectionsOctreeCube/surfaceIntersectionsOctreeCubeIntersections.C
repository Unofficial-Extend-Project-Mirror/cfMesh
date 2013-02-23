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

#include "triSurface.H"
#include "surfaceIntersectionsOctreeCube.H"
#include "helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool surfaceIntersectionsOctreeCube::intersectsTriangle(const label tI) const
{
    const pointField& points = surface_.points();
    const labelledTri& ltri = surface_.localFaces()[tI];

    point max_(-GREAT, -GREAT, -GREAT);
    point min_(GREAT, GREAT, GREAT);
    forAll(ltri, pI)
    {
        max_ = Foam::max(points[ltri[pI]], max_);
        min_ = Foam::min(points[ltri[pI]], min_);
    }
    
    min_ -= point(SMALL, SMALL, SMALL);
    max_ += point(SMALL, SMALL, SMALL);

    return cubeBox_.overlaps(boundBox(min_, max_));
}

bool surfaceIntersectionsOctreeCube::intersectsLine(const point& s, const point& e) const
{
    //- check if the cube contains start point or end point
    const point min = cubeBox_.min() - vector(SMALL,SMALL,SMALL);
    const point max = cubeBox_.max() + vector(SMALL,SMALL,SMALL);

    //- check for intersections of line with the cube faces
    const vector v(e - s);
    if( mag(v.x()) > SMALL )
    {
        if(
            ((s.x() < min.x()) && (e.x() < min.x())) ||
            ((s.x() > max.x()) && (e.x() > max.x()))
        )
            return false;

        {
            //- x-min face
            const vector n(-1, 0, 0);
            const scalar t = (n & (min - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL))
            )
                if(
                    (i.y() > min.y()) && (i.y() < max.y()) &&
                    (i.z() > min.z()) && (i.z() < max.z())
                )
                    return true;
        }
        {
            //- x-max face
            const vector n(1, 0, 0);
            const scalar t = (n & (max - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL))
            )
                if(
                    (i.y() > min.y()) && (i.y() < max.y()) &&
                    (i.z() > min.z()) && (i.z() < max.z())
                )
                    return true;
        }
    }

    if( mag(v.y()) > SMALL)
    {
        if(
            ((s.y() < min.y()) && (e.y() < min.y())) ||
            ((s.y() > max.y()) && (e.y() > max.y()))
        )
            return false;

        {
            //- y-min face
            const vector n(0, -1, 0);
            const scalar t = (n & (min - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL))
            )
                if(
                    (i.x() > min.x()) && (i.x() < max.x()) &&
                    (i.z() > min.z()) && (i.z() < max.z())
                )
                    return true;
        }
        {
            //- y-max face
            const vector n(0, 1, 0);
            const scalar t = (n & (max - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL))
            )
                if(
                    (i.x() > min.x()) && (i.x() < max.x()) &&
                    (i.z() > min.z()) && (i.z() < max.z())
                )
                    return true;
        }
    }
    if( mag(v.z()) > SMALL )
    {
        if(
            ((s.z() < min.z()) && (e.z() < min.z())) ||
            ((s.z() > max.z()) && (e.z() > max.z()))
        )
            return false;

        {
            //- z-min face
            const vector n(0, 0, -1);
            const scalar t = (n & (min - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL)) &&
                (i.x() > min.x()) && (i.x() < max.x()) &&
                (i.y() > min.y()) && (i.y() < max.y())
            )
                return true;
        }
        {
            //- z-min face
            const vector n(0, 0, 1);
            const scalar t = (n & (max - s)) / (n & v);
            const point i = s + t * v;
            if(
                (t > -SMALL) && (t < (1.0+SMALL)) &&
                (i.x() > min.x()) && (i.x() < max.x()) &&
                (i.y() > min.y()) && (i.y() < max.y())
            )
                return true;
        }
    }

    if(
        (s.x() > min.x()) && (s.x() < max.x()) &&
        (s.y() > min.y()) && (s.y() < max.y()) &&
        (s.z() > min.z()) && (s.z() < max.z())
    )
        return true;

    return false;
}

void surfaceIntersectionsOctreeCube::intersectionCandidates
(
    SLList<pointIndexHit>& cp,
    const point& s,
    const point& e
) const
{
    if( intersectsLine(s, e) )
    {
        # ifdef DEBUGSearch
        Info << "Cube " << cubeBox_ << " contains elements "
            << containedElements_ << endl;
        # endif
        point intersection;
        //- add elements within the cube
        for(SLList<label>::const_iterator eIter = containedElements_.begin();
            eIter != containedElements_.end();
            ++eIter
        )
        {
            if(
                help::triLineIntersection
                (
                    surface_,
                    eIter(),
                    s,
                    e,
                    intersection
                )
            )
                cp.append(pointIndexHit(true,intersection, eIter()));
        }
        
        //- check subCubes if there are any
        if( subCubesPtr_ )
        {
            const FixedList<surfaceIntersectionsOctreeCube*, 8>& subCubes_ = *subCubesPtr_;
            forAll(subCubes_, cubeI)
                subCubes_[cubeI]->intersectionCandidates(cp, s, e);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
