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

#include "error.H"
#include "trianglePlaneIntersections.H"
#include "helperFunctions.H"

namespace Foam
{

void trianglePlaneIntersections::calculateIntersections()
{
    const pointField& points = ts_.points();

    const labelledTri& lt = ts_[ltri_];
    direction i(0);
    forAll(lt, pI)
    {
        vector v = points[lt[pI]] - pp_;
        v /= mag(v);

        if( mag(v & n_) < SMALL )
        {
            intersectedPoints_[pI] = true;
            i++;
        }
    }

    // hold the number of intersected edges
    direction j(0);

    if( i == 2 )
    {
        influence_ = TWO_VERTICES;
    }
    else
    {
        forAll(lt, eI)
        {
            const direction next = (eI+1) % 3;

            if( !intersectedPoints_[eI] && !intersectedPoints_[next] )
            {
                vector v(points[lt[next]] - points[lt[eI]]);

                if( mag(n_ & (v/mag(v))) > SMALL )
                {
                    const scalar t((n_ & (pp_ - points[lt[eI]])) / (n_ & v));

                    if( (t > -SMALL) && (t < (1.0+SMALL)) )
                    {
                        edgePoints_[eI] = points[lt[eI]] + t * v;
                        intersectedEdges_[eI] = true;
                        j++;
                    }
                }
            }
        }
    }

    if( i == 1 && j == 0 )
    {
        influence_ = ONE_VERTEX;
    }
    else if( i == 1 && j == 1 )
    {
        influence_ = VERTEX_AND_EDGE;
    }
    else if( i == 0 && j == 2 )
    {
        influence_ = TWO_EDGES;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Constructor
trianglePlaneIntersections::trianglePlaneIntersections
(
    const vector& n,
    const point& pp,
    const triSurface& ts,
    const label lI
)
    :
    n_(n),
    pp_(pp),
    ts_(ts),
    ltri_(lI),
    intersectedPoints_(3, false),
    edgePoints_(3, vector::zero),
    intersectedEdges_(3, false),
    influence_(NONE)
{
    calculateIntersections();
}

trianglePlaneIntersections::~trianglePlaneIntersections()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Member functions
direction trianglePlaneIntersections::triInfluence() const
{
    return influence_;
}

const boolList& trianglePlaneIntersections::intersectedPoints() const
{
    return intersectedPoints_;
}

const pointField& trianglePlaneIntersections::edgePoints() const
{
    return edgePoints_;
}

const boolList& trianglePlaneIntersections::intersectedEdges() const
{
    return intersectedEdges_;
}

label trianglePlaneIntersections::determineRotation
(
    const labelList& intersectedNeighbours,
    const pointField& intersectionPoints,
    const vector& n,
    const face& f,
    const pointField& polyPoints,
    const face& newF,
    const direction npI,
    const DynList<point>& newPoints
) const
{
//     Info << "Intersection points " << intersectionPoints << endl;
//     Info << "Normal " << n << endl;

    bool in0 = help::pointInsideFace(intersectionPoints[0], f, n, polyPoints);
    bool in1 = help::pointInsideFace(intersectionPoints[1], f, n, polyPoints);

//     Info << "in0 " << in0 << " in1 " << in1 << endl;
    
    if( in0 && !in1 )
    {
        return intersectedNeighbours[0];
    }
    else if( !in0 && in1 )
    {
        return intersectedNeighbours[1];
    }
    
    vector ev(intersectionPoints[1] - intersectionPoints[0]);
    ev /= mag(ev);

    const direction prev = (npI-1) >= 0?(npI-1):(newF.size()-1);
    vector pv = newPoints[newF[npI]] - newPoints[newF[prev]];
    pv /= mag(pv);

    vector ref = n ^ pv;
    ref /= mag(ref);

//     Info << "ev " << ev << endl;
//     Info << "ref " << ref << endl;

    if( (ev & ref) >= 0.0 )
    {
        return intersectedNeighbours[1];
    }

    return intersectedNeighbours[0];
}

label trianglePlaneIntersections::determineRotation
(
    const face& f,
    const vector& n,
    const pointField& fp,
    const face& newF,
    const direction npI,
    const DynList<point>& newPoints,
    direction& pointI
) const
{
    const labelList& fe = ts_.faceEdges()[ltri_];
    const labelListList& ef = ts_.edgeFaces();

    labelList neis(2);
    pointField eP(2);
    labelList interPoint(2, -1);
    direction counter(0);

    const pointField& points = ts_.points();

    forAll(intersectedPoints_, pI)
        if( intersectedPoints_[pI] )
        {
            if( pointInsideFace(f, n, fp, pI) )
            {
                pointI = pI;
                return -1;
            }
            
            neis[counter] = -1;
            interPoint[counter] = pI;
            eP[counter++] =
                points[ts_[ltri_][pI]];
        }
    
    forAll(intersectedEdges_, eJ)
        if( intersectedEdges_[eJ] )
        {
            label neighbour = ef[fe[eJ]][0];
            if( neighbour == ltri_ )
                neighbour = ef[fe[eJ]][1];
            
            if( edgePointInsideFace(f, n, fp, eJ) )
                return neighbour;
            
            neis[counter] = neighbour;
            eP[counter++] = edgePoints_[eJ];
        }
    
    //- rotation is still not determined
    if( counter == 2 )
    {
        vector v(eP[1] - eP[0]);
        v /= mag(v);
        
        const direction prev = ((npI-1)>=0?(npI-1):(newF.size()-1));
        vector e(newPoints[newF[npI]] - newPoints[newF[prev]]);
        e /= mag(e);
        
        if( ((e ^ v) & n) >= 0.0 )
        {
            if( interPoint[1] != -1 )
                pointI = interPoint[1];
            return neis[1];
        }
        else
        {
            if( interPoint[0] != -1 )
                pointI = interPoint[0];
            return neis[0];
        }
    }

    return -1;
}

bool trianglePlaneIntersections::edgePointInsideFace
(
    const face& f,
    const vector& n,
    const pointField& fp,
    const direction eI
) const
{
    return help::pointInsideFace(edgePoints_[eI], f, n, fp);
}

bool trianglePlaneIntersections::pointInsideFace
(
    const face& f,
    const vector& n,
    const pointField& fp,
    const direction pI
) const
{
    const pointField& points = ts_.points();
    return help::pointInsideFace(points[ts_[ltri_][pI]], f, n, fp);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
