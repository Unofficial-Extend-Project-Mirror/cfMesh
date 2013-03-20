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

Class
    trianglePlaneIntersections

Description
    Calculates intersections between a plane a triangle

SourceFiles
    trianglePlaneIntersections.C

\*---------------------------------------------------------------------------*/

#ifndef trianglePlaneIntersections_H
#define trianglePlaneIntersections_H

#include "triSurf.H"
#include "DynList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class trianglePlaneIntersections
{
    // Private members
        const vector& n_;
        const point& pp_;

        const triSurf& ts_;
        const label ltri_;

        boolList intersectedPoints_;
        pointField edgePoints_;
        boolList intersectedEdges_;

        direction influence_;
    // Private member functions

        void calculateIntersections();



public:

    trianglePlaneIntersections
    (
        const vector&,
        const point&,
        const triSurf&,
        const label
    );

    ~trianglePlaneIntersections();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    enum typeOfIntersection
    {
        NONE = 0,
        ONE_VERTEX = 2,
        VERTEX_AND_EDGE = 4,
        TWO_VERTICES = 8,
        TWO_EDGES = 16
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

        direction triInfluence() const;

        const boolList& intersectedPoints() const;

        const pointField& edgePoints() const;

        const boolList& intersectedEdges() const;

        label determineRotation
        (
            const labelList&,
            const pointField&,
            const vector&,
            const face&,
            const pointField&,
            const face&,
            const direction,
            const DynList<point>&
        ) const;

        label determineRotation
        (
            const face&,
            const vector&,
            const pointField&,
            const face&,
            const direction,
            const DynList<point>&,
            direction&
        ) const;

        bool edgePointInsideFace
        (
            const face&,
            const vector&,
            const pointField&,
            const direction
        ) const;

        bool pointInsideFace
        (
            const face&,
            const vector&,
            const pointField&,
            const direction
        ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
