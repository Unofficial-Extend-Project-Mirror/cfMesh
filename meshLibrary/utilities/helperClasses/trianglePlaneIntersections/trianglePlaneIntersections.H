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

        label influence_;
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

        label triInfluence() const;

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
            const label,
            const DynList<point>&
        ) const;

        label determineRotation
        (
            const face&,
            const vector&,
            const pointField&,
            const face&,
            const label,
            const DynList<point>&,
            label&
        ) const;

        bool edgePointInsideFace
        (
            const face&,
            const vector&,
            const pointField&,
            const label
        ) const;

        bool pointInsideFace
        (
            const face&,
            const vector&,
            const pointField&,
            const label
        ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
