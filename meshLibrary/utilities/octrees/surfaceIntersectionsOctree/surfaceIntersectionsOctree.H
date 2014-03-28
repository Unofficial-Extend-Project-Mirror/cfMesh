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
    surfaceIntersectionsOctree

Description
    Helper class to search on triSurf.

SourceFiles
    surfaceIntersectionsOctree.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIntersectionsOctree_H
#define surfaceIntersectionsOctree_H

#include "LongList.H"
#include "boolList.H"
#include "pointIndexHit.H"
#include "surfaceIntersectionsOctreeCube.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class triSurf;

/*---------------------------------------------------------------------------*\
                Class surfaceIntersectionsOctree Declaration
\*---------------------------------------------------------------------------*/

class surfaceIntersectionsOctree
{
    // Private data

        //- Reference to surface to work on
        const triSurf& surface_;

        //- Root cube of the octree structure
        surfaceIntersectionsOctreeCube initialCube_;

        //- check if the line intersects a triangle at a point
        bool checkPointIntersections
        (
            const point&,
            const pointIndexHit&
        ) const;

        //- check if the line intersects a triangle at an edge
        bool checkEdgeIntersections
        (
            const point&,
            const pointIndexHit&
        ) const;

        //- use brutal force search to check if all intersections have been
        //- found by the octree search
        void checkIntersections
        (
            const point&,
            const point&,
            const SLList<pointIndexHit>&
        ) const;

        //- find leaves for cubes
        void findLeavesForCube
        (
            surfaceIntersectionsOctreeCube*,
            LongList<surfaceIntersectionsOctreeCube*>&
        ) const;

        //- Disallow default bitwise copy construct
        surfaceIntersectionsOctree(const surfaceIntersectionsOctree&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceIntersectionsOctree&);

public:

    // Constructors

        //- Construct from surface. Holds reference to surface!
        surfaceIntersectionsOctree(const triSurf&, const short, const direction);

    // Destructor

        ~surfaceIntersectionsOctree();


    // Member Functions
        //- checks if the point is inside or outside the surface
        bool isPointInside(const point&) const;

        //- returns the intersection point nearest to the starting point
        pointIndexHit intersection(const point&, const point&) const;

        //- find intersections of the line with the surface nearest
        //- to the start and end point
        void intersection
        (
            const point&,
            const point&,
            pointIndexHit&,
            pointIndexHit&
        ) const;

        //- helper function which checks if all normal vectors
        //- of the triSurf point outside
        bool normalOk(const label, bool&) const;

        //- check if the edge intersects the surface in a manner
        //- which guarantees that the resulting mesh is of good quality
        bool validIntersection
        (
            const point&,
            const point&,
            SLList<pointIndexHit>&
        ) const;

        //- find a cube containing the vertex
        surfaceIntersectionsOctreeCube* findLeafContainingVertex
        (
            const point&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
