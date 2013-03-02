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
    surfaceIntersectionsOctreeCube

Description
    A cube for octree needed for evaluating surface intersections 

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef surfaceIntersectionsOctreeCube_H
#define surfaceIntersectionsOctreeCube_H

#include "SLList.H"
#include "boolList.H"
#include "pointIndexHit.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class surfaceIntersectionsOctreeCube Declaration
\*---------------------------------------------------------------------------*/

class surfaceIntersectionsOctreeCube
{
    // Private data

        //- Reference to surface to work on
        const triSurface& surface_;

        //- cube level in the octree structure
        direction level_;

        //- Bound box of the cube
        boundBox cubeBox_;

        //- labels of elements contained within this surfaceIntersectionsOctreeCube
        SLList<label> containedElements_;

        //- list of cube created from the current one
        //- if the size of the list is not zero than the number of contained
        //- elements is zero
        FixedList<surfaceIntersectionsOctreeCube*, 8>* subCubesPtr_;

        //- is the cube outside, data or outside of the domain
        direction cubeType_;

        //- returns cube's edges
        List< List<point> > edges() const;

        //- Disallow default bitwise copy construct
        surfaceIntersectionsOctreeCube(const surfaceIntersectionsOctreeCube&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceIntersectionsOctreeCube&);

public:
    // Constructors

        //- Construct from surface. Holds reference to surface!
        surfaceIntersectionsOctreeCube(const triSurface&, const boundBox&, const direction&);

    // Destructor

        ~surfaceIntersectionsOctreeCube();

    // Enumerators
        enum cubeTypes
        {
            OUTSIDE = 0,
            DATA = 1,
            INSIDE = 2
        };

    // Member Functions
        //- subdivide the octree cube
        void refineTree(const short, const direction);

        //- is cube a leaf
        inline bool isLeaf() const;
        
        //- return and set cube type
        direction cubeType() const
        {
            return cubeType_;
        };
        void setCubeType(const direction t)
        {
            cubeType_ = t;
        }

        //- vertices
        FixedList<point, 8> vertices() const;
    
        //- sub cubes
        inline FixedList<surfaceIntersectionsOctreeCube*, 8>* subCubes() const;

        //- check if the triangle intersects the cube
        bool intersectsTriangle(const label) const;

        //- return boundBox
        inline const boundBox& bb() const;

        //- add a triangle into the cube
        inline void append(const label);

        //- check if the cube intersects a line
        bool intersectsLine(const point&, const point&) const;

        //- return contained elements
        inline const SLList<label>& containedElements() const;

        //- add labels of triangles contained within the cube
        void intersectionCandidates
        (
            SLList<pointIndexHit>&,
            const point&,
            const point&
        ) const;

        //- is a vertex inside the cube
        bool isVertexInside(const point&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceIntersectionsOctreeCubeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
