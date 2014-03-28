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

#include "surfaceIntersectionsOctree.H"
#include "triSurf.H"
#include "boundBox.H"
#include "SLList.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
surfaceIntersectionsOctree::surfaceIntersectionsOctree
(
    const triSurf& surface,
    const short maxN,
    const direction maxL
)
:
    surface_(surface),
    initialCube_
    (
        surface_,
        boundBox
        (
            min(surface_.points()) - point(0.001, 0.001, 0.001),
            max(surface_.points()) + point(0.001, 0.001, 0.001)
        ),
        0
    )
{
    Info << "Constructing octree" << endl;
    //- all triangles are within the initial cube
    forAll(surface_, faceI)
        initialCube_.append(faceI);

    //- find the smallest edge. Size of the smallest surfaceIntersectionsOctreeCube should not be
    //- smaller than the smallest edge
    const edgeLongList& edges = surface_.edges();
    const pointField& points = surface_.points();

    scalar dist(GREAT);
    forAll(edges, eI)
        if( edges[eI].mag(points) < dist )
            dist = edges[eI].mag(points);

    (const_cast<triSurf&>(surface_)).clearAddressing();

    label n = label(initialCube_.bb().mag() / dist);

    direction k(0);

    while( (pow(2, label(k)) < n) && (k < maxL) )
    {
        ++k;
    }

    //- refine the tree such that there are not more than maxNInCube triangle
    //- in a single octree cube
    initialCube_.refineTree(maxN, k);

    LongList<surfaceIntersectionsOctreeCube*> leaves;
    findLeavesForCube(&initialCube_, leaves);

    forAll(leaves, leafI)
    {
        surfaceIntersectionsOctreeCube& oc = *leaves[leafI];

        if( oc.containedElements().size() > 0 )
        {
            oc.setCubeType(surfaceIntersectionsOctreeCube::DATA);
        }
        else
        {
            bool inside(true);
            const FixedList<point, 8> vrt = oc.vertices();
            forAll(vrt, vrtI)
                if( !isPointInside(vrt[vrtI]) )
                {
                    inside = false;
                    break;
                }

            if( inside )
                oc.setCubeType(surfaceIntersectionsOctreeCube::INSIDE);
        }
    }

    Info << "Finished constructing octree" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceIntersectionsOctree::~surfaceIntersectionsOctree()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void surfaceIntersectionsOctree::findLeavesForCube
(
    surfaceIntersectionsOctreeCube* oc,
    LongList<surfaceIntersectionsOctreeCube*>& lvs
) const
{
    if( oc->isLeaf() )
    {
        //- cube is a leaf
        lvs.append(oc);
    }
    else
    {
        //- cube is not a leaf
        const FixedList<surfaceIntersectionsOctreeCube*, 8>& sc = *oc->subCubes();

        //- check subCubes
        forAll(sc, scI)
            findLeavesForCube(sc[scI], lvs);
    }
}

surfaceIntersectionsOctreeCube*
surfaceIntersectionsOctree::findLeafContainingVertex(const point& p) const
{
    # ifdef OCTREE_DEBUG
    Info << "Finding leaf for vertex " << p << endl;
    # endif

    surfaceIntersectionsOctreeCube* oc =
        const_cast<surfaceIntersectionsOctreeCube*>(&initialCube_);

    if( !oc->isVertexInside(p) )
    {
        # ifdef OCTREE_DEBUG
        Info << "Vertex " << p << " is not in the initial cube" << endl;
        # endif
        return NULL;
    }

    bool finished(false);

    do
    {
        if( !oc->isLeaf() )
        {
            //- find a subCube containing the vertex;
            const FixedList<surfaceIntersectionsOctreeCube*, 8>& sc =
                *oc->subCubes();

            bool found(false);

            forAll(sc, scI)
                if( sc[scI]->isVertexInside(p) )
                {
                    oc = sc[scI];
                    # ifdef OCTREE_DEBUG
                    Info << "New searching cube is " << *oc << endl;
                    # endif
                    found = true;
                    break;
                }

            if( !found )
                FatalErrorIn
                (
                    "surfaceIntersectionsOctreeCube* meshOctree::"
                    "findLeafContainingVertex(const point& p) const"
                ) << "Vertex is not found in any subCubes!!"
                    << abort(FatalError);
        }
        else
        {
            finished = true;
        }
    } while( !finished );

    return oc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
