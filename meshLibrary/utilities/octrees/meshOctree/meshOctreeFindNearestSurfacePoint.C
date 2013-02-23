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

#include "meshOctree.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctree::findNearestSurfacePoint
(
    point& nearest,
    scalar& distSq,
    label& region,
    const point& p
) const
{
    region = -1;

    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if( cLabel < 0 )
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 1.5 * leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    //- find nearest surface vertex to the point p
    bool found;
	label iterationI(0);
    DynList<const meshOctreeCube*, 256> neighbours;

    do
    {
        found = false;
        boundBox bb(p - sizeVec, p + sizeVec);
        distSq = Foam::sqr(sizeVec.x());

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

        //- find nearest projection
        forAll(neighbours, neiI)
        {
            if( !neighbours[neiI]->hasContainedElements() )
                continue;

            const VRWGraph& ct =
                neighbours[neiI]->slotPtr()->containedTriangles_;
            const constRow el =
                ct[neighbours[neiI]->containedElements()];
            forAll(el, tI)
            {
                const point p0 =
                    help::nearestPointOnTheTriangle(el[tI], surface_, p);

                const scalar dSq = magSqr(p0 - p);
                if( dSq < distSq )
                {
                    distSq = dSq;
                    nearest = p0;
                    region = surface_[el[tI]].region();
                    found = true;
                }
            }
        }

        if( !found )
            sizeVec *= 2.0;

    } while( !found && (iterationI++ < 5) );

    if( (!found || (region < 0)) && !Pstream::parRun() )
        Warning << "Could not find a boundary region for vertex " << p << endl;
}

void meshOctree::findNearestSurfacePointInRegion
(
    point& nearest,
    scalar& distSq,
    const label region,
    const point& p
) const
{
    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if( cLabel < 0 )
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 1.5 * leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    //- find nearest surface vertex to the point p
    bool found;
	label iterationI(0);
    DynList<const meshOctreeCube*, 256> neighbours;

    do
    {
        found = false;
        boundBox bb(p - sizeVec, p + sizeVec);
        distSq = Foam::sqr(sizeVec.x());

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

        //- find nearest projection
        forAll(neighbours, neiI)
        {
            if( !neighbours[neiI]->hasContainedElements() )
                continue;

            const VRWGraph& ct =
                neighbours[neiI]->slotPtr()->containedTriangles_;
            const constRow el =
                ct[neighbours[neiI]->containedElements()];
            forAll(el, tI)
            {
                if( surface_[el[tI]].region() != region )
                    continue;

                const point p0 =
                    help::nearestPointOnTheTriangle(el[tI], surface_, p);

                const scalar dSq = magSqr(p0 - p);
                if( dSq < distSq )
                {
                    distSq = dSq;
                    nearest = p0;
                    found = true;
                }
            }
        }

        if( !found )
            sizeVec *= 2.0;

    } while( !found && (iterationI++ < 5) );

    if( (!found || (region < 0)) && !Pstream::parRun() )
        Warning << "Could not find a boundary region for vertex " << p << endl;
}

bool meshOctree::findNearestEdgePoint
(
    const point& p,
    const DynList<label>& regions,
    point& edgePoint,
    scalar& distSq
) const
{
    //- find the estimate for the searching range
    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if( cLabel < 0 )
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 1.5 * leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    DynList<const meshOctreeCube*, 256> neighbours;

    const pointField& sp = surface_.localPoints();
    const edgeList& edges = surface_.edges();
    const labelListList& edgeFaces = surface_.edgeFaces();

    edgePoint = p;
    bool foundAnEdge(false);
    label iterationI(0);

    do
    {
        boundBox bb(p - sizeVec, p + sizeVec);
        distSq = Foam::sqr(sizeVec.x());

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

        forAll(neighbours, neiI)
        {
            if( !neighbours[neiI]->hasContainedEdges() )
                continue;

            const VRWGraph& containedEdges =
                neighbours[neiI]->slotPtr()->containedEdges_;
            const constRow ce =
                containedEdges[neighbours[neiI]->containedEdges()];

            forAll(ce, eI)
            {
                //- find if the edge is in correct patches
                bool correctPatches(true);
                const labelList& ef = edgeFaces[ce[eI]];
                forAll(ef, efI)
                {
                    if( !regions.contains(surface_[ef[efI]].region()) )
                    {
                        correctPatches = false;
                        break;
                    }
                }

                if( !correctPatches )
                    continue;

                const point s = sp[edges[ce[eI]].start()];
                const point e = sp[edges[ce[eI]].end()];
                const point np = help::nearestPointOnTheEdgeExact(s, e, p);

                if( magSqr(np - p) < distSq )
                {
                    distSq = magSqr(np - p);
                    edgePoint = np;
                    foundAnEdge = true;
                }
            }
        }

        if( !foundAnEdge )
            sizeVec *= 2.0;

    } while( !foundAnEdge && (++iterationI < 5) );

    return foundAnEdge;
}

bool meshOctree::findNearestVertexToTheEdge
(
    const FixedList<point, 2>& edgePoints,
    const FixedList<label, 2>& edgePointRegions,
    point& nearest,
    scalar& distSq
) const
{
    const point c = 0.5 * (edgePoints[0] + edgePoints[1]);
    const scalar dst = mag(edgePoints[0] - edgePoints[1]);
    vector sizeVec(dst, dst, dst);

    boundBox bb(c - 1.5 * sizeVec, c + 1.5 * sizeVec);

    DynList<const meshOctreeCube*, 256> leavesInBox;
    findLeavesContainedInBox(bb, leavesInBox);

    const labelListList& edgeFaces = surface_.edgeFaces();
    const pointField& points = surface_.localPoints();
    const edgeList& surfaceEdges = surface_.edges();

    distSq = VGREAT;

    bool found(false);

    forAll(leavesInBox, leafI)
    {
        if( !leavesInBox[leafI]->hasContainedEdges() )
            continue;

        const VRWGraph& containedEdges =
            leavesInBox[leafI]->slotPtr()->containedEdges_;
        const constRow edges =
            containedEdges[leavesInBox[leafI]->containedEdges()];

        forAll(edges, eI)
        {
            const labelList& ef = edgeFaces[edges[eI]];
            if( ef.size() != 2 )
                continue;

            if(
                (
                    (surface_[ef[0]].region() == edgePointRegions[0]) &&
                    (surface_[ef[1]].region() == edgePointRegions[1])
                ) ||
                (
                    (surface_[ef[1]].region() == edgePointRegions[0]) &&
                    (surface_[ef[0]].region() == edgePointRegions[1])
                )
            )
            {
                const edge& edg = surfaceEdges[edges[eI]];

                point nearestOnEdge, nearestOnLine;
                if(
                    help::nearestEdgePointToTheLine
                    (
                        points[edg[0]],
                        points[edg[1]],
                        edgePoints[0],
                        edgePoints[1],
                        nearestOnEdge,
                        nearestOnLine
                    )
                )
                {
                    if( magSqr(nearestOnEdge - nearestOnLine) < distSq )
                    {
                        nearest = nearestOnEdge;
                        distSq = magSqr(nearestOnEdge - nearestOnLine);
                        found = true;
                    }
                }
            }
        }
    }

    return found;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
