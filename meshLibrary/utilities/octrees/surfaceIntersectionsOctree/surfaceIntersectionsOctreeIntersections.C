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
#include "triSurface.H"
#include "boundBox.H"
#include "DynList.H"
#include "helperFunctions.H"
//#include "DimSpace.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * *Private member functions * * * * * * * * * * * * * * * * * //

bool surfaceIntersectionsOctree::checkPointIntersections
(
    const point& s,
    const pointIndexHit& ph
) const
{
    const labelledTri& ltri = surface_[ph.index()];
    const pointField& points = surface_.points();

    forAll(ltri, pI)
    {
        if( mag(points[ltri[pI]] - ph.hitPoint()) < SMALL )
        {
            //- check visibility of the point s from the ltri
            bool visible(false);
            vector n = ltri.normal(points);
            n /= mag(n);
            if( ((s - ph.hitPoint()) & n) >= 0.0 )
                visible = true;

            //- check visibility of other triangles sharing the intersected
            //- vertex
            const labelList& pp = surface_.pointFaces()[ltri[pI]];

            forAll(pp, nI)
            {
                bool vis(false);
                vector nn = surface_[pp[nI]].normal(points);
                nn /= mag(nn);
                if( ((s - ph.hitPoint()) & nn) >= 0.0 )
                    vis = true;

                if( vis != visible )
                    return false;
            }
        }
    }

    return true;
}

bool surfaceIntersectionsOctree::checkEdgeIntersections
(
    const point& s,
    const pointIndexHit& ph
) const
{
    const pointField& points = surface_.points();
    const labelledTri& ltri = surface_[ph.index()];
    vector n = ltri.normal(points);

    const edgeList edges = ltri.edges();

    forAll(edges, eI)
    {
        vector lv = n ^ edges[eI].vec(points);
        lv /= mag(lv);

        vector pv = ph.hitPoint() - points[ltri[eI]];
        pv /= mag(pv);

        const scalar d = pv & lv;

        if( (d > -SMALL) && (d < SMALL) )
        {
            const label edgeI = surface_.faceEdges()[ph.index()][eI];

            const labelList& ef = surface_.edgeFaces()[edgeI];

            if( ef.size() != 2 )
                FatalErrorIn
                (
                    "bool surfaceIntersectionsOctree::checkEdgeIntersections"
                    "("
                    "const point& s,"
                    "const pointIndexHit& ph"
                    ") const"
                ) << "Surface is not closed!!" << exit(FatalError);

            //- check visibility of the point s from the first triangle
            bool vis1(false);
            vector n1 = surface_[ef[0]].normal(points);
            n1 /= mag(n1);
            if( ((s - ph.hitPoint()) & n1) >= 0.0 )
                vis1 = true;

            //- check visibility of the point s from the second triangle
            bool vis2(false);
            vector n2 = surface_[ef[1]].normal(points);
            n2 /= mag(n2);
            if( ((s - ph.hitPoint()) & n2) >= 0.0 )
                vis2 = true;

            if( vis1 != vis2 )
                return false;
        }
    }
    return true;
}

void surfaceIntersectionsOctree::checkIntersections
(
    const point& p,
    const point& e,
    const SLList<pointIndexHit>& intersections
) const
{
    Info << "Intersection candidates are " << intersections << endl;

    short nIntersections(0);
    point intersection;
    SLList<label> ints;

    forAll(surface_, faceI)
        if(
            help::triLineIntersection
            (
                surface_,
                faceI,
                p,
                e,
                intersection
            )
        )
        {
            ints.append(faceI);
            Info << "Intersection " << nIntersections << " is "
                << intersection << endl;
            nIntersections++;
        }

    Info << "nIntersections " << nIntersections << endl;
    for(SLList<label>::const_iterator pIter = ints.begin();
        pIter != ints.end();
        ++pIter
    )
    {
        bool found(false);
        for(
            SLList<pointIndexHit>::const_iterator iIter =
                intersections.begin();
            iIter != intersections.end();
            ++iIter
        )
            if( pIter() == iIter().index() )
            {
                found = true;
                break;
            }

        if( !found )
            FatalErrorIn
            (
                "bool surfaceIntersectionsOctree::isPointInside(const point& p) const"
            ) << "Intersection " << pIter() << " was not found by the octree!"
                << abort(FatalError);     
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool surfaceIntersectionsOctree::isPointInside(const point& p) const
{
    surfaceIntersectionsOctreeCube* oc = findLeafContainingVertex(p);
    if( !oc )
        return false;

    if( oc->cubeType() & surfaceIntersectionsOctreeCube::INSIDE )
    {
        return true;
    }
    else if( oc->cubeType() & surfaceIntersectionsOctreeCube::OUTSIDE )
    {
        return false;
    }

    const point& max = initialCube_.bb().max();
    const point& min = initialCube_.bb().min();

//     if(
//         (p.x() < min.x()) || (p.x() > max.x()) ||
//         (p.y() < min.y()) || (p.y() > max.y()) ||
//         (p.z() < min.z()) || (p.z() > max.z())
//     )
//         return false;

    const point c = (min+max) / 2.0;
    point e(p);

    if( p.x() > c.x() )
    {
        e.x() = max.x();
    }
    else
    {
        e.x() = min.x();
    }

    SLList<pointIndexHit> intersections;
    initialCube_.intersectionCandidates(intersections, p, e);

    # ifdef DEBUGSearch
    checkIntersections(p, e, intersections);
    # endif

    if( intersections.size() == 0 )
    {
        return false;
    }
    else
    {
        pointIndexHit nearest(false, p, -1);
    
        scalar dist(GREAT);

        for(SLList<pointIndexHit>::const_iterator iIter = intersections.begin();
            iIter != intersections.end();
            ++iIter
        )
            if( iIter().hit() )
                if(
                    (mag(iIter().rawPoint() - p) < dist) &&
                    checkPointIntersections(p, iIter()) &&
                    checkEdgeIntersections(p, iIter())
                )
                {
                    dist = mag(iIter().rawPoint() - p);
                    nearest = iIter();
                }

        if( !nearest.hit() )
            return false;

        const pointField& points = surface_.points();
        const labelledTri& ltri = surface_[nearest.index()];
        vector n = ltri.normal(points);
        n /= mag(n);

        const scalar d = ((p - points[ltri[0]]) & n);

        if( d <= -SMALL )
        {
            return true;
        }
        else if( d < SMALL && d > -SMALL )
        {
            //- vertex is at the boundary
            //- move it outside for consistency reasons
            point& mp = const_cast<point&>(p);
            mp += 10*SMALL*n;
        }

        return false;
    }
}

pointIndexHit surfaceIntersectionsOctree::intersection(const point& s, const point& e) const
{
    SLList<pointIndexHit> intersections;
    initialCube_.intersectionCandidates(intersections, s, e);

    DynList<pointIndexHit> hits;

    for(SLList<pointIndexHit>::const_iterator pIter = intersections.begin();
        pIter != intersections.end();
        ++pIter
    )
    {
        bool store(true);

        forAll(hits, hI)
            if(
                (pIter().index() == hits[hI].index()) ||
                (mag(pIter().rawPoint() - hits[hI].rawPoint()) < SMALL)
            )
            {
                store = false;
                break;
            }

        if( store ) hits.append(pIter());
    }

    if( hits.size() > 1 )
        FatalErrorIn
        (
            "pointIndexHit surfaceIntersectionsOctree::intersection"
            "("
            "const point& s,"
            "const point& e,"
            ") const"
        ) << "There are more than one intersection of the edge "
            << s << " " << e << " with the boundary!!"
            << exit(FatalError);

    pointIndexHit nearest(false, s, -1);
    
    scalar dist(GREAT);

    for(SLList<pointIndexHit>::const_iterator iIter = intersections.begin();
        iIter != intersections.end();
        ++iIter
    )
        if(
            (mag(iIter().rawPoint() - s) < dist) &&
            checkPointIntersections(s, iIter()) &&
            checkEdgeIntersections(s, iIter())                
        )
        {
            dist = mag(iIter().rawPoint() - s);
            nearest = iIter();
        }

    # ifdef DEBUGSearch
    checkIntersections(s, e, intersections);
    # endif

    if(
        !nearest.hit() ||
        nearest.index() < 0 ||
        nearest.index() >= surface_.localFaces().size()
    )
    {
        Info << "intersections " << intersections << endl;
        Info << "nearest intersection " << nearest << endl;
        FatalErrorIn
        (
            "pointIndexHit surfaceIntersectionsOctree::intersection"
            "("
            "const point& s,"
            "const point& e,"
            ") const"
        ) << "Intersection is invalid!!"
            << exit(FatalError);
    }

    return nearest;
}

void surfaceIntersectionsOctree::intersection
(
    const point& s,
    const point& e,
    pointIndexHit& hs,
    pointIndexHit& he
) const
{
    SLList<pointIndexHit> intersections;
    initialCube_.intersectionCandidates(intersections, s, e);

    DynList<pointIndexHit> hits;

    for(SLList<pointIndexHit>::const_iterator pIter = intersections.begin();
        pIter != intersections.end();
        ++pIter
    )
    {
        bool store(true);

        forAll(hits, hI)
            if(
                (pIter().index() == hits[hI].index()) ||
                (mag(pIter().rawPoint() - hits[hI].rawPoint()) < SMALL)
            )
            {
                store = false;
                break;
            }

        if( store ) hits.append(pIter());
    }

    if( hits.size() > 2 )
        FatalErrorIn
        (
            "void surfaceIntersectionsOctree::intersection"
            "("
            "const point& s,"
            "const point& e,"
            "pointIndexHit& hs,"
            "pointIndexHit& he"
            ") const"
        ) << "There are more than two intersections of the edge "
            << s << " " << e << " with the boundary!!" << exit(FatalError);

    hs.setMiss();
    he.setMiss();

    scalar distS(GREAT);
    scalar distE(GREAT);

    for(SLList<pointIndexHit>::const_iterator iIter = intersections.begin();
        iIter != intersections.end();
        ++iIter
    )
        if( 
            checkPointIntersections(s, iIter()) &&
            checkEdgeIntersections(s, iIter())
        )
        {
            if( mag(iIter().rawPoint() - s) < distS )
            {
                distS = mag(iIter().rawPoint() - s);
                hs = iIter();
            }

            if( mag(iIter().rawPoint() - e) < distE )
            {
                distE = mag(iIter().rawPoint() - e);
                he = iIter();
            }
        }

    if( intersections.size() > 0 )
    {
        if(
            !hs.hit() ||
            hs.index() < 0 ||
            hs.index() >= surface_.localFaces().size()
        )
        {
            Info << "intersections " << intersections << endl;
            Info << "nearest intersection to the starting point " << hs << endl;
            FatalErrorIn
            (
                "pointIndexHit surfaceIntersectionsOctree::intersection"
                "("
                "const point& s,"
                "const point& e,"
                ") const"
            ) << "Intersection is invalid!!"
                << exit(FatalError);
        }
        
        if(
            !he.hit() ||
            he.index() < 0 ||
            he.index() >= surface_.localFaces().size()
        )
        {
            Info << "intersections " << intersections << endl;
            Info << "nearest intersection to the end point " << he << endl;
            FatalErrorIn
            (
                "pointIndexHit surfaceIntersectionsOctree::intersection"
                "("
                "const point& s,"
                "const point& e,"
                ") const"
            ) << "Intersection is invalid!!"
                << exit(FatalError);
        }
    }

    # ifdef DEBUGSearch
    checkIntersections(s, e, intersections);
    # endif
}

bool surfaceIntersectionsOctree::normalOk(const label faceI, bool& sure) const
{
    const boundBox& bb = initialCube_.bb();
    const point max = bb.max();
    const point min = bb.min();

    const pointField& points = surface_.points();
    const vector normal = surface_[faceI].normal(points);
    const point s = surface_[faceI].centre(points);

    point e;

    bool set(false);
    {
        //- x-min face
        const vector n(-1, 0, 0);
        const scalar t = (n & (min - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.y() > (min.y()-SMALL)) && (i.y() < (max.y()+SMALL)) &&
            (i.z() > (min.z()-SMALL)) && (i.z() < (max.z()+SMALL))
        )
        {
            e = i;
            set = true;
        }
    }

    if( !set )
    {
        //- x-max face
        const vector n(1, 0, 0);
        const scalar t = (n & (max - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.y() > (min.y()-SMALL)) && (i.y() < (max.y()+SMALL)) &&
            (i.z() > (min.z()-SMALL)) && (i.z() < (max.z()+SMALL))
        )
        {
            e = i;
            set = true;
        }         
    }
    if( !set )
    {
        //- y-min face
        const vector n(0, -1, 0);
        const scalar t = (n & (min - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.x() > (min.x()-SMALL)) && (i.x() < (max.x()+SMALL)) &&
            (i.z() > (min.z()-SMALL)) && (i.z() < (max.z()+SMALL))
        )
        {
            e = i;
            set = true;
        }
    }
    if( !set )
    {
        //- y-max face
        const vector n(0, 1, 0);
        const scalar t = (n & (max - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.x() > (min.x()-SMALL)) && (i.x() < (max.x()+SMALL)) &&
            (i.z() > (min.z()-SMALL)) && (i.z() < (max.z()+SMALL))
        )
        {
            e = i;
            set = true;
        }
    }
    if( !set )
    {
        //- z-min face
        const vector n(0, 0, -1);
        const scalar t = (n & (min - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.x() > (min.x()-SMALL)) && (i.x() < (max.x()+SMALL)) &&
            (i.y() > (min.y()-SMALL)) && (i.y() < (max.y()+SMALL))
        )
        {
            e = i;
            set = true;
        }
    }
    if( !set )
    {
        //- z-min face
        const vector n(0, 0, 1);
        const scalar t = (n & (max - s)) / (n & normal);
        const point i = s + t * normal;
        if(
            (t < SMALL) &&
            (i.x() > (min.x()-SMALL)) && (i.x() < (max.x()+SMALL)) &&
            (i.y() > (min.y()-SMALL)) && (i.y() < (max.y()+SMALL))
        )
        {
            e = i;
            set = true;
        }
    }

    SLList<pointIndexHit> intersections;
    intersections.append(pointIndexHit(true,s, faceI));
    initialCube_.intersectionCandidates(intersections, s, e);

    # ifdef DEBUGSearch
    checkIntersections(s, e, intersections);
    # endif

    List<pointIndexHit> hits(intersections);
    intersections.clear();

    forAll(hits, hitI)
    {
        bool found(false);

        for(label hitJ=hitI+1;hitJ<hits.size();hitJ++)
            if
            (
                (hits[hitI].index() == hits[hitJ].index()) ||
                (mag(hits[hitI].rawPoint() - hits[hitJ].rawPoint()) < SMALL)
            )
                found = true;
            
        if(
            !checkEdgeIntersections(s, hits[hitI]) ||
            !checkPointIntersections(s, hits[hitI])
        )
        {
            sure = false;
            return true;
        }
        else
        {
            sure = true;
        }
            
        
        if( !found ) intersections.append(hits[hitI]);
    }

Info << "Intersections " << intersections << endl;

    if( intersections.size() == 1 )
    {
        return false;
    }
    else
    {
        return !(intersections.size() % 2);
    }
}

bool surfaceIntersectionsOctree::validIntersection
(
    const point& s,
    const point& e,
    SLList<pointIndexHit>& intersections
) const
{
    intersections.clear();
    initialCube_.intersectionCandidates(intersections, s, e);

    //- check for duplicate points
    DynList<pointIndexHit> hits;
    for(SLList<pointIndexHit>::const_iterator pIter = intersections.begin();
        pIter != intersections.end();
        ++pIter
    )
    {
        bool store(true);

        forAll(hits, hI)
            if(
                hits[hI].index() == pIter().index() ||
                mag(hits[hI].rawPoint() - pIter().rawPoint()) < 1e-10
            )
            {
                store = false;
                break;
            }

        if( store ) hits.append(pIter());
    }
    
    intersections.clear();
    forAll(hits, hI)
        intersections.append(hits[hI]);

    //- check hits
    for(SLList<pointIndexHit>::const_iterator pIter = intersections.begin();
        pIter != intersections.end();
        ++pIter
    )
    {
        const pointField& points = surface_.points();
        const labelledTri& ltri = surface_[pIter().index()];
        vector n = ltri.normal(points);

        const edgeList edges = ltri.edges();

        forAll(edges, eI)
        {
            vector lv = n ^ edges[eI].vec(points);
            lv /= mag(lv);
            
            vector pv = pIter().hitPoint() - points[ltri[eI]];
            if( mag(pv) < SMALL )
            {
                //- intersection is located at a vertex of the triangle
                return false;
            }
            pv /= mag(pv);
            
            const scalar d = pv & lv;
            
            if( (d > -SMALL) && (d < SMALL) )
            {
                const labelList& fe = surface_.faceEdges()[pIter().index()];
                const labelList& ef = surface_.edgeFaces()[fe[eI]];
                label neighbourTri(-1);
                if( ef[0] == pIter().index() )
                {
                    neighbourTri = ef[1];
                }
                else
                {
                    neighbourTri = ef[0];
                }
                //- intersection is on the edge of the triangle
                if( surface_[neighbourTri].region() != ltri.region() )
                    return false;
            }
        }
    }

    if( intersections.size() > 2 ) return false;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
