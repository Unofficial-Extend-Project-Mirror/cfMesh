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
#include "polyMeshGenModifier.H"
#include "edgeExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceOptimizer.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
#include "DynList.H"
#include "labelPair.H"
#include "labelledScalar.H"
#include "labelledPoint.H"
#include "refLabelledPoint.H"
#include "HashSet.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceClassifyEdges.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGEdgeExtractor

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

bool edgeExtractor::findCornerCandidates()
{
    bool changed(false);

    const triSurf& surf = meshOctree_.surface();
    const pointField& sp = surf.points();

    //- create the surface partitioner and find the corners on the surface mesh
    triSurfacePartitioner surfPartitioner(surf);

    const labelList& corners = surfPartitioner.corners();

    Map<label> surfacePointToCorner;
    forAll(corners, cornerI)
        surfacePointToCorner.insert(corners[cornerI], cornerI);

    List<labelledScalar> nearestPoint
    (
        corners.size(),
        labelledScalar(0, VGREAT)
    );

    //- calculate the search range
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    scalarList searchRange(bPoints.size(), 0.0);
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const point c = bf.centre(points);

        forAll(bf, pI)
        {
            const scalar d = 2.0 * Foam::mag(c - points[bf[pI]]);
            const label bpI = bp[bf[pI]];

            searchRange[bpI] = Foam::max(d, searchRange[bpI]);
        }
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.beAtProcs();
        const DynList<label>& bpNeiProcs = mse.bpNeiProcs();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();

        std::map<label, LongList<labelledScalar> > exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], LongList<labelledScalar>())
            );

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(iter.key(), searchRange[bpI])
                );
            }
        }

        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i].scalarLabel()];

            searchRange[bpI] =
                Foam::max(searchRange[bpI], receivedData[i].value());
        }
    }

    DynList<label> containedTriangles;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    # endif
    forAll(bPoints, bpI)
    {
        const point& p = points[bPoints[bpI]];
        const scalar s = searchRange[bpI];

        boundBox bb(p-vector(s, s, s), p+vector(s, s, s));
        meshOctree_.findTrianglesInBox(bb, containedTriangles);

        //- find the nearest corner on the surface mesh
        forAll(containedTriangles, i)
        {
            const labelledTri& tri = surf[containedTriangles[i]];

            forAll(tri, pI)
            {
                const label spI = tri[pI];
                if( !surfacePointToCorner.found(spI) )
                    continue;

                const label cornerI = surfacePointToCorner[spI];

                const scalar d = Foam::magSqr(sp[spI] - points[bPoints[bpI]]);

                if( nearestPoint[cornerI].value() > d )
                {
                    //- update the nearest point to the found corner
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        nearestPoint[cornerI] = labelledScalar(bpI, d);
                    }
                }
            }
        }
    }

    # ifdef DEBUGEdgeExtractor
    Info << "Nearest points to corners " << nearestPoint << endl;
    # endif

    return changed;
}

bool edgeExtractor::checkCorners()
{
    bool changed(false);

    # ifdef DEBUGEdgeExtractor
    const triSurf* surfPtr = surfaceWithPatches();
    surfPtr->writeSurface("checkCornersBefore.fms");
    deleteDemandDrivenData(surfPtr);
    # endif

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const triSurfacePartitioner& sPartitioner = this->partitioner();
    const List<labelHashSet>& patchPatches = sPartitioner.patchPatches();

    //- allocate a copy of boundary patches
    labelList newBoundaryPatches(facePatch_.size());

    label nCorrected;
    Map<label> otherProcNewPatch;

    label nIteration(0);

    do
    {
        nCorrected = 0;
        newBoundaryPatches = facePatch_;

        //- check whether there exist situations where a boundary face
        //- is surrounded by more faces in different patches than the
        //- faces in the current patch
        if( Pstream::parRun() )
        {
            findOtherFacePatchesParallel
            (
                otherProcNewPatch,
                &newBoundaryPatches
            );
        }

        //- update the information which edges are assigned as feature edges
        meshSurfacePartitioner mPart(mse, newBoundaryPatches);

        std::map<std::pair<label, label>, label> patchNeighbourStatistics;

        //- find corners in the current constelation
        const std::map<label, DynList<label> >& cornerPatches =
            mPart.cornerPatches();

        std::map<label, label> cornerIndex;
        std::map<label, std::pair<point, scalar> > nearestToCorner;

        typedef std::map<label, DynList<label> > mapType;
        typedef std::map<label, std::pair<point, scalar> > distMapType;

        # ifdef DEBUGEdgeExtractor
        Info << "Finding corners " << endl;
        # endif

        //- find nearest corners in the surface mesh
        forAllConstIter(mapType, cornerPatches, it)
        {
            const label bpI = it->first;
            const DynList<label>& neiPatches = it->second;

            const point& p = points[bPoints[bpI]];
            point pMap;
            scalar dSq;
            label nsp;

            if( !meshOctree_.findNearestCorner(pMap, dSq, nsp, p, neiPatches) )
            {
                nsp = -1;
                meshOctree_.findNearestPointToPatches(pMap, dSq, p, neiPatches);
            }

            nearestToCorner[bpI] = std::make_pair(pMap, dSq);
            cornerIndex[bpI] = nsp;
        }

        std::map<label, DynList<label> > bpMappedAtSurfaceCorner;
        for
        (
            std::map<label, label>::const_iterator mIt=cornerIndex.begin();
            mIt!=cornerIndex.end();
            ++mIt
        )
        {
            if( mIt->second < 0 )
            {
                //- the corner does not exist in the surface mesh
                //- this may be a situation where parts of the surface are in
                //- close proximity and are not topologically connected
                Warning << "Should not get in here. Not implementes" << endl;
            }
            else
            {
                bpMappedAtSurfaceCorner[mIt->second].append(mIt->first);
            }
        }

        # ifdef DEBUGEdgeExtractor
        for(mapType::const_iterator mIt=bpMappedAtSurfaceCorner.begin();mIt!=bpMappedAtSurfaceCorner.end();++mIt)
            if( mIt->second.size() > 1 )
                Info << "Surface corner " << mIt->first << " mapped mesh points " << mIt->second << endl;
        # endif

        //- for all edge nodes find their nearest counterparts in the surface
        # ifdef DEBUGEdgeExtractor
        Info << "Finding nearest edges" << endl;
        # endif

        const labelHashSet& featureEdge = mPart.featureEdges();

        std::map<label, std::pair<point, scalar> > nearestToEdgePoint;
        std::map<label, label> edgePointIndex;
        std::map<std::pair<label, label>, labelLongList> edgesSharedByPatches;

        forAllConstIter(labelHashSet, featureEdge, it)
        {
            const label beI = it.key();

            //- get the patches bounded by this feature edge
            DynList<label> patches(2);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = newBoundaryPatches[edgeFaces(beI, 1)];
            }
            else if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = otherProcNewPatch[beI];
            }

            //- store the edge into the right group
            std::pair<label, label> patchPair
            (
                Foam::min(patches[0], patches[1]),
                Foam::max(patches[0], patches[1])
            );
            edgesSharedByPatches[patchPair].append(beI);

            //- check if some points have already been checked
            const edge& e = edges[beI];
            const label bps = bp[e.start()];
            const label bpe = bp[e.end()];

            DynList<label> checkPoints;
            if
            (
                (cornerPatches.find(bps) == cornerPatches.end()) &&
                (nearestToEdgePoint.find(bps) == nearestToEdgePoint.end())
            )
                checkPoints.append(bps);

            if
            (
                (cornerPatches.find(bpe) == cornerPatches.end()) &&
                (nearestToEdgePoint.find(bpe) == nearestToEdgePoint.end())
            )
                checkPoints.append(bpe);

            forAll(checkPoints, i)
            {
                const point& p = points[bPoints[checkPoints[i]]];
                point pMap;
                scalar dSq;
                label nse;

                if( !meshOctree_.findNearestEdgePoint(pMap, dSq, nse, p, patches) )
                {
                    nse = -1;
                    meshOctree_.findNearestPointToPatches(pMap, dSq, p, patches);
                }

                nearestToEdgePoint[checkPoints[i]] = std::make_pair(pMap, dSq);
                edgePointIndex[checkPoints[i]] = nse;
            }
        }

        labelHashSet invalidFeatureEdges;
        for
        (
            std::map<std::pair<label, label>, labelLongList>::const_iterator mIt=edgesSharedByPatches.begin();
            mIt!=edgesSharedByPatches.end();
            ++mIt
        )
        {
            const bool validConnection =
                patchPatches[mIt->first.first].found(mIt->first.second);

            # ifdef DEBUGEdgeExtractor
            Info << "Patches " << mIt->first.first << " and " << mIt->first.second
                 << " share " << mIt->second << " edges " << validConnection << endl;
            # endif

            if( !validConnection )
            {
                //- find surface facets in the vicinity of the edge
                //- and check whether there exist various disconnected surface
                //- parts in the vicinity of this group of edge

                const DynList<label>& invalidEdges = mIt->second;

                forAll(invalidEdges, i)
                    invalidFeatureEdges.insert(invalidEdges[i]);
            }
        }

        # ifdef DEBUGEdgeExtractor
        Info << "Invalid feature edges " << invalidFeatureEdges << endl;
        ::exit(0);
        # endif

        //- check the vicinity of the corner point and check whether
        //- it shall be replaced by some other point in the vicinity
        std::map<label, DynList<labelPair> > facesAtCornerNewPatches;
        forAllConstIter(distMapType, nearestToCorner, it)
        {
            const label bpI = it->first;

            //- find all points connected to the corner point via a face
            labelHashSet nearPoints;
            std::map<label, DynList<label> > facesContainingPoint;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const face& bf = bFaces[bfI];

                forAll(bf, pI)
                {
                    nearPoints.insert(bf[pI]);
                    facesContainingPoint[bf[pI]].append(bfI);
                }
            }

            # ifdef DEBUGEdgeExtractor
            forAllConstIter(labelHashSet, nearPoints, iterN)
            {
                Info << "Near point " << iterN.key() << endl;
                Info << "Faces containing near point are " << facesContainingPoint[iterN.key()] << endl;
            }
            # endif

            //- find the nearest point to the location where the corner
            //- shall be mapped onto the surface mesh
            label bestPoint(-1);
            scalar minDistSq(VGREAT);
            forAllConstIter(labelHashSet, nearPoints, iter)
            {
                const point& p = points[iter.key()];

                const scalar dSq = magSqr(p - nearestToCorner[bpI].first);

                if( dSq < minDistSq )
                {
                    minDistSq = dSq;
                    bestPoint = iter.key();
                }
            }

            if( bestPoint == bPoints[bpI] )
            {
                # ifdef DEBUGEdgeExtractor
                Info << "Current corner is nearest" << endl;
                # endif

                continue;
            }

            # ifdef DEBUGEdgeExtractor
            surfPtr = surfaceWithPatches(bpI);
            surfPtr->writeSurface("corner_"+help::scalarToText(bpI)+"_iter_"+help::scalarToText(nIteration)+".fms");
            deleteDemandDrivenData(surfPtr);
            Info << "Best candidate " << bestPoint << endl;
            # endif

            //- sort faces and edges at the corner in the counter-clockwise order
            DynList<label> pFaces, pEdges;

            DynList<label> front;
            front.append(0);

            while( front.size() )
            {
                const label fI = front.removeLastElement();
                const label bfI = pointFaces(bpI, fI);

                pFaces.append(bfI);

                const face& bf = bFaces[bfI];
                const label pos = bf.which(bPoints[bpI]);

                const label beI = faceEdges(bfI, bf.rcIndex(pos));

                pEdges.append(beI);

                label nei = edgeFaces(beI, 0);
                if( nei == bfI )
                    nei = edgeFaces(beI, 1);

                if( pFaces.contains(nei) || (nei < 0) )
                    continue;

                front.append(pointFaces.containsAtPosition(bpI, nei));
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Boundary point " << bpI
                 << " pFaces " << pFaces
                 << " pEdges " << pEdges << endl;
            forAll(pFaces, i)
            {
                Info << "Face " << pFaces[i] << " has nodes " << bFaces[pFaces[i]] << endl;
                Info << "Face " << pFaces[i] << " is in patch " << facePatch_[pFaces[i]] << endl;
                Info << "Edge " << pEdges[i] << " has nodes " << edges[pEdges[i]] << endl;
            }
            # endif

            //- find the best fitting edge to move towards the corner
            label bestEdge(-1);
            scalar bestAlignment(0.0);

            vector dv = it->second.first - points[bPoints[bpI]];
            dv /= (mag(dv) + VSMALL);

            if( facesContainingPoint[bestPoint].size() == 1 )
            {
                const label bfI = facesContainingPoint[bestPoint][0];

                //- only the edges of this face can be candidates
                forAll(pEdges, i)
                {
                    const label beI = pEdges[i];
                    const edge& e = edges[pEdges[i]];

                    if( !(edgeFaces(beI, 0) == bfI || edgeFaces(beI, 1) == bfI) )
                        continue;

                    vector ev
                    (
                        points[e.otherVertex(bPoints[bpI])] - points[bPoints[bpI]]
                    );
                    ev /= (mag(ev) + VSMALL);

                    const scalar metric = 0.5 * (1.0 + (dv & ev));

                    if( metric > bestAlignment )
                    {
                        bestAlignment = metric;
                        bestEdge = pEdges[i];
                    }
                }
            }
            else
            {
                //- the edge containing the best point is the candidate
                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if( e.otherVertex(bestPoint) != -1 )
                    {
                        //- select the edge which contains best fitting point
                        bestEdge = pEdges[i];
                        break;
                    }
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Best edge " << bestEdge
                 << " has alignment " << bestAlignment << endl;
            # endif

            //- find groups of sorted faces at this corner which are bounded by
            //- existing feature edges and the new candidate
            DynList<label> faceInGroup(pFaces.size(), -1);
            label groupI(0);
            forAll(pFaces, i)
            {
                if( faceInGroup[i] != -1 )
                    continue;

                DynList<label> front;
                front.append(i);
                faceInGroup[i] = groupI;

                while( front.size() != 0 )
                {
                    const label currIndex = front.removeLastElement();

                    const label nbeI = pEdges[currIndex];
                    const label pbeI = pEdges[pFaces.rcIndex(currIndex)];
                    if( (nbeI != bestEdge) && !featureEdge.found(nbeI) )
                    {
                        const label nfI = pFaces.fcIndex(currIndex);

                        if( faceInGroup[nfI] == -1 )
                        {
                            faceInGroup[nfI] = groupI;
                            front.append(nfI);
                        }
                    }
                    else if( (pbeI != bestEdge) && !featureEdge.found(pbeI) )
                    {
                        const label pfI = pFaces.rcIndex(currIndex);

                        if( faceInGroup[pfI] == -1 )
                        {
                            faceInGroup[pfI] = groupI;
                            front.append(pfI);
                        }
                    }
                }

                ++groupI;
            }

            # ifdef DEBUGEdgeExtractor
            Info << "faceInGroup " << faceInGroup << endl;
            # endif

            //- check which group of faces shall change patch in order to make
            //- the best fitting edge a feature edge
            DynList<labelPair> groupPairs, groupPatches;
            labelPair groupsForChanging;
            forAll(pEdges, i)
            {
                const label beI = pEdges[i];

                if( beI == bestEdge )
                {
                    label pos = pFaces.containsAtPosition(edgeFaces(beI, 0));
                    groupsForChanging.first() = faceInGroup[pos];
                    pos = pFaces.containsAtPosition(edgeFaces(beI, 1));
                    groupsForChanging.second() = faceInGroup[pos];
                }
                else if( featureEdge.found(beI) )
                {
                    labelPair lp, lpp;
                    label pos = pFaces.containsAtPosition(edgeFaces(beI, 0));
                    lpp.first() = facePatch_[pFaces[pos]];
                    lp.first() = faceInGroup[pos];
                    pos = pFaces.containsAtPosition(edgeFaces(beI, 1));
                    lpp.second() = facePatch_[pFaces[pos]];
                    lp.second() = faceInGroup[pos];

                    groupPairs.append(lp);
                    groupPatches.append(lpp);
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "group pairs " << groupPairs << endl;
            Info << "Group patches " << groupPatches << endl;
            Info << "Groups for changing " << groupsForChanging << endl;
            # endif

            //- check which groups shall change their patch
            //- desired patches at the best edge are determined by finding
            //- the best alignment between the bestEdge and the feature edges
            scalar maxAlignment(0.0);
            label bestGroupsOfPatches(-1);
            forAll(groupPatches, i)
            {
                const labelPair& pp = groupPatches[i];

                const point& bes = points[edges[bestEdge].start()];
                const point& bee = points[edges[bestEdge].end()];

                DynList<label> patches(2);
                patches[0] = pp.first();
                patches[1] = pp.second();

                point mps, mpe;
                scalar dSqS, dSqE;
                label nse;

                meshOctree_.findNearestEdgePoint(mps, dSqS, nse, bes, patches);
                meshOctree_.findNearestEdgePoint(mpe, dSqE, nse, bee, patches);

                const scalar align = 0.5 * (1.0 + ((bee - bes) & (mpe - mps)));

                if( align > maxAlignment )
                {
                    maxAlignment = align;
                    bestGroupsOfPatches = i;
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Best groups of patches " << bestGroupsOfPatches << endl;
            # endif

            //- assign new patches
            FixedList<label, 2> newPatchForGroup(-1);
            DynList<labelPair>& newPatchForFace = facesAtCornerNewPatches[bpI];
            forAll(groupPatches, i)
            {
                if( i == bestGroupsOfPatches )
                    continue;

                const labelPair& gp = groupPairs[i];
                const labelPair& gpp = groupPatches[i];

                if
                (
                    gp.first() == groupsForChanging.first() ||
                    gp.first() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.second();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if( faceInGroup[j] != gp.first() )
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            # ifdef DEBUGEdgeExtractor
                            Info << "2. beI " << beI << endl;
                            # endif

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if( nei == bfI )
                                nei = edgeFaces(beI, 1);

                            const label posNei = pFaces.containsAtPosition(nei);
                            if( (posNei < 0) || (faceInGroup[posNei] != faceInGroup[j]) )
                            {
                                if( facePatch_[bfI] != facePatch_[nei] )
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];
                                    //- calculate deformation energy
                                    //- of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
                                    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold += 1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                                if( otherPatch != facePatch_[nei] )
                                {
                                    //- calculate deformation energy
                                    //- of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
                                    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew += 1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                            }
                        }
                    }

                    # ifdef DEBUGEdgeExtractor
                    Info << "1. Eold " << Eold << " Enew " << Enew << endl;
                    # endif

                    if( Enew <= Eold )
                    {
                        newPatchForGroup[0] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if( faceInGroup[j] != gp.first() )
                                continue;

                            newPatchForFace.append(labelPair(pFaces[j], otherPatch));
                        }
                    }
                }
                else if
                (
                    gp.second() == groupsForChanging.first() ||
                    gp.second() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.first();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if( faceInGroup[j] != gp.second() )
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            # ifdef DEBUGEdgeExtractor
                            Info << "2. beI " << beI << endl;
                            # endif

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if( nei == bfI )
                                nei = edgeFaces(beI, 1);

                            const label posNei = pFaces.containsAtPosition(nei);
                            if( (posNei < 0) || (faceInGroup[posNei] != faceInGroup[j]) )
                            {
                                if( facePatch_[bfI] != facePatch_[nei] )
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];

                                    //- calculate deformation energy
                                    //- of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
                                    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold += 1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                                if( otherPatch != facePatch_[nei] )
                                {
                                    //- calculate deformation energy
                                    //- of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
                                    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew += 1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                            }
                        }
                    }

                    # ifdef DEBUGEdgeExtractor
                    Info << "2. Eold " << Eold << " Enew " << Enew << endl;
                    # endif

                    if( Enew <= Eold )
                    {
                        newPatchForGroup[1] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if( faceInGroup[j] != gp.second() )
                                continue;

                            newPatchForFace.append(labelPair(pFaces[j], otherPatch));
                        }
                    }
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "New patches for boundary faces " << newPatchForFace << endl;
            # endif
        }

        labelHashSet changedPatch;
        for(std::map<label, DynList<labelPair> >::const_iterator it=facesAtCornerNewPatches.begin();it!=facesAtCornerNewPatches.end();++it)
        {
            const DynList<labelPair>& lp = it->second;
            forAll(lp, i)
            {
                if( changedPatch.found(lp[i].first()) )
                    FatalError << "Face " << lp[i].first() << " is already modified" << abort(FatalError);

                changedPatch.insert(lp[i].first());

                newBoundaryPatches[lp[i].first()] = lp[i].second();
                ++nCorrected;
            }
        }

        //- find the nearest edge points
/*        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            DynList<label> neiPatch(bf.size(), newBoundaryPatches[bfI]);
            DynList<label> ePoints, cPoints;

            forAll(bf, eI)
            {
                const label beI = faceEdges(bfI, eI);
                const label bpI = bp[bf[eI]];

                if( edgePointIndex.find(bpI) != edgePointIndex.end() )
                    ePoints.append(eI);
                if( cornerIndex.find(bpI) != cornerIndex.end() )
                    cPoints.append(eI);

                if( featureEdge.found(beI) )
                {
                    if( edgeFaces.sizeOfRow(beI) == 2 )
                    {
                        label otherFace = edgeFaces(beI, 0);
                        if( otherFace == bfI )
                            otherFace = edgeFaces(beI, 1);

                        neiPatch[eI] = newBoundaryPatches[otherFace];
                    }
                    else if( edgeFaces.sizeOfRow(beI) == 1 )
                    {
                        neiPatch[eI] = otherProcNewPatch[beI];
                    }
                }
            }

            if( cPoints.size() == 0 )
                continue;

            Info << "Face " << bfI << " corner points " << cPoints << endl;
            Info << "Feature edges " << neiPatch << endl;

            DynList<label> nearestPoint(cPoints.size(), -1);
            forAll(cPoints, cpI)
            {
                const label bpI = bp[bf[cPoints[cpI]]];
                //const label nearestSurfaceCorner = cornerIndex[bpI];

                const vector cornerPoint = nearestToCorner[bpI].first;
                //const vector dSq = nearestToCorner[bpI].second;

                //- find face point nearest to the surface corner
                scalar minDist(VGREAT);
                label nearPoint(-1);
                forAll(bf, pI)
                {
                    const point& p = points[bf[pI]];

                    const scalar distSq = magSqr(p - cornerPoint);

                    if( distSq < minDist )
                    {
                        minDist = distSq;
                        nearPoint = pI;
                    }
                }

                nearestPoint[cpI] = nearPoint;
            }

            Info << "Nearest points to corners " << nearestPoint << endl;
            forAll(cPoints, cpI)
            {
                //- check if the corner is closest to itself
                if( nearestPoint[cpI] == cPoints[cpI] )
                    continue;

                //- find the shortest path from the current corner point
                //- to the point nearest to the corner
                scalar distPos(0.0), distNeg(0.0);
                bool otherPosCorner(false), otherNegCorner(false);

                //- calculate squared length of the path in the positive direction
                label index(cPoints[cpI]);
                while( index != nearestPoint[cpI] )
                {
                    const label pos = cPoints.containsAtPosition(index);
                    if( (pos >= 0) && (cPoints[pos] != cPoints[cpI]) )
                        otherPosCorner = true;

                    const edge e = bf.faceEdge(index);
                    distPos += magSqr(points[e[0]] - points[e[1]]);

                    index = bf.fcIndex(index);
                }

                //- calculate squared length of the path in the negative direction
                index = cPoints[cpI];
                while( index != nearestPoint[cpI] )
                {
                    index = bf.rcIndex(index);

                    const label pos = cPoints.containsAtPosition(index);
                    if( (pos >= 0) && (nearestPoint[pos] != nearestPoint[cpI]) )
                        otherNegCorner = true;

                    const edge e = bf.faceEdge(index);
                    distNeg += magSqr(points[e[0]] - points[e[1]]);
                }

                Info << "Dist pos " << distPos << endl;
                Info << "Other corner pos " << otherPosCorner << endl;

                Info << "Dist neg " << distNeg << endl;
                Info << "Other corner neg " << otherNegCorner << endl;

                if( !otherPosCorner && distPos < distNeg )
                {
                    newBoundaryPatches[bfI] = neiPatch[bf.rcIndex(cPoints[cpI])];
                    ++nCorrected;
                }
                else if( !otherNegCorner && (distNeg < distPos) )
                {
                    newBoundaryPatches[bfI] = neiPatch[cPoints[cpI]];
                    ++nCorrected;
                }
            }

        }
        */

        reduce(nCorrected, sumOp<label>());

        if( nCorrected )
        {
            changed = true;
            facePatch_ = newBoundaryPatches;
        }

        if( nIteration++ > 5 )
        {
            Info << "Too many iterations" << endl;
            break;
        }

        # ifdef DEBUGEdgeExtractor
        surfPtr = surfaceWithPatches();
        surfPtr->writeSurface("checkCornersAfterIteration_"+help::scalarToText(nIteration)+".fms");
        deleteDemandDrivenData(surfPtr);
        # endif

    } while( nCorrected != 0 );



    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

