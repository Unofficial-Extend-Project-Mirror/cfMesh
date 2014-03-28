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
    const List<DynList<label> >& cornerPatches =
        surfPartitioner.cornerPatches();

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

    Info << "Nearest points to corners " << nearestPoint << endl;

    return changed;
}

bool edgeExtractor::checkCorners()
{
    bool changed(false);

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const triSurfacePartitioner& sPart = this->partitioner();
    const labelList& surfEdgeIndex = sPart.edgePartitions();

    //- allocate a copy of boundary patches
    labelList newBoundaryPatches(facePatch_.size());

    label nCorrected;
    Map<label> otherProcNewPatch;

    label iter(0);

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
                &facePatch_
            );
        }

        //- update the information which edges are assigned as feature edges
        meshSurfacePartitioner mPart(mse, facePatch_);

        //- find corners in the current constelation and their nearest
        //- counterparts in the suface mesh
        const std::map<label, DynList<label> >& cornerPatches =
            mPart.cornerPatches();

        std::map<label, label> cornerIndex;
        std::map<label, std::pair<point, scalar> > nearestToCorner;

        typedef std::map<label, DynList<label> > mapType;
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

        //- for all edge nodes find their nearest counterparts in the surface
        const labelHashSet& featureEdge = mPart.featureEdges();

        std::map<label, std::pair<point, scalar> > nearestToEdgePoint;
        std::map<label, label> edgePointIndex;

        forAllConstIter(labelHashSet, featureEdge, it)
        {
            const label beI = it.key();

            //- get the patches bounded by this feature edge
            DynList<label> patches(2);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                patches[0] = facePatch_[edgeFaces(beI, 0)];
                patches[1] = facePatch_[edgeFaces(beI, 1)];
            }
            else if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                patches[0] = facePatch_[edgeFaces(beI, 0)];
                patches[1] = otherProcNewPatch[beI];
            }

            //- check if some points have alreay been checked
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


        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            DynList<label> neiPatch(bf.size(), facePatch_[bfI]);

            forAll(bf, eI)
            {
                const label beI = faceEdges(bfI, eI);

                if( featureEdge.found(beI) )
                {
                    if( edgeFaces.sizeOfRow(beI) == 2 )
                    {
                        label otherFace = edgeFaces(beI, 0);
                        if( otherFace == bfI )
                            otherFace = edgeFaces(beI, 1);

                        neiPatch[eI] = facePatch_[otherFace];
                    }
                    else if( edgeFaces.sizeOfRow(beI) == 1 )
                    {
                        neiPatch[eI] = otherProcNewPatch[beI];
                    }
                }
            }
        }

        reduce(nCorrected, sumOp<label>());

        //::exit(0);

        if( nCorrected )
        {
            changed = true;
            facePatch_ = newBoundaryPatches;
        }

        if( ++iter > 0 )
            break;

    } while( nCorrected != 0 );

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

