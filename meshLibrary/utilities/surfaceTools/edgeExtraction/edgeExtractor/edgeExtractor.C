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
// Private member functions

void edgeExtractor::calculateValence()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    pointValence_.setSize(mse.boundaryPoints().size());
    pointValence_ = 0;

    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();

    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            ++pointValence_[bp[bf[pI]]];
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& bpNeiProcs = mse.bpNeiProcs();

        std::map<label, LongList<labelPair> > exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], LongList<labelPair>())
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
                    labelPair(iter.key(), pointValence_[bpI])
                );
            }
        }

        LongList<labelPair> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelPair& lp = receivedData[i];

            pointValence_[globalToLocal[lp.first()]] += lp.second();
        }
    }
}

void edgeExtractor::calculateSingleCellEdge()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& faceCells = mse.faceOwners();

    //- find the number of boundary faces for each cell in the mesh
    edgeType_.setSize(edgeFaces.size());
    edgeType_ = NONE;

    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) == 2 )
        {
            const label c0 = faceCells[edgeFaces(eI, 0)];
            const label c1 = faceCells[edgeFaces(eI, 1)];

            if( c0 == c1 )
                edgeType_[eI] |= SINGLECELLEDGE;
        }
    }
}

void edgeExtractor::findFeatureEdgesNearEdge()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const edgeList& edges = mse.edges();

    featureEdgesNearEdge_.setSize(edges.size());
    labelLongList nFeatureEdgesAtEdge(edges.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        labelLongList localData;
        DynList<label> nearEdges;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];
            const vector c = e.centre(points);
            const scalar d = 1.5 * e.mag(points);

            const boundBox bb(c - vector(d, d, d), c + vector(d, d, d));

            //- get the edges near the current edge
            meshOctree_.findEdgesInBox(bb, nearEdges);
            forAllReverse(nearEdges, i)
            {
                const label pos = nearEdges.containsAtPosition(nearEdges[i]);

                if( pos < i )
                    nearEdges.removeElement(i);
            }

            localData.append(edgeI);
            nFeatureEdgesAtEdge[edgeI] = nearEdges.size();
            forAll(nearEdges, i)
                localData.append(nearEdges[i]);
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        featureEdgesNearEdge_.setSizeAndRowSize(nFeatureEdgesAtEdge);

        # pragma omp barrier
        # else
        featureEdgesNearEdge_.setSizeAndRowSize(nFeatureEdgesAtEdge);
        # endif

        //- copy the data to the graph
        label counter(0);
        while( counter < localData.size() )
        {
            const label edgeI = localData[counter++];

            const label size = nFeatureEdgesAtEdge[edgeI];

            for(label i=0;i<size;++i)
                featureEdgesNearEdge_(edgeI, i) = localData[counter++];
        }
    }

    //Info << "Feature edges near edge " << featureEdgesNearEdge_ << endl;
    //::exit(0);
}

void edgeExtractor::markPatchPoints(boolList& patchPoint)
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& bp = mse.bp();

    patchPoint.setSize(bPoints.size());
    patchPoint = true;

    std::map<label, label> otherProcPatch;
    if( Pstream::parRun() )
    {
        const Map<label>& otherProc = mse.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndEdgeAddressing();

        //- create communication matrix
        std::map<label, labelLongList> exchangeData;
        const DynList<label>& neiProcs = mse.beNeiProcs();
        forAll(neiProcs, procI)
            exchangeData.insert
            (
                std::make_pair(neiProcs[procI], labelLongList())
            );

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                labelLongList& dts = exchangeData[otherProc[beI]];
                //- send data as follows:
                //- 1. global edge label
                //- 2. patch of the attached boundary face
                dts.append(it.key());
                dts.append(facePatch_[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocal[receivedData[counter++]];
            const label fPatch = receivedData[counter++];

            otherProcPatch[beI] = fPatch;
        }
    }

    //- set the patchPoint to false for all vertices at feature edges
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(edgeFaces, beI)
    {
        if( edgeFaces.sizeOfRow(beI) == 2 )
        {
            //- an ordinary edge
            if( facePatch_[edgeFaces(beI, 0)] != facePatch_[edgeFaces(beI, 1)] )
            {
                const edge& e = edges[beI];
                patchPoint[bp[e.start()]] = false;
                patchPoint[bp[e.end()]] = false;
            }
        }
        else if( edgeFaces.sizeOfRow(beI) == 1 )
        {
            //- an edge at a parallel interface
            const label otherPatch = otherProcPatch[beI];

            if( facePatch_[edgeFaces(beI, 0)] != otherPatch )
            {
                const edge& e = edges[beI];
                patchPoint[bp[e.start()]] = false;
                patchPoint[bp[e.end()]] = false;
            }
        }
        else
        {
            //- this is a non-manifold edge
            const edge& e = edges[beI];
            patchPoint[bp[e.start()]] = false;
            patchPoint[bp[e.end()]] = false;
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that the information is spread to all processors
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const labelList& globalPointLabel =
            mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();


        std::map<label, labelLongList> sendData;
        forAll(neiProcs, i)
            sendData.insert(std::make_pair(neiProcs[i], labelLongList()));

        forAll(bpAtProcs, bpI)
        {
            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc != Pstream::myProcNo() )
                    sendData[neiProc].append(globalPointLabel[bpI]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(sendData, receivedData);

        forAll(receivedData, i)
                patchPoint[globalToLocal[receivedData[i]]] = false;
    }
}

const meshSurfaceEngine& edgeExtractor::surfaceEngine() const
{
    if( !surfaceEnginePtr_ )
    {
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            if( !surfaceEnginePtr_ )
            {
                surfaceEnginePtr_ = new meshSurfaceEngine(mesh_);
            }
        }
    }

    return *surfaceEnginePtr_;
}

const triSurfacePartitioner& edgeExtractor::partitioner() const
{
    if( !surfPartitionerPtr_ )
    {
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            if( !surfPartitionerPtr_ )
            {
                surfPartitionerPtr_ =
                    new triSurfacePartitioner(meshOctree_.surface());
            }
        }
    }

    return *surfPartitionerPtr_;
}

const triSurfaceClassifyEdges& edgeExtractor::edgeClassifier() const
{
    if( !surfEdgeClassificationPtr_ )
    {
        surfEdgeClassificationPtr_ =
            new triSurfaceClassifyEdges(meshOctree_);
    }

    return *surfEdgeClassificationPtr_;
}

void edgeExtractor::findFaceCandidates
(
    labelLongList& faceCandidates,
    const labelLongList* facePatchPtr,
    const Map<label>* otherFacePatchPtr
) const
{
    faceCandidates.clear();
    if( !facePatchPtr )
        facePatchPtr = &facePatch_;

    const labelLongList& fPatches = *facePatchPtr;

    if( !otherFacePatchPtr )
    {
        Map<label> otherFacePatch;
        findOtherFacePatchesParallel(otherFacePatch, &fPatches);

        otherFacePatchPtr = &otherFacePatch;
    }

    const Map<label>& otherFacePatch = *otherFacePatchPtr;

    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    # ifdef USE_OMP
    # pragma omp parallel if( faceEdges.size() > 1000 )
    # endif
    {
        # ifdef USE_OMP
        labelLongList procCandidates;
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(faceEdges, bfI)
        {
            DynList<label> allNeiPatches;
            forAllRow(faceEdges, bfI, eI)
            {
                const label beI = faceEdges(bfI, eI);

                if( edgeFaces.sizeOfRow(beI) == 2 )
                {
                    label fNei = edgeFaces(beI, 0);
                    if( fNei == bfI )
                        fNei = edgeFaces(faceEdges(bfI, eI), 1);

                    allNeiPatches.appendIfNotIn(fPatches[fNei]);
                }
                else if( edgeFaces.sizeOfRow(beI) == 1 )
                {
                    allNeiPatches.appendIfNotIn(otherFacePatch[beI]);
                }
            }

            if( allNeiPatches.size() > 1 )
            {
                //- this face is probably near some feature edge
                # ifdef USE_OMP
                procCandidates.append(bfI);
                # else
                faceCandidates.append(bfI);
                # endif
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        {
            forAll(procCandidates, i)
                faceCandidates.append(procCandidates[i]);
        }
        # endif
    }
}

void edgeExtractor::findOtherFacePatchesParallel
(
    Map<label>& otherFacePatch,
    const labelLongList* facePatchPtr
) const
{
    otherFacePatch.clear();

    if( !facePatchPtr )
        facePatchPtr = &facePatch_;

    const labelLongList& fPatches = *facePatchPtr;

    if( Pstream::parRun() )
    {
        const meshSurfaceEngine& mse = this->surfaceEngine();
        const VRWGraph& edgeFaces = mse.edgeFaces();
        const Map<label>& otherProc = mse.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndEdgeAddressing();

        //- create communication matrix
        std::map<label, labelLongList> exchangeData;
        const DynList<label>& neiProcs = mse.beNeiProcs();
        forAll(neiProcs, procI)
            exchangeData.insert
            (
                std::make_pair(neiProcs[procI], labelLongList())
            );

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                labelLongList& dts = exchangeData[otherProc[beI]];
                //- send data as follows:
                //- 1. global edge label
                //- 2. patch of the attached boundary face
                dts.append(it.key());
                dts.append(fPatches[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocal[receivedData[counter++]];
            const label fPatch = receivedData[counter++];

            otherFacePatch.insert(beI, fPatch);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

edgeExtractor::edgeExtractor
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    surfaceEnginePtr_(NULL),
    meshOctree_(octree),
    surfPartitionerPtr_(NULL),
    surfEdgeClassificationPtr_(NULL),
    pointValence_(),
    facePatch_(),
    edgeType_(),
    featureEdgesNearEdge_()
{
    calculateValence();

    calculateSingleCellEdge();

    findFeatureEdgesNearEdge();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Destructor

edgeExtractor::~edgeExtractor()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
    deleteDemandDrivenData(surfEdgeClassificationPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void edgeExtractor::moveVerticesTowardsDiscontinuities(const label nIterations)
{
    Info << "Reducing Hausdorff distance:" << flush;

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointFaces = mse.pointFaces();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();

    meshSurfaceEngineModifier modifier(mse);

    vectorField faceCentreDisplacement(bFaces.size());
    List<labelledPoint> pointDisplacements(bPoints.size());

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            //- find displacements of face centres
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 40)
            # endif
            forAll(bFaces, bfI)
            {
                const vector centre = bFaces[bfI].centre(points);

                point newP;
                scalar distSq;
                label patchI;
                meshOctree_.findNearestSurfacePoint
                (
                    newP,
                    distSq,
                    patchI,
                    centre
                );

                faceCentreDisplacement[bfI] = newP - centre;
            }

            //- initialise displacements
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 40)
            # endif
            forAll(pointDisplacements, bpI)
                pointDisplacements[bpI] = labelledPoint(0, vector::zero);

            # ifdef USE_OMP
            # pragma omp barrier
            # endif

            //- calculate displacements of boundary points as the average
            //- of face centre displacements
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 40)
            # endif
            forAll(pointFaces, bpI)
            {
                forAllRow(pointFaces, bpI, pfI)
                {
                    pointDisplacements[bpI].coordinates() +=
                        faceCentreDisplacement[pointFaces(bpI, pfI)];
                    ++pointDisplacements[bpI].pointLabel();
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = mse.bpNeiProcs();
            const VRWGraph& bpAtProcs = mse.bpAtProcs();

            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(neiProcs, i)
                exchangeData[i] = LongList<refLabelledPoint>();

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
                        refLabelledPoint(iter.key(), pointDisplacements[bpI])
                    );
                }
            }

            LongList<refLabelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const label globalLabel = receivedData[i].objectLabel();
                const labelledPoint& lp = receivedData[i].lPoint();

                const label bpI = globalToLocal[globalLabel];

                pointDisplacements[bpI].coordinates() += lp.coordinates();
                pointDisplacements[bpI].pointLabel() += lp.pointLabel();
            }
        }

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40)
        # endif
        forAll(pointDisplacements, bpI)
        {
            const labelledPoint& lp = pointDisplacements[bpI];
            const point mp =
                points[bPoints[bpI]] + lp.coordinates() / lp.pointLabel();

            //Info << "Original point " << bPoints[bpI] << " has coordinates "
            //        << points[bPoints[bpI]] << endl;
            //Info << "Displacement vector " << lp.coordinates() / lp.pointLabel() << endl;
            //Info << "Moved point " << mp << endl;

            point newPoint;
            label patchI;
            scalar distSq;

            meshOctree_.findNearestSurfacePoint(newPoint, distSq, patchI, mp);

            //Info << "Mapped point " << newPoint << nl << endl;

            modifier.moveBoundaryVertexNoUpdate(bpI, newPoint);
        }

        //- update geometry
        modifier.updateGeometry();
        modifier.syncVerticesAtParallelBoundaries();

        Info << '.' << flush;
    }

    Info << endl;
}

void edgeExtractor::distributeBoundaryFaces()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    //const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const pointFieldPMG& points = mse.points();

    //- set the size of the facePatch list
    facePatch_.setSize(bFaces.size());

    //- check if the mesh already has patches
    if( mesh_.boundaries().size() > 1 )
        Warning << "Mesh patches are already assigned!" << endl;

    //- set size of patchNames, newBoundaryFaces_ and newBoundaryOwners_
    const triSurf& surface = meshOctree_.surface();
    const label nPatches = surface.patches().size();

    //- find the region for face by finding the patch nearest
    //- to the face centre
    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 100 ) schedule(dynamic, 40)
    # endif
    forAll(bFaces, bfI)
    {
        const point c = bFaces[bfI].centre(points);

        label fPatch;
        point p;
        scalar distSq;

        meshOctree_.findNearestSurfacePoint(p, distSq, fPatch, c);

        if( (fPatch > -1) && (fPatch < nPatches) )
        {
            facePatch_[bfI] = fPatch;
        }
        else
        {
            FatalErrorIn
            (
                "void meshSurfaceEdgeExtractorNonTopo::"
                "distributeBoundaryFaces()"
            ) << "Cannot distribute a face " << bFaces[bfI] << " into any "
                << "surface patch!. Exiting.." << exit(FatalError);
        }
    }
}

void edgeExtractor::findEdgeCandidates()
{
    const triSurf& surface = meshOctree_.surface();
    const vectorField& sp = surface.points();
    const VRWGraph& facetEdges = surface.facetEdges();
    const VRWGraph& edgeFacets = surface.edgeFacets();

    const triSurfacePartitioner& partitioner = this->partitioner();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const VRWGraph& faceEdges = mse.faceEdges();

    Map<label> otherFacePatch;
    findOtherFacePatchesParallel(otherFacePatch, &facePatch_);
    labelLongList faceCandidates;
    findFaceCandidates(faceCandidates, &facePatch_, &otherFacePatch);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) \
    if( faceCandidates.size() > 100 )
    # endif
    forAll(faceCandidates, fcI)
    {
        const label bfI = faceCandidates[fcI];

        forAllRow(faceEdges, bfI, i)
        {
            const label eI = faceEdges(bfI, i);
            edgeType_[eI] |= CANDIDATE;
        }
    }

    //- find distances of all vertices supporting CANDIDATE edges
    //- from feature edges separating various patches
    const VRWGraph& pEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();

    List<List<labelledScalar> > featureEdgesNearPoint(bPoints.size());

    DynList<label> containedTriangles(100);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    # endif
    forAll(pEdges, bpI)
    {
        // TODO rewrite for execution on distributed machines
        bool check(false);
        forAllRow(pEdges, bpI, peI)
        {
            const label eI = pEdges(bpI, peI);

            if( edgeType_[eI] & CANDIDATE )
            {
                check = true;
                break;
            }
        }

        if( check )
        {
            //- check the squared distance from the nearest feature edge
            scalar rSq(0.0);
            forAllRow(pEdges, bpI, peI)
            {
                const label eI = pEdges(bpI, peI);
                const edge& e = edges[eI];
                const scalar dSq = magSqr(points[e.start()] - points[e.end()]);

                rSq = Foam::max(rSq, dSq);
            }

            rSq *= 2.0;
            const scalar r = Foam::sqrt(rSq);

            //- create a boundaing box used for searching neighbour edges
            const point& p = points[bPoints[bpI]];
            boundBox bb(p - point(r, r, r), p + point(r, r, r));

            //- find the surface triangles in the vicinity of the point
            //- check for potential feature edges
            containedTriangles.clear();
            meshOctree_.findTrianglesInBox(bb, containedTriangles);

            DynList<label> featureEdgeCandidates;

            forAll(containedTriangles, ctI)
            {
                const label tI = containedTriangles[ctI];

                forAllRow(facetEdges, tI, feI)
                {
                    const label seI = facetEdges(tI, feI);

                    if( edgeFacets.sizeOfRow(seI) == 2 )
                    {
                        const label p0 = surface[edgeFacets(seI, 0)].region();
                        const label p1 = surface[edgeFacets(seI, 1)].region();

                        if( p0 != p1 )
                        {
                            featureEdgeCandidates.appendIfNotIn(seI);
                        }
                    }
                    else
                    {
                        featureEdgeCandidates.appendIfNotIn(seI);
                    }
                }
            }

            //- check the distance of the vertex from the candidates
            List<labelledScalar>& featureEdgeDistances =
                featureEdgesNearPoint[bpI];
            featureEdgeDistances.setSize(featureEdgeCandidates.size());
            forAll(featureEdgeCandidates, i)
            {
                const label seI = featureEdgeCandidates[i];

                const point s = sp[edges[seI].start()];
                const point e = sp[edges[seI].end()];
                const point np = help::nearestPointOnTheEdgeExact(s, e, p);

                featureEdgeDistances[i] = labelledScalar(seI, magSqr(np - p));
            }

            //- find nearest edges
            sort(featureEdgeDistances);
        }
    }

    //- start post-processing gathered data
    const labelList& edgePartition = partitioner.edgePartitions();

    List<List<labelledScalar> > edgePartitionAndWeights(edges.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) \
    if( edges.size() > 1000 )
    # endif
    forAll(edgeType_, edgeI)
    {
        if( edgeType_[edgeI] & CANDIDATE )
        {
            const edge& e = edges[edgeI];

            const List<labelledScalar>& sc =
                featureEdgesNearPoint[bp[e.start()]];
            const List<labelledScalar>& ec =
                featureEdgesNearPoint[bp[e.end()]];

            //- find the feature-edge partition for which the sum of
            //- node weights is minimal.
            DynList<labelledScalar> weights;
            forAll(sc, i)
            {
                const label sPart = edgePartition[sc[i].scalarLabel()];

                forAll(ec, j)
                {
                    const label ePart = edgePartition[ec[j].scalarLabel()];

                    if( (sPart >= 0) && (sPart == ePart) )
                    {
                        weights.append
                        (
                            labelledScalar
                            (
                                sPart,
                                sc[i].value() + ec[j].value()
                            )
                        );
                    }
                }
            }

            //- store the data
            edgePartitionAndWeights[edgeI].setSize(weights.size());
            forAll(edgePartitionAndWeights[edgeI], epI)
                edgePartitionAndWeights[edgeI][epI] = weights[epI];

            //- sort the data according to the weights
            stableSort(edgePartitionAndWeights[edgeI]);
        }
    }

    Info << "Edge partitions and weights " << edgePartitionAndWeights << endl;
}

bool edgeExtractor::checkConcaveEdgeCells()
{
    bool changed(false);

    const triSurf& surf = meshOctree_.surface();
    const VRWGraph& edgeFacets = surf.edgeFacets();

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const label bndStartFace = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const labelList& faceCells = mse.faceOwners();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    //- analyse the surface mesh and find out which edges are concave or convex
    const triSurfaceClassifyEdges& edgeClassifier = this->edgeClassifier();
    const List<direction>& edgeType = edgeClassifier.edgeTypes();

    //- create a copy of facePatch array for local modifications
    labelLongList newBoundaryPatches(facePatch_);

    //- start checking the surface of the mesh
    label nChanged;

    boolList patchPoint(mse.boundaryPoints().size(), false);

    do
    {
        nChanged = 0;

        //- check which surface points are surrounded by boundary faces
        //- in the same surface patch
        markPatchPoints(patchPoint);

        //- check whether exist edges of a single cell which shall be projected
        //- onto a concave edge
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40) reduction(+ : nChanged)
        # endif
        forAll(edgeType_, beI)
        {
            if( !(edgeType_[beI] & SINGLECELLEDGE) )
                continue;

            //- check if all faces are assigned to the same patch
            const label firstPatch = newBoundaryPatches[edgeFaces(beI, 0)];
            const label secondPatch = newBoundaryPatches[edgeFaces(beI, 1)];

            if( firstPatch == secondPatch )
                continue;

            const cell& c = cells[faceCells[edgeFaces(beI, 0)]];

            //- find edges within the bounding box determined by the cell
            point pMin(VGREAT, VGREAT, VGREAT);
            point pMax(-VGREAT, -VGREAT, -VGREAT);
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                {
                    pMin = Foam::min(pMin, points[f[pI]]);
                    pMax = Foam::max(pMax, points[f[pI]]);
                }
            }

            const point cc = 0.5 * (pMin + pMax);
            const point diff = pMax - pMin;
            boundBox bb(cc-diff, cc+diff);
            DynList<label> containedEdges;
            meshOctree_.findEdgesInBox(bb, containedEdges);

            //- check if there exists concave edges boundaing patches
            //- assigned to boundary faces of the current cell
            forAll(containedEdges, ceI)
            {
                const label eI = containedEdges[ceI];

                if( edgeFacets.sizeOfRow(eI) != 2 )
                    continue;
                if( !(edgeType[eI] & triSurfaceClassifyEdges::FEATUREEDGE) )
                    continue;

                if( edgeType[eI] & triSurfaceClassifyEdges::CONCAVEEDGE )
                {
                    const label patch0 = surf[edgeFacets(eI, 0)].region();
                    const label patch1 = surf[edgeFacets(eI, 1)].region();

                    if
                    (
                        ((firstPatch == patch0) && (secondPatch == patch1)) ||
                        ((firstPatch == patch1) && (secondPatch == patch0))
                    )
                    {
                        DynList<DynList<label>, 2> facesInPatch;
                        facesInPatch.setSize(2);

                        DynList<label, 2> nFacesInPatch;
                        nFacesInPatch.setSize(2);
                        nFacesInPatch = 0;

                        DynList<bool, 2> hasPatchPoints;
                        hasPatchPoints.setSize(2);
                        hasPatchPoints = false;

                        forAll(c, fI)
                        {
                            if( c[fI] < bndStartFace )
                                continue;

                            const label bfI = c[fI] - bndStartFace;
                            const face& bf = bFaces[bfI];

                            if( newBoundaryPatches[bfI] == patch1 )
                            {
                                facesInPatch[1].append(bfI);
                                ++nFacesInPatch[1];

                                forAll(bf, pI)
                                {
                                    if( patchPoint[bp[bf[pI]]] )
                                        hasPatchPoints[1] = true;
                                }
                            }
                            else if( newBoundaryPatches[bfI] == patch0 )
                            {
                                facesInPatch[0].append(bfI);
                                ++nFacesInPatch[0];

                                forAll(bf, pI)
                                {
                                    if( patchPoint[bp[bf[pI]]] )
                                        hasPatchPoints[0] = true;
                                }
                            }
                        }

                        if( nFacesInPatch[1] > nFacesInPatch[0] )
                        {
                            //- there exist more faces in patch 1
                            //- assign all boundary faces to the same patch
                            forAll(facesInPatch[0], i)
                                newBoundaryPatches[facesInPatch[0][i]] = patch1;
                            ++nChanged;
                            break;
                        }
                        else if( nFacesInPatch[0] > nFacesInPatch[1] )
                        {
                            //- there exist more faces in patch 0
                            //- assign all boundary faces to the same patch
                            forAll(facesInPatch[1], i)
                                newBoundaryPatches[facesInPatch[1][i]] = patch0;
                            ++nChanged;
                            break;
                        }
                        else
                        {
                            if( hasPatchPoints[0] && !hasPatchPoints[1] )
                            {
                                //- transfer all faces to patch 1
                                forAll(facesInPatch[0], i)
                                    newBoundaryPatches[facesInPatch[0][i]] =
                                        patch1;
                                ++nChanged;
                                break;
                            }
                            else if( !hasPatchPoints[0] && hasPatchPoints[1] )
                            {
                                //- transfer all faces to patch 0
                                forAll(facesInPatch[1], i)
                                    newBoundaryPatches[facesInPatch[1][i]] =
                                        patch0;
                                ++nChanged;
                                break;
                            }
                            else
                            {
                                //- just transfer all faces to the same patch
                                forAll(facesInPatch[1], i)
                                    newBoundaryPatches[facesInPatch[1][i]] =
                                        patch0;
                                ++nChanged;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if( Pstream::parRun() )
            reduce(nChanged, sumOp<label>());

        if( nChanged )
            changed = true;

    } while( nChanged != 0 );

    //- transfer the information back to facePatch
    facePatch_.transfer(newBoundaryPatches);

    return changed;
}

bool edgeExtractor::checkFacePatches()
{
    bool changed(false);

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    //- allocate a copy of boundary patches
    labelLongList newBoundaryPatches(facePatch_);

    label nCorrected;
    Map<label> otherProcNewPatch;

    do
    {
        nCorrected = 0;

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

        //- find the faces which have neighbouring faces in other patches
        labelLongList candidates;
        findFaceCandidates(candidates, &newBoundaryPatches, &otherProcNewPatch);

        //- go through the list of faces and check if they shall remain
        //- in the current patch
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40) \
        if( candidates.size() > 1000 ) reduction(+ : nCorrected)
        # endif
        forAll(candidates, i)
        {
            const label bfI = candidates[i];

            DynList<label> allNeiPatches;
            DynList<label> neiPatches;
            neiPatches.setSize(faceEdges.sizeOfRow(bfI));

            forAllRow(faceEdges, bfI, eI)
            {
                const label beI = faceEdges(bfI, eI);

                if( edgeFaces.sizeOfRow(beI) == 2 )
                {
                    label fNei = edgeFaces(beI, 0);
                    if( fNei == bfI )
                        fNei = edgeFaces(faceEdges(bfI, eI), 1);

                    allNeiPatches.appendIfNotIn(newBoundaryPatches[fNei]);
                    neiPatches[eI] = newBoundaryPatches[fNei];
                }
                else if( edgeFaces.sizeOfRow(beI) == 1 )
                {
                    allNeiPatches.appendIfNotIn(otherProcNewPatch[beI]);
                    neiPatches[eI] = otherProcNewPatch[beI];
                }
            }

            //- check if some faces have to be distributed to another patch
            //- in order to reduce the number of feature edges
            if(
                (
                    (allNeiPatches.size() == 1) &&
                    (allNeiPatches[0] == newBoundaryPatches[bfI])
                )
            )
                continue;

            Map<label> nNeiInPatch(allNeiPatches.size());
            forAll(allNeiPatches, i)
                nNeiInPatch.insert(allNeiPatches[i], 0);
            forAll(neiPatches, eI)
                ++nNeiInPatch[neiPatches[eI]];

            label newPatch(-1), nNeiEdges(0);
            forAllConstIter(Map<label>, nNeiInPatch, it)
            {
                if( it() > nNeiEdges )
                {
                    newPatch = it.key();
                    nNeiEdges = it();
                }
                else if
                (
                    (it() == nNeiEdges) && (it.key() == newBoundaryPatches[bfI])
                )
                {
                    newPatch = it.key();
                }
            }

            if( (newPatch < 0) || (newPatch == newBoundaryPatches[bfI]) )
                continue;

            boolList sharedEdge(bFaces[bfI].size(), false);
            forAll(neiPatches, eI)
                if( neiPatches[eI] == newPatch )
                    sharedEdge[eI] = true;

            if( help::areElementsInChain(sharedEdge) )
            {
                //- change the patch to the newPatch
                ++nCorrected;
                newBoundaryPatches[bfI] = newPatch;
            }
        }

        reduce(nCorrected, sumOp<label>());

        if( nCorrected )
            changed = true;

    } while( nCorrected != 0 );

    //- transfer the new patches back
    facePatch_.transfer(newBoundaryPatches);

    return changed;
}

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

    DynList<label> containedTriangles(200);
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

void edgeExtractor::projectDeterminedFeatureVertices()
{
    VRWGraph pointPatches;
    pointPatches.setSizeAndRowSize(pointValence_);
    forAll(pointPatches, pI)
        pointPatches.setRowSize(pI, 0);

    const meshSurfaceEngine& mse = surfaceEngine();
    const pointFieldPMG& points = mse.mesh().points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    //- calculate patches for each point
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            pointPatches.appendIfNotIn(bp[bf[pI]], facePatch_[bfI]);
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& bpNeiProcs = mse.bpNeiProcs();

        std::map<label, LongList<labelPair> > exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], LongList<labelPair>())
            );

        //- collet the data distributed to others
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                LongList<labelPair>& data = exchangeData[neiProc];

                forAllRow(pointPatches, bpI, ppI)
                    data.append(labelPair(it.key(), pointPatches(bpI, ppI)));
            }
        }

        //- exchange information
        LongList<labelPair> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- unify the data
        forAll(receivedData, i)
        {
            const labelPair& lp = receivedData[i];

            pointPatches.appendIfNotIn(globalToLocal[lp.first()], lp.second());
        }
    }

    meshSurfaceEngineModifier surfMod(mse);

    for(label iterI=0;iterI<2;++iterI)
    {
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 10)
        # endif
        forAll(pointPatches, bpI)
        {
            if( pointPatches.sizeOfRow(bpI) < 2 )
                continue;

            const point& p = points[bPoints[bpI]];

            point newP(vector::zero);
            label counter(0);

            forAllRow(pointPatches, bpI, i)
            {
                point pp;
                scalar dSq;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    pp,
                    dSq,
                    pointPatches(bpI, i),
                    p
                );

                newP += pp;
                ++counter;
            }

            newP /= counter;
            surfMod.moveBoundaryVertexNoUpdate(bpI, newP);
        }

        surfMod.syncVerticesAtParallelBoundaries();
        surfMod.updateGeometry();
    }
}

bool edgeExtractor::untangleSurface()
{
    bool changed(false);

    meshSurfaceEngine& mse =
        const_cast<meshSurfaceEngine&>(this->surfaceEngine());
    meshSurfaceOptimizer optimizer(mse, meshOctree_);
    changed = optimizer.preOptimizeSurface();

    return changed;
}

void edgeExtractor::extractEdges()
{
//    bool changed;
//    label nIter(0);

    distributeBoundaryFaces();

    # ifdef DEBUGEdgeExtractor
    const triSurf* sPtr = surfaceWithPatches();
    sPtr->writeSurface("initialDistributionOfPatches.stl");
    deleteDemandDrivenData(sPtr);
    # endif

//    do
//    {
//        changed = false;

        Info << "Checking cells near concave edges" << endl;
        if( checkConcaveEdgeCells() )
        {
            # ifdef DEBUGEdgeExtractor
            Info << "Changes due to concave edge cells" << endl;
            fileName sName("changedConcaveCells"+help::scalarToText(nIter)+".stl");
            sPtr = surfaceWithPatches();
            sPtr->writeSurface(sName);
            deleteDemandDrivenData(sPtr);
            # endif

            //changed = true;
        }

        Info << "Checking patch in the neighbourhood of each face" << endl;
        if( checkFacePatches() )
        {
            # ifdef DEBUGEdgeExtractor
            Info << "Changes due to face patches" << endl;
            fileName sName("checkFacePatches"+help::scalarToText(nIter)+".stl");
            sPtr = surfaceWithPatches();
            sPtr->writeSurface(sName);
            deleteDemandDrivenData(sPtr);
            # endif

            //changed = true;
        }

        //findEdgeCandidates();
/*
        if( findCornerCandidates() )
        {
            # ifdef DEBUGEdgeExtractor
            Info << "Changes due to corner candidates" << endl;
            fileName sName("findCornerCandidates"+help::scalarToText(nIter)+".stl");
            sPtr = surfaceWithPatches();
            sPtr->writeSurface(sName);
            deleteDemandDrivenData(sPtr);
            # endif

            changed = true;
        }

        projectDeterminedFeatureVertices();

        if( untangleSurface() )
        {
            # ifdef DEBUGEdgeExtractor
            Info << "Changes due to untangling" << endl;
            fileName sName("untangleSurface"+help::scalarToText(nIter)+".stl");
            sPtr = surfaceWithPatches();
            sPtr->writeSurface(sName);
            deleteDemandDrivenData(sPtr);
            # endif

            changed = true;
        }
        else
        {
            break;
        }
        */

//    } while( changed && ++nIter < 3 );

    # ifdef DEBUGEdgeExtractor
    const triSurf* sPtr = surfaceWithPatches();
    sPtr->writeSurface("finalDistributionOfPatches.stl");
    deleteDemandDrivenData(sPtr);
    # endif
}

const triSurf* edgeExtractor::surfaceWithPatches() const
{
    //- allocate the memory for the surface mesh
    triSurf* surfPtr = new triSurf();

    //- surface of the volume mesh
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const pointFieldPMG& points = mesh_.points();

    //- modifier of the new surface mesh
    triSurfModifier surfModifier(*surfPtr);
    surfModifier.patchesAccess() = meshOctree_.surface().patches();
    pointField& sPts = surfModifier.pointsAccess();
    sPts.setSize(mse.boundaryPoints().size());

    //- copy points
    forAll(bp, pointI)
    {
        if( bp[pointI] < 0 )
            continue;

        sPts[bp[pointI]] = points[pointI];
    }

    //- create the triangulation of the volume mesh surface
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        labelledTri tri;
        tri.region() = facePatch_[bfI];
        tri[0] = bp[bf[0]];

        for(label i=bf.size()-2;i>0;--i)
        {
            tri[1] = bp[bf[i]];
            tri[2] = bp[bf[i+1]];

            surfPtr->appendTriangle(tri);
        }
    }

    return surfPtr;
}

void edgeExtractor::updateMeshPatches()
{
    const triSurf& surface = meshOctree_.surface();
    const label nPatches = surface.patches().size();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwner = mse.faceOwners();

    wordList patchNames(nPatches);
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners(bFaces.size());

    //- set patchNames
    forAll(surface.patches(), patchI)
        patchNames[patchI] = surface.patches()[patchI].name();

    //- append boundary faces
    forAll(bFaces, bfI)
    {
        newBoundaryFaces.appendList(bFaces[bfI]);
        newBoundaryOwners[bfI] = faceOwner[bfI];
    }

    //- replace the boundary with the new patches
    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        facePatch_
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
