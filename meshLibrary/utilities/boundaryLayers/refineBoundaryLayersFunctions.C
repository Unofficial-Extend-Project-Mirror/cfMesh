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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "VRWGraphList.H"

#include "labelledPair.H"
#include "labelledScalar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class meshBndLayerNeighbourOperator
{
    const polyMeshGen& mesh_;
    const meshSurfaceEngine& mse_;

public:

    meshBndLayerNeighbourOperator
    (
        const polyMeshGen& mesh,
        const meshSurfaceEngine& mse
    )
    :
        mesh_(mesh),
        mse_(mse)
    {}

    label size() const
    {
        return mesh_.cells().size();
    }

    void operator()(const label cellI, DynList<label>& neighbourCells) const
    {
        neighbourCells.clear();

        const edgeList& edges = mse_.edges();
        const VRWGraph& bpEdges = mse_.boundaryPointEdges();
        const labelList& bp = mse_.bp();

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const faceListPMG& faces = mesh_.faces();
        const cell& c = mesh_.cells()[cellI];

        forAll(c, fI)
        {
            //- find the neighbour cell over this face
            label nei = owner[c[fI]];

            if( nei == cellI )
                nei = neighbour[c[fI]];

            if( nei < 0 )
                continue;

            //- face must be a quad
            const face& f = faces[c[fI]];

            if( f.size() != 4 )
                continue;

            forAll(f, eI)
            {
                const edge e = f.faceEdge(eI);

                const label bpI = bp[e.start()];
                const label bpJ = bp[e.end()];

                if( (bpI < 0) || (bpJ < 0) )
                    continue;

                forAllRow(bpEdges, bpI, peI)
                {
                    const label beI = bpEdges(bpI, peI);

                    if( edges[beI] == e )
                    {
                        neighbourCells.append(nei);
                        break;
                    }
                }
            }
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >& neiGroups,
        const labelListType& elementInGroup,
        const DynList<label>& localGroupLabel
    ) const
    {
        const PtrList<writeProcessorPatch>& procBoundaries =
            mesh_.procBoundaries();
        const labelList& owner = mesh_.owner();

        //- send the data to other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelList groupOwner(procBoundaries[patchI].patchSize());
            for(label faceI=0;faceI<size;++faceI)
            {
                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                {
                    groupOwner[faceI] = -1;
                    continue;
                }

                groupOwner[faceI] = localGroupLabel[groupI];
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                groupOwner.byteSize()
            );

            toOtherProc << groupOwner;
        }

        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            forAll(receivedData, faceI)
            {
                if( receivedData[faceI] < 0 )
                    continue;

                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                    continue;

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(receivedData[faceI]);
            }
        }
    }
};

class meshBndLayerSelectorOperator
{
    const polyMeshGen& mesh_;

public:

    meshBndLayerSelectorOperator(const polyMeshGen& mesh)
    :
        mesh_(mesh)
    {}

    bool operator()(const label cellI) const
    {
        const faceListPMG& faces = mesh_.faces();

        const cell& c = mesh_.cells()[cellI];
        const PtrList<writePatch>& boundaries = mesh_.boundaries();
        const label start = boundaries[0].patchStart();

        label nBndFaces(0), baseFace(-1), otherBase(-1), nQuads(0);
        forAll(c, fI)
        {
            if( faces[c[fI]].size() == 4 )
                ++nQuads;

            if( (c[fI] - start) >= 0 )
            {
                baseFace = fI;
                ++nBndFaces;
            }
            else if( faces[c[fI]].size() != 4 )
            {
                otherBase = fI;
            }
        }

        if( nBndFaces != 1 )
            return false;

        bool isPrism(false);

        if(
            (nQuads == c.size()) ||
            (
                (c.size() - nQuads) == 2 &&
                (baseFace != -1) && (otherBase != -1) &&
                (faces[c[baseFace]].size() == faces[c[otherBase]].size()) &&
                !help::shareAnEdge(faces[c[baseFace]], faces[c[otherBase]])
            )
        )
            isPrism = true;

        if( isPrism )
            return true;

        return false;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::analyseLayers()
{
    Info << "Analysing mesh for bnd layer existence" << endl;

    //- find layers in patch
    labelListPMG cellInLayer;
    const label nGroups =
        help::groupMarking
        (
            cellInLayer,
            meshBndLayerNeighbourOperator(mesh_, surfaceEngine()),
            meshBndLayerSelectorOperator(mesh_)
        );

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& faceOwner = mse.faceOwners();

    const PtrList<writePatch>& boundaries = mesh_.boundaries();

    //- create patch name to index addressing
    std::map<word, label> patchNameToIndex;
    forAll(boundaries, patchI)
        patchNameToIndex[boundaries[patchI].patchName()] = patchI;

    //- check layer labels over a patch
    List<DynList<label> > groupsAtPatch(boundaries.size());
    forAll(facePatch, bfI)
    {
        const label cellI = faceOwner[bfI];

        if( cellInLayer[cellI] < 0 )
            continue;

        groupsAtPatch[facePatch[bfI]].appendIfNotIn(cellInLayer[cellI]);
    }

    //- set the information which patches have an extruded layer
    labelList groupIDs(nGroups, -1);

    layerAtPatch_.setSize(boundaries.size());
    layerAtPatch_ = -1;

    label nValidLayers(0);
    forAll(groupsAtPatch, patchI)
    {
        if( groupsAtPatch[patchI].size() == 1 )
        {
            const label groupI = groupsAtPatch[patchI][0];

            if( groupIDs[groupI] == -1 )
                groupIDs[groupI] = nValidLayers++;

            layerAtPatch_[patchI] = groupIDs[groupI];
        }
    }

    //- set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nValidLayers);
    forAll(layerAtPatch_, patchI)
    {
        if( layerAtPatch_[patchI] < 0 )
            continue;

        patchesInLayer_[layerAtPatch_[patchI]].append
        (
            boundaries[patchI].patchName()
        );
    }

    //- set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(boundaries, patchI)
    {
        const label layerI = layerAtPatch_[patchI];

        if( layerI < 0 )
            continue;

        label maxNumLayers(1);
        bool found(false);
        forAll(patchesInLayer_[layerI], lpI)
        {
            const word pName = patchesInLayer_[layerI][lpI];

            std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if( it != numLayersForPatch_.end() )
            {
                found = true;

                //- check if the layer is interrupted at this patch
                if(
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the numbe of layers and lock this location
                    nLayersAtPatch[patchNameToIndex[pName]] = it->second;
                    protectedValue[patchNameToIndex[pName]] = true;
                }
                else
                {
                    //- take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it->second);
                }
            }
        }

        if( !found )
            maxNumLayers = globalNumLayers_;

        //- set the number of layer to all patches which are not protected
        forAll(patchesInLayer_[layerI], lpI)
        {
            const label ptchI = patchNameToIndex[patchesInLayer_[layerI][lpI]];

            if( !protectedValue[ptchI] )
                nLayersAtPatch[ptchI] = maxNumLayers;
        }
    }

    Info << "nLayersAtPatch " << nLayersAtPatch << endl;

    //- set the number of boundary layers which shall be generated above
    //- each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = 0;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        if( layerAtPatch_[patchI] < 0 )
            continue;

        if( nLayersAtPatch[patchI] < 0 )
        {
            nLayersAtBndFace_[bfI] = globalNumLayers_;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];
        }
    }

    Info << "nLayersAtBndFace_ " << nLayersAtBndFace_ << endl;
}

void refineBoundaryLayers::calculateAddressing
(
    const label bfI,
    label& baseFace,
    DynList<edge, 48>& edges,
    DynList<DynList<label, 2>, 48>& edgeFaces,
    DynList<DynList<label, 10>, 24>& faceEdges
) const
{
    const label nInternalFaces = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh_.faces();
    const cell& c = mesh_.cells()[faceOwner[bfI]];

    faceEdges.setSize(c.size());
    baseFace = -1;
    forAll(c, fI)
    {
        if( c[fI] - nInternalFaces == bfI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }
}

bool refineBoundaryLayers::findHairsForFace
(
    const label bfI,
    DynList<edge>& hairEdges
) const
{
    const label nInternalFaces = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh_.faces();
    const cell& c = mesh_.cells()[faceOwner[bfI]];

    //- check cell topology
    DynList<edge, 48> edges;
    DynList<DynList<label, 2>, 48> edgeFaces;
    DynList<DynList<label, 10>, 24> faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if( c[fI] - nInternalFaces == bfI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }

    if( (baseFace < 0) || ((c.size() - faces[c[baseFace]].size()) != 2) )
        return false;

    //- check if all faces attached to the base face are quads
    bool isPrism(true);

    const face& bf = faces[c[baseFace]];
    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][(pI+bf.size()-1)%bf.size()];

        if( edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2 )
        {
            isPrism = false;
            break;
        }

        //- find the face attached to the edge after the current point
        label otherNextFace = edgeFaces[nextEdge][0];
        if( otherNextFace == baseFace )
            otherNextFace = edgeFaces[nextEdge][1];

        //- find the face attached to the edge before the current point
        label otherPrevFace = edgeFaces[prevEdge][0];
        if( otherPrevFace == baseFace )
            otherPrevFace = edgeFaces[prevEdge][1];

        label commonEdge;
        for(commonEdge=0;commonEdge<edges.size();++commonEdge)
            if(
                edgeFaces[commonEdge].contains(otherNextFace) &&
                edgeFaces[commonEdge].contains(otherPrevFace)
            )
                break;

        if( commonEdge == edges.size() )
        {
            isPrism = false;
            break;
        }

        //- there exists a common edge which shall be used as a hair
        if( edges[commonEdge].start() == bf[pI] )
        {
            hairEdges.append(edges[commonEdge]);
        }
        else
        {
            hairEdges.append(edges[commonEdge].reverseEdge());
        }
    }

    return isPrism;
}

bool refineBoundaryLayers::findSplitEdges()
{
    bool validLayer(true);

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& bp = mse.bp();

    # ifdef USE_OMP
    # pragma omp parallel if( bFaces.size() > 1000 )
    # endif
    {
        edgeListPMG localEdges;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(bFaces, bfI)
        {
            if( layerAtPatch_[facePatch[bfI]] < 0 )
                continue;

            //- find hair edges for this face
            DynList<edge> hairEdges;
            if( !findHairsForFace(bfI, hairEdges) )
            {
                validLayer = false;
                continue;
            }

            const face& bf = bFaces[bfI];
            forAll(bf, pI)
            {
                //- check if every the hair shall be store or not
                //- only a hair edge from a face with the smallest label
                //- out of all faces at a points is stored
                const label bpI = bp[bf[pI]];

                bool store(true);
                forAllRow(pFaces, bpI, pfI)
                {
                    const face& obf = bFaces[pFaces(bpI, pfI)];

                    if(
                        (obf.which(hairEdges[pI].end()) < 0) &&
                        (pFaces(bpI, pfI) < bfI)
                    )
                    {
                        store = false;
                        break;
                    }
                }

                if( store )
                {
                    //- hair edge shall be stored
                    localEdges.append(hairEdges[pI]);
                }
            }
        }

        # ifdef USE_OMP
        //- find the starting element for this thread
        label startEl;
        # pragma omp critical
        {
            startEl = splitEdges_.size();

            splitEdges_.setSize(startEl+localEdges);
        }

        //- copy the local data to splitEdges_
        forAll(localEdges, i)
            splitEdges_[startEl++] = localEdges[i];
        # else
        //- just transfer the data to splitEdges_
        splitEdges_.transfer(localEdges);
        # endif
    }

    //- create point to split edges addressing
    splitEdgesAtPoint_.reverseAddressing(splitEdges_);

    reduce(validLayer, minOp<bool>());

    return validLayer;
}

void refineBoundaryLayers::generateNewVertices()
{
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& pointPatches = mse.pointPatches();
    const labelList& bp = mse.bp();

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelListPMG nNodesAtEdge(splitEdges_.size());

    //- count the number of vertices for each split edge
    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # else
    const label nThreads = 1;
    # endif

    DynList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];
        nPoints = 0;

        //- start counting vertices at each thread
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            //- get the requested number of boundary layers
            label nLayers(globalNumLayers_);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);

            const label bpI = bp[e.start()];

            forAllRow(pointPatches, bpI, ppI)
            {
                const word& patchName =
                    boundaries[pointPatches(bpI, ppI)].patchName();

                //- overrride the global value with the maximum number of layers
                //- at this edge
                const std::map<word, label>::const_iterator it =
                    numLayersForPatch_.find(patchName);
                if( it != numLayersForPatch_.end() )
                    nLayers = Foam::max(nLayers, it->second);

                //- override with the maximum ratio
                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if( rIt != thicknessRatioForPatch_.end() )
                    ratio = rIt->second;

                //- override with the minimum thickness set for this edge
                const std::map<word, scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if( tIt != maxThicknessForPatch_.end() )
                    thickness = Foam::min(thickness, tIt->second);
            }

            //- store the information
            firstLayerThickness[seI] = thickness;
            thicknessRatio[seI] = ratio;
            nNodesAtEdge[seI] = nLayers + 1;
            nPoints += nLayers - 1;
        }
    }

    if( Pstream::parRun() )
    {
        //- transfer the information over all processor for edges
        //- at inter-processor boundaries
        const labelListPMG& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const VRWGraph& edgeAtProcs = mesh_.addressingData().edgeAtProcs();
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& neiProcs = mesh_.addressingData().edgeNeiProcs();
        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();

        //- exchange point number of layers
        std::map<label, LongList<labelPair> > exchangeNumLayers;
        std::map<label, LongList<labelledScalar> > exchangeThickness;
        std::map<label, LongList<labelledScalar> > exchangeRatio;
        forAll(neiProcs, i)
        {
            exchangeNumLayers.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeThickness.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
            exchangeRatio.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
        }

        //- exchange the number of layers
        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            const label s = se.start();
            label edgeI(-1);
            forAllRow(pointEdges, s, peI)
            {
                const label eI = pointEdges(s, peI);

                if( edges[eI] == se )
                {
                    edgeI = eI;
                    break;
                }
            }

            const label geI = globalEdgeLabel[edgeI];

            if( globalToLocal.found(geI) )
            {
                forAllRow(edgeAtProcs, edgeI, i)
                {
                    const label neiProc = edgeAtProcs(edgeI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeNumLayers[neiProc].append
                    (
                        labelPair(geI, nNodesAtEdge[seI])
                    );
                    exchangeThickness[neiProc].append
                    (
                        labelledScalar(geI, firstLayerThickness[seI])
                    );
                    exchangeRatio[neiProc].append
                    (
                        labelledScalar(geI, thicknessRatio[seI])
                    );
                }
            }
        }

        //- exchange number of layers
        LongList<labelPair> receivedNumLayers;
        help::exchangeMap(exchangeNumLayers, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            const label eI = globalToLocal[lp.first()];
            nNodesAtEdge[eI] = std::max(nNodesAtEdge[eI], lp.second());
        }

        //- exchange thickness ratio
        LongList<labelledScalar> receivedScalar;
        help::exchangeMap(exchangeRatio, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            const label eI = globalToLocal[ls.scalarLabel()];
            thicknessRatio[eI] = std::max(thicknessRatio[eI], ls.value());
        }

        //- exchange ,aximum thickness of the first layer
        receivedScalar.clear();
        help::exchangeMap(exchangeThickness, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            const label eI = globalToLocal[ls.scalarLabel()];
            firstLayerThickness[eI] =
                std::min(firstLayerThickness[eI], ls.value());
        }
    }

    //- allocate the memory
    newVerticesForSplitEdge_.setSizeAndRowSize(nNodesAtEdge);

    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    points.setSize(numPoints);

    Info << "Generating split vertices" << endl;

    //- generate vertices at split edges
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;

            scalar firstThickness = magv / nLayers;
            if( thicknessRatio[seI] > (1. + SMALL) )
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayers)) /
                        (1.0 - thicknessRatio[seI])
                    );
            }

            firstThickness = std::min(firstLayerThickness[seI], firstThickness);
            firstThickness /= (magv + VSMALL);

            //- generate vertices for this edge
            newVerticesForSplitEdge_(seI, 0) = e.start();

            scalar param = firstThickness;

            for(label pI=1;pI<nLayers;++pI)
            {
                //- generate the new vertex
                const point newP = points[e.start()] + param * v;

                param += firstThickness * Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_(seI, pI) = nPoints;
                points[nPoints++] = newP;
            }

            newVerticesForSplitEdge_(seI, nLayers) = e.end();
        }
    }

    Info << "Finished generating vertices at split edges" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
