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

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"
#include "helperFunctions.H"
#include "labelledScalar.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "OFstream.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::writeVTK
(
    const fileName& fName,
    const pointField& origin,
    const vectorField& vecs
)
{
    if( origin.size() != vecs.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeVTK(const fileName&,"
            " const pointField&, const vectorField&)"
        ) << "Sizes do not match" << abort(FatalError);

    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*origin.size() << " float\n";
    forAll(origin, pI)
    {
        const point& p = origin[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = p + vecs[pI];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << vecs.size()
         << " " << 3*vecs.size() << nl;
    forAll(vecs, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType,
    const vectorField& vecs
) const
{
    if( vecs.size() != hairEdges_.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeHairEdges"
            "(const fileName&, const direction, const vectorField&) const"
        ) << "Sizes do not match" << abort(FatalError);

    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField copyVecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            copyVecs[counter] = vecs[heI] * he.mag(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName, pts, copyVecs);
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType
) const
{
    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField vecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            vecs[counter] = he.vec(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName,pts, vecs);
}

const meshSurfaceEngine& boundaryLayerOptimisation::meshSurface() const
{
    if( !meshSurfacePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfaceEngine&"
                " boundaryLayerOptimisation::meshSurface()"
            ) << "Cannot generate meshSurfaceEngine" << abort(FatalError);
        # endif

        meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
    }

    return *meshSurfacePtr_;
}

const meshSurfacePartitioner&
boundaryLayerOptimisation::surfacePartitioner() const
{
    if( !partitionerPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfacePartitioner& "
                "boundaryLayerOptimisation::surfacePartitioner()"
            ) << "Cannot generate meshSurfacePartitioner" << abort(FatalError);
        # endif

        partitionerPtr_ = new meshSurfacePartitioner(meshSurface());
    }

    return *partitionerPtr_;
}

void boundaryLayerOptimisation::calculateNormalVectors
(
    const direction eType,
    pointNormalsType& pointPatchNormal
) const
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const vectorField& fNormals = mse.faceNormals();

    //- calculate point normals with respect to all patches at a point
    pointPatchNormal.clear();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( !(hairEdgeType_[hairEdgeI] & eType) )
            continue;

        const label bpI = bp[hairEdges_[hairEdgeI][0]];

        //- create an entry in a map
        patchNormalType* patchNormalPtr(NULL);
        # ifdef USE_OMP
        # pragma omp critical
            patchNormalPtr = &pointPatchNormal[bpI];
        # else
        patchNormalPtr = &pointPatchNormal[bpI];
        # endif
        patchNormalType& patchNormal = *patchNormalPtr;

        //- sum normals of faces attached to a point
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI].first = fNormals[bfI];
                patchNormal[patchI].second = mag(fNormals[bfI]);
            }
            else
            {
                patchNormal[patchI].first += fNormals[bfI];
                patchNormal[patchI].second += mag(fNormals[bfI]);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- gather information about face normals on other processors
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, LongList<refLabelledPointScalar> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const patchNormalType& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(patchNormalType, patchNormal, pIt)
                        exchangeData[neiProc].append
                        (
                            refLabelledPointScalar
                            (
                                it.key(),
                                labelledPointScalar
                                (
                                    pIt->first,
                                    pIt->second.first,
                                    pIt->second.second
                                )
                            )
                        );
                }
            }
        }

        LongList<refLabelledPointScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const refLabelledPointScalar& rlps = receivedData[i];
            const label bpI = globalToLocal[rlps.objectLabel()];

            patchNormalType& patchNormal = pointPatchNormal[bpI];

            const labelledPointScalar& lps = rlps.lps();
            patchNormal[lps.pointLabel()].first += lps.coordinates();
            patchNormal[lps.pointLabel()].second += lps.scalarValue();
        }
    }

    //- finally, calculate normal vectors
    # ifdef USE_OMP
    # pragma omp parallel
    # pragma omp single nowait
    # endif
    forAllIter(pointNormalsType, pointPatchNormal, it)
    {
        # ifdef USE_OMP
        # pragma omp task firstprivate(it)
        {
        # endif

        patchNormalType& patchNormal = it->second;

        forAllIter(patchNormalType, patchNormal, pIt)
        {
            pIt->second.first /= pIt->second.second;
            pIt->second.first /= (mag(pIt->second.first) + VSMALL);
        }

        # ifdef USE_OMP
        }
        # endif
    }
}

void boundaryLayerOptimisation::calculateHairEdges()
{
    const meshSurfaceEngine& mse = meshSurface();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const labelList& faceOwner = mse.faceOwners();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();

    const meshSurfacePartitioner& mPart = surfacePartitioner();

    //- detect layers in the mesh
    const detectBoundaryLayers detectLayers(mPart);

    hairEdges_ = detectLayers.hairEdges();
    hairEdgesAtBndPoint_ = detectLayers.hairEdgesAtBndPoint();

    //- mark boundary faces which are base face for the boundary layer
    const labelList& layerAtBndFace = detectLayers.faceInLayer();
    isBndLayerBase_.setSize(bFaces.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(layerAtBndFace, bfI)
    {
        if( layerAtBndFace[bfI] < 0 )
        {
            isBndLayerBase_[bfI] = false;
        }
        else
        {
            isBndLayerBase_[bfI] = true;
        }
    }

    # ifdef DEBUGLayer
    const label bndLayerFaceId = mesh_.addFaceSubset("bndLayerFaces");
    const label startBndFaceI = mesh_.boundaries()[0].patchStart();
    forAll(isBndLayerBase_, bfI)
        if( isBndLayerBase_[bfI] )
            mesh_.addFaceToSubset(bndLayerFaceId, startBndFaceI+bfI);
    # endif

    //- check if a face is an exiting face for a bnd layer
    isExitFace_.setSize(isBndLayerBase_.size());
    isExitFace_ = false;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeFaces, edgeI)
    {
        //- avoid edges at inter-processor boundaries
        if( edgeFaces.sizeOfRow(edgeI) != 2 )
            continue;

        const label f0 = edgeFaces(edgeI, 0);
        const label f1 = edgeFaces(edgeI, 1);

        //- both faces have to be part of the same cell
        if( faceOwner[f0] != faceOwner[f1] )
            continue;

        //- check if the feature edge is attachd to exittin faces
        if
        (
            (isBndLayerBase_[f0] && (bFaces[f1].size() == 4)) &&
            (isBndLayerBase_[f1] && (bFaces[f0].size() == 4))
        )
        {
            isExitFace_[f0] = true;
            isExitFace_[f1] = true;
        }
    }

    # ifdef DEBUGLayer
    const label exittingFaceId = mesh_.addFaceSubset("exittingFaces");
    forAll(isExitFace_, bfI)
        if( isExitFace_[bfI] )
            mesh_.addFaceToSubset(exittingFaceId, startBndFaceI+bfI);
    # endif

    //- classify hair edges as they require different tretment
    //- in the smoothing procedure
    hairEdgeType_.setSize(hairEdges_.size());

    const labelHashSet& corners = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const labelHashSet& featureEdges = mPart.featureEdges();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        hairEdgeType_[hairEdgeI] = NONE;

        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        //- check if it is a boundary edge
        forAllRow(bpEdges, bpI, peI)
        {
            const label beI = bpEdges(bpI, peI);

            const edge& be = edges[bpEdges(bpI, peI)];

            if( be == e )
            {
                hairEdgeType_[hairEdgeI] |= BOUNDARY;

                if( featureEdges.found(beI) )
                    hairEdgeType_[hairEdgeI] |= FEATUREEDGE;
            }

            if( corners.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATCORNER;
            }
            else if( edgePoints.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATEDGE;
            }
        }

        if( !(hairEdgeType_[hairEdgeI] & BOUNDARY) )
            hairEdgeType_[hairEdgeI] |= INSIDE;
    }

    thinnedHairEdge_.setSize(hairEdges_.size());

    //- calculate which other hair edges influence a hair edges
    //- and store it in a graph
    hairEdgesNearHairEdge_.setSize(hairEdges_.size());

    const cellListPMG& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    VRWGraph bpFacesHelper(bpEdges.size());
    forAll(faceOwner, bfI)
    {
        const label cellI = faceOwner[bfI];

        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];
                if( bpI < 0 )
                    continue;

                bpFacesHelper.appendIfNotIn(bpI, c[fI]);
            }
        }
    }

    forAll(hairEdges_, hairEdgeI)
    {
        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        DynList<label> neiHairEdges;

        //- find mesh faces comprising of the current hair edge
        forAllRow(bpFacesHelper, bpI, pfI)
        {
            const face& f = faces[bpFacesHelper(bpI, pfI)];

            //- face must be a quad
            if( f.size() != 4 )
                continue;

            //- check if the current face comprises of the hair edge
            label faceEdge(-1);
            forAll(f, eI)
                if( f.faceEdge(eI) == e )
                {
                    faceEdge = eI;
                    break;
                }

            if( faceEdge != -1 )
            {
                //- check if the opposite edge is also a hair edge
                const label eJ = (faceEdge+2) % 4;

                const edge fe = f.faceEdge(eJ);

                for(label i=0;i<2;++i)
                {
                    const label bpJ = bp[fe[i]];

                    if( bpJ >= 0 )
                    {
                        forAllRow(hairEdgesAtBndPoint_, bpJ, pI)
                        {
                            const label heJ = hairEdgesAtBndPoint_(bpJ, pI);
                            if( hairEdges_[heJ] == fe )
                                neiHairEdges.append(heJ);
                        }
                    }
                }
            }
        }

        hairEdgesNearHairEdge_.setRow(hairEdgeI, neiHairEdges);
    }

    # ifdef DEBUGLayer
    const label hairEdgesId = mesh_.addPointSubset("hairEdgePoints");
    const label bndHairEdgeId = mesh_.addPointSubset("bndHairEdgePoints");
    const label featureHairEdgeId = mesh_.addPointSubset("featureEdgePoints");
    const label cornerHairEdgeId = mesh_.addPointSubset("cornerHairEdgePoints");
    const label hairEdgeAtEdgeId = mesh_.addPointSubset("hairEdgeAtEdgePoints");

    forAll(hairEdgeType_, heI)
    {
        const edge& e = hairEdges_[heI];

        mesh_.addPointToSubset(hairEdgesId, e.start());
        mesh_.addPointToSubset(hairEdgesId, e.end());
        if( hairEdgeType_[heI] & FEATUREEDGE)
        {
            mesh_.addPointToSubset(featureHairEdgeId, e.start());
            mesh_.addPointToSubset(featureHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & BOUNDARY)
        {
            mesh_.addPointToSubset(bndHairEdgeId, e.start());
            mesh_.addPointToSubset(bndHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATCORNER)
        {
            mesh_.addPointToSubset(cornerHairEdgeId, e.start());
            mesh_.addPointToSubset(cornerHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATEDGE)
        {
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.start());
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.end());
        }
    }

    mesh_.write();
    # endif
}

void boundaryLayerOptimisation::calculateHairVectorsAtTheBoundary
(
    vectorField& hairVecs
)
{
    //- set the size of hairVecs
    hairVecs.setSize(hairEdges_.size());

    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const direction hairType = hairEdgeType_[hairEdgeI];

        if( hairType & BOUNDARY )
        {
            const edge& he = hairEdges_[hairEdgeI];
            vector& hv = hairVecs[hairEdgeI];
            hv = vector::zero;

            if( hairType & FEATUREEDGE )
            {
                //- do not modify hair vectors at feature edges
                hv = he.vec(points);
            }
            else if( hairType & (ATEDGE|ATCORNER) )
            {
                //- this is a case of O-layer at a corner or feature edge

                //- find the surface edges corresponding to the hair edge
                label beI(-1);

                const label bps = bp[he.start()];
                forAllRow(bpEdges, bps, bpeI)
                {
                    const label beJ = bpEdges(bps, bpeI);

                    if( edges[beJ] == he )
                    {
                        beI = beJ;
                        continue;
                    }
                }

                if( beI < 0 )
                    FatalErrorIn
                    (
                        "boundaryLayerOptimisation::"
                        "calculateHairVectorsAtTheBoundary(vectorField&)"
                    ) << "Cannot find hair edge "
                      << hairEdgeI << abort(FatalError);

                //- find the vector at the same angle from both feature edges
                forAllRow(edgeFaces, beI, befI)
                {
                    const face& bf = bFaces[edgeFaces(beI, befI)];
                    const vector fNormal = bf.normal(points);

                    const label pos = bf.which(he.start());

                    if( pos < 0 )
                        FatalErrorIn
                        (
                            "boundaryLayerOptimisation::"
                            "calculateHairVectorsAtTheBoundary(vectorField&)"
                        ) << "Cannot find hair edge "
                          << hairEdgeI << " in face " << bf
                          << abort(FatalError);

                    if( he.end() == bf.prevLabel(pos) )
                    {
                        const edge fe = bf.faceEdge(pos);
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                    else
                    {
                        const edge fe = bf.faceEdge(bf.rcIndex(pos));
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                }
            }
            else
            {
                FatalErrorIn
                (
                    "boundaryLayerOptimisation::"
                    "calculateHairVectorsAtTheBoundary(vectorField&)"
                ) << "Invalid hair type " << label(hairType)
                  << abort(FatalError);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- collect data at inter-processor boundaries
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();

        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
        const VRWGraph& edgesAtProcs =
            mesh_.addressingData().edgeAtProcs();
        const labelLongList& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const Map<label>& globalToLocalEdge =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& eNeiProcs =
            mesh_.addressingData().edgeNeiProcs();

        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(eNeiProcs, i)
            exchangeData[eNeiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(hairEdgesAtBndPoint_, bpI, i)
            {
                const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                const edge& he = hairEdges_[hairEdgeI];

                const direction eType = hairEdgeType_[hairEdgeI];

                //- filte out edges which are not relevant
                if( !(eType & BOUNDARY) )
                    continue;

                forAllRow(pointEdges, he.start(), peI)
                {
                    const label edgeI = pointEdges(he.start(), peI);

                    if( he == edges[edgeI] )
                    {
                        labelledPoint lp
                        (
                            globalEdgeLabel[edgeI],
                            hairVecs[hairEdgeI]
                        );

                        forAllRow(edgesAtProcs, edgeI, j)
                        {
                            const label neiProc = edgesAtProcs(edgeI, j);

                            if( neiProc == Pstream::myProcNo() )
                                continue;

                            LongList<labelledPoint>& dts =
                                exchangeData[neiProc];

                            dts.append(lp);
                        }
                    }
                }
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPoint& lp = receivedData[i];
            const label edgeI = globalToLocalEdge[lp.pointLabel()];
            const edge& e = edges[edgeI];

            bool found(false);
            for(label pI=0;pI<2;++pI)
            {
                const label bpI = bp[e[pI]];

                if( bpI < 0 )
                    continue;

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    if( hairEdges_[hairEdgeI] == e )
                    {
                        hairVecs[hairEdgeI] += lp.coordinates();

                        found = true;
                        break;
                    }
                }

                if( found )
                    break;
            }
        }
    }

    //- calculate new normal vectors
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairVecs, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
            hairVecs[hairEdgeI] /= (mag(hairVecs[hairEdgeI]) + VSMALL);
    }

    # ifdef DEBUGLayer
    Info << "Saving bnd hair vectors" << endl;
    writeHairEdges("bndHairVectors.vtk", (BOUNDARY | ATEDGE), hairVecs);
    # endif
}

void boundaryLayerOptimisation::optimiseHairNormalsAtTheBoundary()
{
    if( nSmoothNormals_ == 0 )
        return;

    pointFieldPMG& points = mesh_.points();

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();

    //- calculate hair vectors
    //- they point in the normal direction to the surface
    vectorField hairVecs(hairEdges_.size());
    calculateHairVectorsAtTheBoundary(hairVecs);

    //- smooth the variation of normals to reduce the twisting of faces
    label nIter(0);

    while( nIter++ < nSmoothNormals_ );
    {
        vectorField newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            if( eType & BOUNDARY )
            {
                if( eType & (FEATUREEDGE | ATCORNER) )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal += hairVecs[hairEdgeI];
                }
                else if( eType & ATEDGE )
                {
                    const edge& he = hairEdges_[hairEdgeI];

                    DynList<label, 2> edgeFaces;
                    const label bps = bp[he.start()];
                    forAllRow(pFaces, bps, pfI)
                    {
                        const label bfI = pFaces(bps, pfI);
                        const face& bf = bFaces[bfI];

                        forAll(bf, eI)
                        {
                            if( bf.faceEdge(eI) == he )
                            {
                                edgeFaces.append(bfI);
                            }
                        }
                    }

                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        //- check if the neighbour hair edge shares a boundary
                        //- face with the current hair edge
                        bool useNeighbour(false);
                        const edge& nhe = hairEdges_[hairEdgeJ];
                        forAll(edgeFaces, efI)
                        {
                            const face& bf = bFaces[edgeFaces[efI]];

                            forAll(bf, eI)
                            {
                                if( bf.faceEdge(eI) == nhe )
                                {
                                    useNeighbour = true;
                                    break;
                                }
                            }
                        }

                        if( useNeighbour )
                            newNormal += hairVecs[hairEdgeJ];
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsAtTheBoundary(const label)"
                    ) << "Cannot smooth hair with type " << label(eType)
                      << abort(FatalError);
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                    const edge& he = hairEdges_[hairEdgeI];

                    const direction eType = hairEdgeType_[hairEdgeI];

                    //- filte out edges which are not relevant
                    if( !(eType & BOUNDARY) )
                        continue;
                    if( !(eType & ATEDGE) )
                        continue;
                    if( eType & (FEATUREEDGE|ATCORNER) )
                        continue;

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledPoint lp
                            (
                                globalEdgeLabel[edgeI],
                                hairVecs[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledPoint>& dts =
                                    exchangeData[neiProc];

                                dts.append(lp);
                            }
                        }
                    }
                }
            }

            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledPoint& lp = receivedData[i];
                const label edgeI = globalToLocalEdge[lp.pointLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            hairVecs[hairEdgeI] += lp.coordinates();

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, heI)
        {
            if( hairEdgeType_[heI] & BOUNDARY )
            {
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
                newNormals[heI] = 0.5 * (newNormals[heI] + hairVecs[heI]);
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
            }
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "bndHairVectors_"+help::scalarToText(nIter)+".vtk",
                (BOUNDARY | ATEDGE),
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
        {
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            points[he.end()] = points[he.start()] + hv * he.mag(points);
        }
    }
}

void boundaryLayerOptimisation::optimiseHairNormalsInside()
{
    if( nSmoothNormals_ == 0 )
        return;

    pointFieldPMG& points = mesh_.points();

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();

    //- they point in the normal direction to the surface
    vectorField hairVecs(hairEdges_.size());

    if( reCalculateNormals_ )
    {
        //- calculate point normals with respect to all patches at a point
        pointNormalsType pointPatchNormal;
        calculateNormalVectors(INSIDE, pointPatchNormal);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            const direction hairType = hairEdgeType_[hairEdgeI];

            if( hairType & INSIDE )
            {
                vector& hv = hairVecs[hairEdgeI];
                hv = vector::zero;

                const label bpI = bp[hairEdges_[hairEdgeI].start()];

                label counter(0);
                const patchNormalType& patchNormals = pointPatchNormal[bpI];
                forAllConstIter(patchNormalType, patchNormals, pIt)
                {
                    hv -= pIt->second.first;
                    ++counter;
                }

                if( counter == 0 )
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsInside()"
                    ) << "No valid patches for boundary point "
                      << bp[hairEdges_[hairEdgeI].start()] << abort(FatalError);
                }

                hv /= (mag(hv) + VSMALL);
            }
            else
            {
                //- initialise boundary hair vectors. They influence internal
                //- hairs connected to them
                vector hvec = hairEdges_[hairEdgeI].vec(points);
                hvec /= (mag(hvec) + VSMALL);
                hairVecs[hairEdgeI] = hvec;
            }
        }
    }
    else
    {
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            //- initialise boundary hair vectors.
            vector hvec = hairEdges_[hairEdgeI].vec(points);
            hvec /= (mag(hvec) + VSMALL);
            hairVecs[hairEdgeI] = hvec;
        }
    }

    # ifdef DEBUGLayer
    writeHairEdges("insideHairVectors.vtk", (INSIDE|BOUNDARY), hairVecs);
    # endif

    //- smooth the variation of normals to reduce twisting of faces
    label nIter(0);

    while( nIter++ < nSmoothNormals_ )
    {
        vectorField newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            const edge& he = hairEdges_[hairEdgeI];
            const vector& heVec = hairVecs[hairEdgeI];

            if( eType & INSIDE )
            {
                if( eType & ATCORNER )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal = hairVecs[hairEdgeI];
                }
                else if( eType & ATEDGE )
                {
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        const edge& nhe = hairEdges_[hairEdgeJ];
                        const vector& nhVec = hairVecs[hairEdgeJ];

                        vector n = nhVec ^ (points[nhe[0]] - points[he[0]]);
                        n /= (mag(n) + VSMALL);

                        vector newVec = heVec - (heVec & n) * n;
                        newVec /= (mag(newVec) + VSMALL);

                        scalar weight = 1.0;

                        if( Pstream::parRun() )
                        {
                            //- edges at inter-processor boundaries contribute
                            //- at two sides are given weight 0.5
                            const edge be(he[0], nhe[0]);
                            const label bpI = bp[he[0]];

                            forAllRow(bpEdges, bpI, bpeI)
                            {
                                const edge& bndEdge = edges[bpEdges(bpI, bpeI)];

                                if( bndEdge == be )
                                {
                                    weight = 0.5;
                                    break;
                                }
                            }
                        }

                        newNormal += weight * newVec;
                    }
                }
                else
                {
                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        scalar weight = 1.0;

                        if( Pstream::parRun() )
                        {
                            //- edges at inter-processor boundaries contribute
                            //- at two sides are given weight 0.5
                            const edge& nhe = hairEdges_[hairEdgeJ];
                            const edge be(he[0], nhe[0]);
                            const label bpI = bp[he[0]];

                            forAllRow(bpEdges, bpI, bpeI)
                            {
                                const edge& bndEdge = edges[bpEdges(bpI, bpeI)];

                                if( bndEdge == be )
                                {
                                    weight = 0.5;
                                    break;
                                }
                            }
                        }

                        newNormal += weight * hairVecs[hairEdgeJ];
                    }
                }
            }
            else
            {
                //- copy the existing hair vector
                newNormal = hairVecs[hairEdgeI];
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                    const edge& he = hairEdges_[hairEdgeI];

                    const direction eType = hairEdgeType_[hairEdgeI];

                    //- handle boundary points, only
                    if( !(eType & INSIDE) || (eType & ATCORNER) )
                        continue;

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledPoint lp
                            (
                                globalEdgeLabel[edgeI],
                                hairVecs[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledPoint>& dts =
                                    exchangeData[neiProc];

                                dts.append(lp);
                            }
                        }
                    }
                }
            }

            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledPoint& lp = receivedData[i];
                const label edgeI = globalToLocalEdge[lp.pointLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            hairVecs[hairEdgeI] += lp.coordinates();

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, hairEdgeI)
        {
            if( hairEdgeType_[hairEdgeI] & INSIDE )
                newNormals[hairEdgeI] /= (mag(newNormals[hairEdgeI]) + VSMALL);
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "insideHairVectors_"+help::scalarToText(nIter)+".vtk",
                INSIDE,
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & INSIDE )
        {
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            points[he.end()] = points[he.start()] + hv * he.mag(points);
        }
    }
}

scalar boundaryLayerOptimisation::calculateThickness
(
    const label heI,
    const label heJ
) const
{
    const pointFieldPMG& points = mesh_.points();

    //- references to hair edges
    const edge& he = hairEdges_[heI];
    const edge& nhe = hairEdges_[heJ];

    //- distance vector between the surface points of hair edges
    const point& sp = points[he[0]];
    const point& ep = points[nhe[0]];
    const vector dv = ep - sp;
    const scalar magDv = mag(dv);

    //- calculate layer thickness
    const scalar currThickness = he.mag(points);
    scalar retThickness = currThickness;

    const scalar currNeiThickness = nhe.mag(points);
    scalar suggestedNeiThickness = currNeiThickness;

    //- calculate layer height at the current point
    const point npAlpha = help::nearestPointOnTheEdge(sp, ep, points[he[1]]);
    const scalar currHeight = mag(npAlpha - points[he[1]]);
    scalar retHeight = currHeight;
    const scalar cosAlpha = sign((npAlpha - sp) & dv) * mag(npAlpha - sp);
    const scalar alpha =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosAlpha / (currThickness + VSMALL))
            )
        );

    //- calculate the height of the layer at the neighbour
    //- point
    const point npBeta = help::nearestPointOnTheEdge(ep, sp, points[nhe[1]]);
    const scalar currNeiHeight = mag(npBeta - points[nhe[1]]);
    scalar suggestedNeiHeight = currNeiHeight;
    const scalar cosBeta = sign((npBeta - ep) & -dv) * mag(npBeta - ep);
    const scalar beta =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosBeta / (currNeiThickness + VSMALL))
            )
        );

    //- check if the current thickness is Ok for the local curvature
    if( (alpha + beta) < M_PI )
    {
        const scalar gamma = M_PI - (alpha + beta);
        const scalar sinGamma = Foam::max(SMALL, Foam::sin(gamma));
        const scalar sinAlpha = Foam::max(SMALL, Foam::sin(alpha));
        const scalar sinBeta = Foam::max(SMALL, Foam::sin(beta));

        //- max allowed thickness and layer height due to curvature
        retThickness =
            Foam::min
            (
                retThickness,
                featureSizeFactor_ * magDv * sinBeta / sinGamma
            );
        retHeight *= (retThickness / (currThickness + VSMALL));

        //- max allowed neighbour hair thickness
        //- and layer height due to curvature
        suggestedNeiThickness =
            Foam::min
            (
                suggestedNeiThickness,
                featureSizeFactor_ * magDv * sinAlpha / sinGamma
            );
        suggestedNeiHeight *=
            (suggestedNeiThickness / (currNeiThickness + VSMALL));
    }

    //- check the height variation
    const scalar tanVal = (retHeight - suggestedNeiHeight) / (magDv + VSMALL);

    if( tanVal > relThicknessTol_ )
    {
        retHeight = suggestedNeiHeight + relThicknessTol_ * magDv;

        retThickness = (retHeight / currHeight) * currThickness;
    }

    return retThickness;
}

scalar boundaryLayerOptimisation::calculateThicknessOverCell
(
    const label heI,
    const label cellI,
    const label baseFaceI
) const
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

    const face& bf = faces[baseFaceI];

    const edge& he = hairEdges_[heI];

    const point& sp = points[he[0]];
    const point& ep = points[he[1]];

    scalar maxThickness = he.mag(points);

    //- the base face must not contain the hair edge
    //- this is the case at exitting layers
    forAll(bf, eI)
        if( bf.faceEdge(eI) == he )
            return maxThickness;

    forAll(c, fI)
    {
        if( c[fI] == baseFaceI )
            continue;

        const face& f = faces[c[fI]];

        if( help::shareAnEdge(bf, f) && (f.which(he.start()) == -1) )
        {
            point intersection;

            if( !help::lineFaceIntersection(sp, ep, f, points, intersection) )
                continue;

            const scalar maxDist = featureSizeFactor_ * mag(intersection - sp);

            maxThickness =
                Foam::min(maxThickness, maxDist);
        }
    }

    return maxThickness;
}

void boundaryLayerOptimisation::optimiseThicknessVariation
(
    const direction edgeType
)
{
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const label start = mesh_.nInternalFaces();
    const labelList& faceOwner = mse.faceOwners();

    vectorField hairDirections(hairEdges_.size());
    scalarField hairLength(hairEdges_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        vector n = hairEdges_[hairEdgeI].vec(points);

        hairLength[hairEdgeI] = (Foam::mag(n) + VSMALL);
        hairDirections[hairEdgeI] = n / hairLength[hairEdgeI];
    }

    //- check if the hair edge intersects some other face in the cells
    //- attached to the hair edge
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & edgeType )
        {
            const label bpI = bp[hairEdges_[hairEdgeI].start()];

            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);
                const label baseFaceI = start + bfI;
                const label cOwn = faceOwner[bfI];

                const scalar maxThickness =
                    calculateThicknessOverCell
                    (
                        hairEdgeI,
                        cOwn,
                        baseFaceI
                    );

                if( hairLength[hairEdgeI] > maxThickness )
                {
                    //- make the hair edge shorter
                    hairLength[hairEdgeI] = maxThickness;

                    thinnedHairEdge_[hairEdgeI] = true;
                }
            }
        }
    }

    //- reduce thickness of the layer
    //- such that the variation of layer thickness
    //- It is an iterative process where the layer is thinned in the regions
    //- where the tangent is greater than the tolerance value or the curvature
    //- permits thicker boundary layers.
    boolList activeHairEdge(hairEdges_.size(), true);
    bool changed;
    label nIter(0);
    do
    {
        changed = false;

        boolList modifiedHairEdge(hairEdges_.size(), false);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            const scalar magN = hairLength[hairEdgeI];

            if( magN < VSMALL )
                FatalErrorIn
                (
                    "void boundaryLayerOptimisation::optimiseThicknessVariation"
                    "(const direction, const label, const scalar, const scalar)"
                ) << "Zero layer thickness at hair edge " << hairEdgeI
                  << ". Exitting..." << exit(FatalError);

            if( hairEdgeType_[hairEdgeI] & edgeType )
            {
                forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                {
                    const label hairEdgeJ =
                        hairEdgesNearHairEdge_(hairEdgeI, nheI);

                    if( !activeHairEdge[hairEdgeJ] )
                        continue;

                    const scalar maxThickness =
                        calculateThickness
                        (
                            hairEdgeI,
                            hairEdgeJ
                        );

                    if( hairLength[hairEdgeI] > maxThickness )
                    {
                        //- make the hair edge shorter
                        hairLength[hairEdgeI] = maxThickness;

                        changed = true;
                        thinnedHairEdge_[hairEdgeI] = true;
                        modifiedHairEdge[hairEdgeI] = true;
                    }
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledScalar> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    if( !(hairEdgeType_[hairEdgeI] & edgeType) )
                        continue;

                    const edge& he = hairEdges_[hairEdgeI];

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledScalar lScalar
                            (
                                globalEdgeLabel[edgeI],
                                hairLength[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledScalar>& dts =
                                    exchangeData[neiProc];

                                dts.append(lScalar);
                            }
                        }
                    }
                }
            }

            LongList<labelledScalar> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledScalar& lScalar = receivedData[i];
                const label edgeI = globalToLocalEdge[lScalar.scalarLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            if( lScalar.value() < hairLength[hairEdgeI] )
                            {
                                hairLength[hairEdgeI] = lScalar.value();
                                changed = true;
                                thinnedHairEdge_[hairEdgeI] = true;
                                modifiedHairEdge[hairEdgeI] = true;
                            }

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- reduce the information over all processors
        reduce(changed, maxOp<bool>());

        if( !changed )
            break;

        //- move boundary vertices to the new positions
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            if( hairEdgeType_[hairEdgeI] & edgeType )
            {
                const edge& he = hairEdges_[hairEdgeI];
                const vector& hv = hairDirections[hairEdgeI];
                const point& s = points[he.start()];

                points[he.end()] = s + hairLength[hairEdgeI] * hv;

                activeHairEdge[hairEdgeI] = modifiedHairEdge[hairEdgeI];
            }
        }
    } while( changed && (++nIter < 1000) );
}

bool boundaryLayerOptimisation::optimiseLayersAtExittingFaces()
{
    bool modified(false);

    //- find edge points inside the mesh with more than one hair edge
    //- attached to it
    labelList nEdgesAtPoint(mesh_.points().size(), 0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, heI)
    {
        # ifdef USE_OMP
        # pragma omp atomic
        # endif
        ++nEdgesAtPoint[hairEdges_[heI].end()];
    }

    //- find the points with more than one hair edge which was modified
    //- in the previous procedure
    boolList thinnedPoints(mesh_.points().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(thinnedHairEdge_, heI)
    {
        if
        (
            thinnedHairEdge_[heI] &&
            (nEdgesAtPoint[hairEdges_[heI].end()] > 2)
        )
        {
            modified = true;
            thinnedPoints[hairEdges_[heI].end()] = true;
        }
    }

    reduce(modified, maxOp<bool>());

    if( !modified )
        return false;

    Info << "Hair edges at exitting faces shall "
         << "be modified due to inner constraints" << endl;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::optimiseLayer()
{
    //- create surface smoother
    meshSurfaceOptimizer surfOpt(meshSurface());

    //- lock exitting faces and feature edges
    labelLongList lockedFaces;
    forAll(isExitFace_, bfI)
        if( isExitFace_[bfI] )
            lockedFaces.append(bfI);
    surfOpt.lockBoundaryFaces(lockedFaces);
    surfOpt.lockFeatureEdges();

    label nIter(0);
    do
    {
        thinnedHairEdge_ = false;

        //- calculate normals at the boundary
        optimiseHairNormalsAtTheBoundary();

        //- smoothing thickness variation of boundary hairs
        optimiseThicknessVariation(BOUNDARY);

        if( true )
        {
            meshSurfaceEngineModifier bMod(meshSurface());
            bMod.updateGeometry();

            surfOpt.optimizeSurface(2);
            bMod.updateGeometry();
        }

        # ifdef DEBUGLayer
        label counter(0);
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++counter;
        reduce(counter, sumOp<label>());
        Info << "Thinned " << counter << " bnd hair edges" << endl;
        # endif

        //- optimise normals inside the mesh
        optimiseHairNormalsInside();

        //- optimise thickness variation inside the mesh
        optimiseThicknessVariation(INSIDE);

        # ifdef DEBUGLayer
        label intCounter = 0;
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++intCounter;
        Info << "Thinned " << (intCounter - counter)
             << " inner hair edges" << endl;
        # endif
    } while( optimiseLayersAtExittingFaces() && (++nIter < maxNumIterations_) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
