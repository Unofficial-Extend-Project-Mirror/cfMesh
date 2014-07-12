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

#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::calculateHairEdges()
{
    const edgeList& edges = meshSurface_.edges();
    const VRWGraph& bpEdges = meshSurface_.boundaryPointEdges();
    const labelList& bp = meshSurface_.bp();

    //- create information about feature edges and corners at the surface
    meshSurfacePartitioner mPart(meshSurface_);

    //- detect layers in the mesh
    detectBoundaryLayers detectLayers(mPart);

    hairEdges_ = detectLayers.hairEdges();
    hairEdgesAtBndPoint_ = detectLayers.hairEdgesAtBndPoint();

    hairEdgeType_.setSize(hairEdges_.size());
    hairEdgeType_ = INSIDE;

    const labelHashSet& corners = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const labelHashSet& featureEdges = mPart.featureEdges();

    //- classify hair edges
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        //- check if it is a boundary edge
        forAllRow(bpEdges, bpI, peI)
        {
            const label beI = bpEdges(bpI, peI);

            const edge& be = edges[bpEdges(bpI, peI)];

            if( be == e )
            {
                hairEdgeType_[hairEdgeI] = BOUNDARY;

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
    }

    //- calculate which other hair edges influence a hair edges
    //- and store it in a graph
    hairEdgesNearHairEdge_.setSize(hairEdges_.size());

    const labelList& faceOwner = meshSurface_.faceOwners();
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
                        forAllRow(hairEdgesAtBndPoint_, bpJ, hpI)
                        {
                            const label heJ = hairEdgesAtBndPoint_(bpJ, hpI);
                            if( hairEdges_[heJ] == fe )
                                neiHairEdges.append(heJ);
                        }
                    }
                }
            }
        }

        hairEdgesNearHairEdge_.setRow(hairEdgeI, neiHairEdges);
    }
/*
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            Pout << "Hair edges near edge " << hairEdgesNearHairEdge_ << endl;
        }

        returnReduce(1, sumOp<label>());
    }

    //- write hair edges to file
    if( true )
    {
        fileName fName("hairEdges_");
        fName += help::scalarToText(Pstream::myProcNo());
        fName += ".vtk";
        OFstream file(fName);

        //- write the header
        file << "# vtk DataFile Version 3.0\n";
        file << "vtk output\n";
        file << "ASCII\n";
        file << "DATASET POLYDATA\n";

        //- write points
        const pointFieldPMG& points = mesh_.points();
        file << "POINTS " << 2*hairEdges_.size() << " float\n";
        forAll(hairEdges_, eI)
        {
            const point& p0 = points[hairEdges_[eI][0]];
            file << p0.x() << ' ' << p0.y() << ' ' << p0.z() << '\n';
            const point& p1 = points[hairEdges_[eI][1]];
            file << p1.x() << ' ' << p1.y() << ' ' << p1.z() << '\n';
        }

        file << "\nLINES " << hairEdges_.size()
             << ' ' << 3*hairEdges_.size() << nl;
        forAll(hairEdges_, edgeI)
        {
            file << "2 " << (2*edgeI)
                 << token::SPACE << (2*edgeI+1) << nl;
        }
        file << nl;

        if( !file )
            FatalErrorIn
            (
                "void exportFeatureEdges(const triSurf&, const fileName&)"
            ) << "Writting of feature edges failed!" << exit(FatalError);
    }

    //- write hair edges near hair edge
    forAll(hairEdgesNearHairEdge_, eI)
    {
        //- write hair edges to file
        fileName fName("hairEdges_");
        fName += help::scalarToText(Pstream::myProcNo());
        fName += "_point_";
        fName += help::scalarToText(eI);
        fName += ".vtk";
        OFstream file(fName);

        //- write the header
        file << "# vtk DataFile Version 3.0\n";
        file << "vtk output\n";
        file << "ASCII\n";
        file << "DATASET POLYDATA\n";

        DynList<label> nearEdges;
        nearEdges.append(eI);
        forAllRow(hairEdgesNearHairEdge_, eI, i)
            nearEdges.append(hairEdgesNearHairEdge_(eI, i));

        //- write points
        const pointFieldPMG& points = mesh_.points();
        file << "POINTS " << 2*nearEdges.size() << " float\n";
        forAll(nearEdges, i)
        {
            const edge& e = hairEdges_[nearEdges[i]];
            const point& p0 = points[e[0]];
            file << p0.x() << ' ' << p0.y() << ' ' << p0.z() << '\n';
            const point& p1 = points[e[1]];
            file << p1.x() << ' ' << p1.y() << ' ' << p1.z() << '\n';
        }

        file << "\nLINES " << nearEdges.size()
             << ' ' << 3*nearEdges.size() << nl;
        forAll(nearEdges, edgeI)
        {
            file << "2 " << (2*edgeI)
                 << token::SPACE << (2*edgeI+1) << nl;
        }
        file << nl;
    }

    if( true )
    {
        //- write hair edges to file
        fileName fName("hairEdgesAtBnd_");
        fName += help::scalarToText(Pstream::myProcNo());
        fName += ".vtk";
        OFstream file(fName);

        label nBndEdges(0);
        forAll(hairEdgeType_, eI)
            if( hairEdgeType_[eI] & BOUNDARY )
                ++nBndEdges;
        Pout << "Number of bnd edges " << nBndEdges << endl;

        //- write the header
        file << "# vtk DataFile Version 3.0\n";
        file << "vtk output\n";
        file << "ASCII\n";
        file << "DATASET POLYDATA\n";

        //- write points
        const pointFieldPMG& points = mesh_.points();
        file << "POINTS " << 2*nBndEdges << " float\n";
        forAll(hairEdges_, i)
        {
            if( !(hairEdgeType_[i] & BOUNDARY) )
                continue;

            const edge& e = hairEdges_[i];

            const point& p0 = points[e[0]];
            file << p0.x() << ' ' << p0.y() << ' ' << p0.z() << '\n';
            const point& p1 = points[e[1]];
            file << p1.x() << ' ' << p1.y() << ' ' << p1.z() << '\n';
        }

        file << "\nLINES " << nBndEdges
             << ' ' << 3*nBndEdges << nl;
        nBndEdges = 0;
        forAll(hairEdges_, i)
        {
            if( !(hairEdgeType_[i] & BOUNDARY) )
                continue;

            file << "2 " << (2*nBndEdges)
                 << token::SPACE << (2*nBndEdges+1) << nl;
            ++nBndEdges;
        }
        file << nl;
    }
    */
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::optimiseHairNormals(const label nIterations)
{
    pointFieldPMG& points = mesh_.points();

    //- calculate direction of hair vector based on the surface normal
    const labelList& facePatch = meshSurface_.boundaryFacePatches();
    const labelList& bp = meshSurface_.bp();
    const VRWGraph& pointFaces = meshSurface_.pointFaces();
    const vectorField& fNormals = meshSurface_.faceNormals();
    const edgeList& edges = meshSurface_.edges();
    const VRWGraph& bpEdges = meshSurface_.boundaryPointEdges();
    const VRWGraph& edgeFaces = meshSurface_.edgeFaces();

    //- calculate point normals with respect to all patches at a point
    typedef std::map<label, std::pair<point, scalar> > ltvMap;
    typedef std::map<label, ltvMap> lltvMap;
    lltvMap pointPatchNormal;

    forAll(hairEdges_, hairEdgeI)
    {
        const label bpI = bp[hairEdges_[hairEdgeI][0]];

        ltvMap& patchNormal = pointPatchNormal[bpI];

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
        const Map<label>& globalToLocal =
            meshSurface_.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = meshSurface_.bpNeiProcs();
        const VRWGraph& bpAtProcs = meshSurface_.bpAtProcs();

        std::map<label, LongList<refLabelledPointScalar> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const ltvMap& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(ltvMap, patchNormal, pIt)
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

            ltvMap& patchNormal = pointPatchNormal[bpI];

            const labelledPointScalar& lps = rlps.lps();
            patchNormal[lps.pointLabel()].first += lps.coordinates();
            patchNormal[lps.pointLabel()].second += lps.scalarValue();
        }
    }

    forAllIter(lltvMap, pointPatchNormal, it)
    {
        ltvMap& patchNormal = it->second;

        forAllIter(ltvMap, patchNormal, pIt)
            pIt->second.first /= pIt->second.second;
    }

    //- calculate hair vectors
    //- they point in the normal direction to the surface
    vectorField hairVecs(hairEdges_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const direction hairType = hairEdgeType_[hairEdgeI];
        const label bpI = bp[hairEdges_[hairEdgeI][0]];

        vector& hv = hairVecs[hairEdgeI];

        const ltvMap& patchNormals = pointPatchNormal[bpI];

        if( hairType & INSIDE )
        {
            hv = vector::zero;
            label counter(0);
            forAllConstIter(ltvMap, patchNormals, pIt)
            {
                hv -= pIt->second.first;
                ++counter;
            }

            if( counter == 0 )
            {
                FatalErrorIn
                (
                    "void boundaryLayerOptimisation::optimiseHairNormals()"
                ) << "No valid patches for boundary point "
                  << bpI << abort(FatalError);
            }

            hv /= counter;
            hv /= (mag(hv) + VSMALL);
        }
        else if( hairType & BOUNDARY )
        {
            if( hairType & FEATUREEDGE )
            {
                hv = hairEdges_[hairEdgeI].vec(points);
                hv /= (mag(hv)+VSMALL);
            }
            else
            {
                //- find the patch to which the hair edge belongs to
                label edgeInPatch(-1);
                forAllRow(bpEdges, bpI, bpeI)
                {
                    const label beI = bpEdges(bpI, bpeI);

                    if( edges[beI] == hairEdges_[hairEdgeI] )
                    {
                        edgeInPatch = facePatch[edgeFaces(beI, 0)];
                        break;
                    }
                }

                if( edgeInPatch < 0 )
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormals()"
                    ) << "Cannot find the patch" << abort(FatalError);
                }

                //- collect data based on patches at a point
                vector patchNormal;
                hv = vector::zero;
                label counter(0);
                forAllConstIter(ltvMap, patchNormals, pIt)
                {
                    if( pIt->first == edgeInPatch )
                    {
                        patchNormal = pIt->second.first;
                    }
                    else
                    {
                        hv -= pIt->second.first;
                        ++counter;
                    }
                }

                if( counter == 0 )
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormals()"
                    ) << "No valid patches for boundary point "
                      << bpI << abort(FatalError);
                }

                hv /= counter;

                const point& p = points[hairEdges_[hairEdgeI][0]];
                const plane pl(p, patchNormal);
                const point np = pl.nearestPoint(p + hv);
                hv = np - p;

                hv /= (mag(hv) + VSMALL);
            }
        }
    }

    //- smooth the variation of normals to reduce the twisting of faces
    label nIter(0);
    do
    {
        vectorField newNormals(hairEdges_.size());
        labelList nNeighbours(hairEdges_.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            const edge& e = hairEdges_[hairEdgeI];

            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;
            label& nNei = nNeighbours[hairEdgeI];
            nNei = 0;

            const direction eType = hairEdgeType_[hairEdgeI];

            if( eType & INSIDE )
            {
                if( eType & ATEDGE )
                {
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        const edge& ne = hairEdges_[hairEdgeJ];

                        if( hairEdgeType_[hairEdgeJ] & (ATEDGE | ATCORNER) )
                        {
                            const vector distVec
                            (
                                points[e[0]] -
                                points[ne[0]]
                            );

                            vector n = hairVecs[hairEdgeJ] ^ distVec;
                            n /= (mag(n) + VSMALL);

                            //- project the hair vector into the plane
                            //- detemined by the starting poins of the hair
                            //- edge and the normal n
                            const vector hairProj
                            (
                                hairVecs[hairEdgeI] -
                                (hairVecs[hairEdgeI] & n) * n
                            );

                            newNormal += hairProj;
                            ++nNei;
                        }
                    }
                }
                else if( eType == INSIDE )
                {
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        const edge& ne = hairEdges_[hairEdgeJ];

                        const vector distVec
                        (
                            points[e[0]] -
                            points[ne[0]]
                        );

                        vector n = hairVecs[hairEdgeJ] ^ distVec;
                        n /= (mag(n) + VSMALL);

                        //- project the hair vector into the plane detemined
                        //- by the starting poins of the hair edge and
                        //- the normal n
                        const vector hairProj
                        (
                            hairVecs[hairEdgeI] -
                            (hairVecs[hairEdgeI] & n) * n
                        );

                        newNormal += hairProj;
                        ++nNei;
                    }
                }
                else
                {
                    newNormal += hairVecs[hairEdgeI];
                    ++nNei;
                }
            }
            else if( eType & BOUNDARY )
            {
                if( eType & FEATUREEDGE )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal += hairVecs[hairEdgeI];
                    ++nNei;
                }
                else
                {
                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        if( !(hairEdgeType_[hairEdgeJ] & BOUNDARY) )
                            continue;

                        const edge& ne = hairEdges_[hairEdgeJ];

                        const vector distVec
                        (
                            points[e[0]] -
                            points[ne[0]]
                        );

                        vector n = hairVecs[hairEdgeJ] ^ distVec;
                        n /= (mag(n) + VSMALL);

                        //- project the hair vector into the plane detemined
                        //- by the starting points of the hair edge and
                        //- the normal n
                        const vector hairProj
                        (
                            hairVecs[hairEdgeI] -
                            (hairVecs[hairEdgeI] & n) * n
                        );

                        newNormal += hairProj;
                        ++nNei;
                    }
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                meshSurface_.globalToLocalBndPointAddressing();

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

            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                    const edge& he = hairEdges_[hairEdgeI];

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            refLabelledPoint rlp
                            (
                                globalEdgeLabel[edgeI],
                                labelledPoint
                                (
                                    nNeighbours[hairEdgeI],
                                    hairVecs[hairEdgeI]
                                )
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<refLabelledPoint>& dts =
                                    exchangeData[neiProc];

                                dts.append(rlp);
                            }
                        }
                    }
                }
            }

            LongList<refLabelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const refLabelledPoint& rlp = receivedData[i];
                const label edgeI = globalToLocalEdge[rlp.objectLabel()];
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
                            hairVecs[hairEdgeI] += rlp.lPoint().coordinates();
                            nNeighbours[hairEdgeI] +=
                                rlp.lPoint().pointLabel();

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
            if( nNeighbours[hairEdgeI] == 0 )
            {
                bool isBnd(false);
                if( hairEdgeType_[hairEdgeI] & BOUNDARY )
                    isBnd = true;
                FatalError << isBnd << " Edge " << hairEdgeI << " at proc "
                    << Pstream::myProcNo() << " failed" << abort(FatalError);
            }
            newNormals[hairEdgeI] /= nNeighbours[hairEdgeI];
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);
    } while( nIter++ < nIterations );

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const edge& he = hairEdges_[hairEdgeI];

        const vector hv = (hairVecs[hairEdgeI] / mag(hairVecs[hairEdgeI]));

        points[he.end()] = points[he.start()] + hv * he.mag(points);
    }
}

void boundaryLayerOptimisation::optimiseThicknessVariation
(
    const scalar tangentTol
)
{
    pointFieldPMG& points = mesh_.points();

    //- reduce thickness of the layer
    //- such that the variation of layer thickness
    //- It is an iterative process where the layer is thinned in the regions
    //- where the tangent is greater than the tolerance value.
    vectorField hairDirections(hairEdges_.size());
    scalarField hairLength(hairEdges_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        vector n = hairEdges_[hairEdgeI].vec(points);

        hairLength[hairEdgeI] = Foam::max(VSMALL, Foam::mag(n));
        hairDirections[hairEdgeI] = n / hairLength[hairEdgeI];
    }

    bool changed;
    label nIter(0);
    do
    {
        changed = false;

        scalarField newLength(hairLength);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            const scalar magN = hairLength[hairEdgeI];

            if( magN < VSMALL )
                continue;

            forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
            {
                const label hairEdgeJ =
                    hairEdgesNearHairEdge_(hairEdgeI, nheI);

                const scalar magNeiNormal = hairLength[hairEdgeJ];

                const vector distVec
                (
                    points[hairEdges_[hairEdgeI][0]] -
                    points[hairEdges_[hairEdgeJ][0]]
                );

                const scalar magDistVec = mag(distVec);

                const scalar tanVec
                (
                    (magN - magNeiNormal) /
                    (magDistVec + VSMALL)
                );

                if( tanVec > tangentTol )
                {
                    //- make the hair edge shorter
                    newLength[hairEdgeI] =
                        hairLength[hairEdgeJ] + magDistVec * tangentTol;

                    changed = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                meshSurface_.globalToLocalBndPointAddressing();
            const labelList& bp = meshSurface_.bp();

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
            if( hairEdgeType_[hairEdgeI] & BOUNDARY )
            {
                const edge& he = hairEdges_[hairEdgeI];
                const vector& hv = hairDirections[hairEdgeI];
                const point& s = points[he.start()];

                points[he.end()] = s + newLength[hairEdgeI] * hv;
            }
        }

        //- move intenal vertices to the new positions
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            if( hairEdgeType_[hairEdgeI] & INSIDE )
            {
                const edge& he = hairEdges_[hairEdgeI];
                const vector& hv = hairDirections[hairEdgeI];
                const point& s = points[he.start()];

                points[he.end()] = s + newLength[hairEdgeI] * hv;
            }
        }

        //- update hairLength
        hairLength.transfer(newLength);

    } while( changed && (++nIter < 20) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
