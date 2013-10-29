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
#include "demandDrivenData.H"
#include "FixedList.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::refineFace
(
    const face& f,
    const FixedList<label, 2>& nLayersInDirection,
    DynList<DynList<label, 4> >& newFaces
)
{
    //- this face must be a quad
    if( f.size() != 4 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineFace(const face&,"
            " const FixedList<label, 2>&, DynList<DynList<label, 4> >&)"
        ) << "Face " << f << " is not a quad" << endl;
        return;
    }

    //- direction 0 represents edges 0 and 2
    //- direction 1 represents edges 1 and 3
    if( (nLayersInDirection[0] <= 1) && (nLayersInDirection[1] <= 1) )
    {
        newFaces.setSize(1);
        newFaces[0] = f;
        return;
    }

    //- check which face edge is a direction 0 and which one is a direction 1
    label dir0(-1), dir1(-1);
    labelPair dir0Edges(-1, -1), dir1Edges(-1, -1);
    forAll(f, eI)
    {
        const edge e = f.faceEdge(eI);

        label ses(-1), see(-1);
        bool start(false), end(false);
        forAllRow(splitEdgesAtPoint_, e.start(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_(e.start(), i)];

            if( (se.start() == e.start()) && (se.end() == f.prevLabel(eI)) )
            {
                ses = splitEdgesAtPoint_(e.start(), i);
                start = true;
                break;
            }
        }

        forAllRow(splitEdgesAtPoint_, e.end(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_(e.end(), i)];

            if( (se.start() == e.end()) && (se.end() == f[(eI+2)%4]) )
            {
                see = splitEdgesAtPoint_(e.end(), i);
                end = true;
                break;
            }
        }

        if( start && end )
        {
            if( dir0 == -1 )
            {
                dir0 = eI;
                dir0Edges = labelPair(ses, see);
            }
            else if( dir1 == -1 )
            {
                dir1 = eI;
                dir1Edges = labelPair(ses, see);
            }
            else
            {
                FatalErrorIn
                (
                    "void refineBoundaryLayers::refineFace(const face&,"
                    " const FixedList<label, 2>&, DynList<DynList<label, 4> >&)"
                ) << "More than two split directions for a face"
                  << abort(FatalError);
            }
        }
    }

    /*
    Info << "Refining face " << f << endl;
    Info << "Splits in direction " << nLayersInDirection << endl;
    Info << "Here " << endl;
    Info << "Dir0 " << dir0 << endl;
    Info << "dir0Edges " << dir0Edges << endl;
    Info << "Dir1 " << dir1 << endl;
    Info << "dir1Edges " << dir1Edges << endl;
    */

    //- in case of only one refinement direction, it must direction 0
    if( (dir1 != -1) && (dir0 == -1) )
    {
        dir0 = dir1;
        dir0Edges = dir1Edges;
        dir1 = -1;
    }
    else if( (dir0 != -1) && (dir1 != -1) && (dir1 != f.rcIndex(dir0)) )
    {
        //- alternate value to preserve correct face orientation
        const label add = dir0;
        dir0 = dir1;
        dir1 = add;

        const labelPair lpAdd = dir0Edges;
        dir0Edges = dir1Edges;
        dir1Edges = lpAdd;
    }

    //- permutate the number of refinements in each direction
    const label nLayersDir0 = dir0>=0?nLayersInDirection[dir0%2]:1;
    const label nLayersDir1 = dir1>=0?nLayersInDirection[dir1%2]:1;

    /*
    Info << "Face has points " << f << endl;
    Info << "dirEdges0 " << dir0Edges << endl;
    Info << "dir1Edges " << dir1Edges << endl;
    if( dir0 >= 0 )
    {
        Info << "Points on edge " << dir0Edges.first()
             << " are " << newVerticesForSplitEdge_[dir0Edges.first()] << endl;
        Info << "Points on edge " << dir0Edges.second()
             << " are " << newVerticesForSplitEdge_[dir0Edges.second()] << endl;
    }
    if( dir1 >= 0 )
    {
        Info << "Points on edge " << dir1Edges.first()
             << " are " << newVerticesForSplitEdge_[dir1Edges.first()] << endl;
        Info << "Points on edge " << dir1Edges.second()
             << " are " << newVerticesForSplitEdge_[dir1Edges.second()] << endl;
    }
    Info << "nLayersDir0 " << nLayersDir0 << endl;
    Info << "nLayersDir1 " << nLayersDir1 << endl;
    */

    //- map the face onto a matrix for easier orientation
    DynList<DynList<label> > facePoints;
    facePoints.setSize(nLayersDir0+1);
    forAll(facePoints, i)
    {
        facePoints[i].setSize(nLayersDir1+1);
        facePoints[i] = -1;
    }

    //- add points in the matrix
    for(label i=0;i<nLayersDir0;++i)
    {
        facePoints[i][0] = newVerticesForSplitEdge_(dir0Edges.second(), i);
        facePoints[i][nLayersDir1] =
            newVerticesForSplitEdge_(dir0Edges.first(), i);
    }
    facePoints[nLayersDir0][0] = splitEdges_[dir0Edges.second()].end();
    facePoints[nLayersDir0][nLayersDir1] = splitEdges_[dir0Edges.first()].end();

    for(label i=1;i<nLayersDir1;++i)
    {
        facePoints[0][i] = newVerticesForSplitEdge_(dir1Edges.first(), i);
        facePoints[nLayersDir0][i] =
            newVerticesForSplitEdge_(dir1Edges.second(), i);
    }

    //- create missing vertices if there are any
    const pointFieldPMG& points = mesh_.points();

    forAll(facePoints, i)
    {
        forAll(facePoints[i], j)
        {
            if( facePoints[i][j] < 0 )
            {
                const scalar u
                (
                    mag(points[facePoints[i][0]] - points[facePoints[0][0]]) /
                    splitEdges_[dir0Edges.second()].mag(points)
                );

                const scalar v
                (
                    mag(points[facePoints[0][i]] - points[facePoints[0][0]]) /
                    splitEdges_[dir1Edges.first()].mag(points)
                );

                //- calculate the coordinates of the missing point via
                //- transfinite interpolation
                const point newP
                (
                    (1.0 - v) * points[facePoints[i][0]] +
                    v * points[facePoints[i][nLayersDir1]] +
                    (1.0 - u) * points[facePoints[0][j]] +
                    u * points[facePoints[nLayersDir0][j]] -
                    (1.0 - u) * (1.0 - v) * points[facePoints[0][0]] -
                    u * v * points[facePoints[nLayersDir0][0]] -
                    u * (1.0 - v) * points[facePoints[nLayersDir0][nLayersDir1]] -
                    (1.0 - u) * v * points[facePoints[0][nLayersDir1]]
                );

                Info << "Point coordinate " << newP << endl;

                //- add the vertex to the mesh
                facePoints[i][j] = points.size();
                mesh_.appendVertex(newP);
            }
        }
    }

    //Info << "Face points after creating vertices " << facePoints << endl;

    //- Finally, create the faces
    for(label i=0;i<nLayersDir0;++i)
    {
        for(label j=0;j<nLayersDir1;++j)
        {
            if( !((i == (nLayersDir0 - 1)) && (j == (nLayersDir1 - 1))) )
            {
                //- create quad face
                DynList<label, 4> f;
                f.setSize(4);

                f[0] = facePoints[i][j];
                f[1] = facePoints[i+1][j];
                f[2] = facePoints[i+1][j+1];
                f[3] = facePoints[i][j+1];

                newFaces.append(f);
            }
        }
    }

    //- create the last face which may not be a quad
    DynList<label, 4> newF;

    if( dir0 != -1 && dir1 == -1 )
    {
        //- face is split in one direction, only
        label eLabel = dir0Edges.second();
        label size = newVerticesForSplitEdge_.sizeOfRow(eLabel);

        for(label i=nLayersDir0-1;i<size;++i)
            newF.append(newVerticesForSplitEdge_(eLabel, i));

        eLabel = dir0Edges.first();
        size = newVerticesForSplitEdge_.sizeOfRow(eLabel);
        for(label i=size;i>=nLayersDir0;--i)
            newF.append(newVerticesForSplitEdge_(eLabel, i-1));
    }
    else if( dir0 != -1 && dir1 != -1 )
    {
        //- face is split in both directions
        newF.append(facePoints[nLayersDir0-1][nLayersDir1-1]);

        //- add additional points on edge
        label eLabel = dir1Edges.second();
        label size = newVerticesForSplitEdge_.sizeOfRow(eLabel) - 1;

        for(label i=nLayersDir1-1;i<size;++i)
            newF.append(newVerticesForSplitEdge_(eLabel, i));

        //- add other corner
        newF.append(facePoints[nLayersDir0][nLayersDir1]);

        //- add additional points on edge
        eLabel = dir0Edges.first();
        size = newVerticesForSplitEdge_.sizeOfRow(eLabel) - 1;
        for(label i=size;i>=nLayersDir0;--i)
            newF.append(newVerticesForSplitEdge_(eLabel, i-1));
    }

    newFaces.append(newF);

    //Info << "Input face " << f << endl;
    //Info << "Decomposed faces are " << newFaces << endl;
    //if( (nLayersInDirection[0] > 1) && (nLayersInDirection[1] > 1) )
    //::exit(1);
}

void refineBoundaryLayers::generateNewFaces()
{
    //- generate new boundary and inter-processor faces
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatches = mse.boundaryFacePatches();
    const edgeList& edges = mse.edges();
    const labelList& bp = mse.bp();
    const VRWGraph& bfEdges = mse.faceEdges();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const VRWGraph& beFaces = mse.edgeFaces();

    //- mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    //- container for faces
    facesFromFace_.setSize(faces.size());
    newFaces_.clear();

    //- split internal faces
    for(label faceI=0;faceI<nInternalFaces;++faceI)
    {
        const face& f = faces[faceI];

        //- only quad faces can be split
        if( f.size() != 4 )
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(f);
            continue;
        }

        //- check if there exist an edge of the face at the boundary
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(f, eI)
        {
            const edge fe = f.faceEdge(eI);

            const label bps = bp[fe.start()];

            if( bps < 0 )
                continue;

            forAllRow(bpEdges, bps, bpsI)
            {
                const label beI = bpEdges(bps, bpsI);

                if( edges[beI] == fe )
                {
                    //- this edge is attached to the boundary
                    //- get the number of layers for neighbouring cells
                    const label nSplits0 = nLayersAtBndFace_[beFaces(beI, 0)];
                    const label nSplits1 = nLayersAtBndFace_[beFaces(beI, 1)];

                    //- set the number of layers for the given direction
                    const label dir = eI % 2;
                    nRefinementInDirection[dir] =
                        Foam::max
                        (
                            nRefinementInDirection[dir],
                            Foam::max(nSplits0, nSplits1)
                        );
                }
            }
        }

        //- refine the face
        DynList<DynList<label, 4> > newFacesForFace;

        refineFace(f, nRefinementInDirection, newFacesForFace);

        //- store decomposed faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(newFacesForFace[fI]);
        }
    }

    //- refine boundary faces where needed
    //- it is required in locations where two or three layers intersect
    const label startingBoundaryFace = mesh_.boundaries()[0].patchStart();
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const label faceI = startingBoundaryFace + bfI;

        //- only quad faces can be split
        if( bf.size() != 4 )
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(bf);
            continue;
        }

        //- check whether this face shall be refined and in which directions
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(bf, eI)
        {
            const label beI = bfEdges(bfI, eI);

            //- get the neighbour face over the edge
            label neiFace = beFaces(beI, 0);

            if( beFaces.sizeOfRow(beI) != 2 )
                continue;

            if( neiFace == bfI )
                neiFace = beFaces(beI, 1);

            //- faces cannot be in the same layer
            if(
                layerAtPatch_[facePatches[neiFace]] ==
                layerAtPatch_[facePatches[bfI]]
            )
                continue;

            //- set the refinement direction for this face
            nRefinementInDirection[eI%2] = nLayersAtBndFace_[neiFace];
        }

        //- refine the face
        DynList<DynList<label, 4> > newFacesForFace;

        refineFace(bf, nRefinementInDirection, newFacesForFace);

        //- store the refined faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(newFacesForFace[fI]);
        }
    }

    if( Pstream::parRun() )
    {
        //- refine faces at interprocessor boundaries
        const PtrList<writeProcessorPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- exchange information about the number of splits
        //- to other processors
        std::map<label, DynList<labelPair, 2> > localSplits;
        forAll(procBoundaries, patchI)
        {
            labelListPMG sendData;

            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;
                const face& f = faces[faceI];

                forAll(f, eI)
                {
                    const edge fe = f.faceEdge(eI);

                    const label bps = bp[fe.start()];

                    if( bps < 0 )
                        continue;

                    forAllRow(bpEdges, bps, bpeI)
                    {
                        const label beI = bpEdges(bps, bpeI);

                        if( edges[beI] == fe )
                        {
                            //- this edge is attached to the boundary
                            //- get the number of layers for neighbouring cell
                            const label nSplits0 =
                                nLayersAtBndFace_[beFaces(beI, 0)];

                            //- add the data to the list for sending
                            const label dir = !(eI % 2);

                            //- add face label, direction
                            //- and the number of splits
                            sendData.append(fI);
                            sendData.append(dir);
                            sendData.append(nSplits0);
                            localSplits[faceI].append(labelPair(dir, nSplits0));
                        }
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        //- receive data from other procesors
        forAll(procBoundaries, patchI)
        {
            //- get the data sent from the neighbour processor
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const label fI = receivedData[counter++];
                const label dir = receivedData[counter++];
                const label nSplits = receivedData[counter++];

                DynList<labelPair, 2>& currentSplits = localSplits[start+fI];
                forAll(currentSplits, i)
                {
                    if( currentSplits[i].first() == dir )
                        currentSplits[i].second() =
                            Foam::max(currentSplits[i].second(), nSplits);
                }
            }
        }

        //- perform splitting
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;

                std::map<label, DynList<labelPair, 2> >::const_iterator it =
                    localSplits.find(faceI);

                if( it == localSplits.end() )
                {
                    //- this face is not split
                    facesFromFace_.append(faceI, newFaces_.size());
                    newFaces_.appendList(faces[faceI]);
                    continue;
                }

                //- split the face and add the faces to the list
                DynList<DynList<label, 4> > facesFromFace;
                if( procBoundaries[patchI].owner() )
                {
                    //- this processor own this patch
                    FixedList<label, 2> nLayersInDirection;
                    const DynList<labelPair, 2>& dirSplits = it->second;
                    forAll(dirSplits, i)
                        nLayersInDirection[dirSplits[i].first()] =
                            dirSplits[i].second();

                    refineFace(faces[faceI], nLayersInDirection, facesFromFace);

                    //- add faces
                    forAll(facesFromFace, i)
                    {
                        facesFromFace_.append(faceI, newFaces_.size());
                        newFaces_.appendList(facesFromFace[i]);
                    }
                }
                else
                {
                    //- reverse the face before splitting
                    FixedList<label, 2> nLayersInDirection;
                    const DynList<labelPair, 2>& dirSplits = it->second;
                    forAll(dirSplits, i)
                        nLayersInDirection[(dirSplits[i].first()+1)%2] =
                            dirSplits[i].second();

                    const face rFace = faces[faceI].reverseFace();
                    refineFace(rFace, nLayersInDirection, facesFromFace);

                    forAll(facesFromFace, i)
                    {
                        const DynList<label, 4>& df = facesFromFace[i];
                        DynList<label, 4> rFace = help::reverseFace(df);

                        facesFromFace_.append(faceI, newFaces_.size());
                        newFaces_.appendList(rFace);
                    }
                }
            }
        }
    }

    Info << "facesFromFace_ " << facesFromFace_ << endl;
    Info << "newFaces_ " << newFaces_ << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

