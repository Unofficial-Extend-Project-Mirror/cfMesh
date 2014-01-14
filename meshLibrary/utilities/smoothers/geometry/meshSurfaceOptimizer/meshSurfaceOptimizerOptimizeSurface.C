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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshOctree.H"
#include "triangle.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "labelledPoint.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label meshSurfaceOptimizer::findInvertedVertices
(
    boolList& smoothVertex
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();

    if( smoothVertex.size() != bPoints.size() )
    {
        smoothVertex.setSize(bPoints.size());
        smoothVertex = true;
    }

    label nInvertedTria(0);

    //- check the vertices at the surface
    //- mark the ones where the mesh is tangled
    meshSurfaceCheckInvertedVertices vrtCheck(surfaceEngine_, &smoothVertex);
    const labelHashSet& inverted = vrtCheck.invertedVertices();

    smoothVertex = false;
    forAll(bPoints, bpI)
    {
        if( inverted.found(bPoints[bpI]) )
        {
            ++nInvertedTria;
            smoothVertex[bpI] = true;
        }
    }

    if( Pstream::parRun() )
        reduce(nInvertedTria, sumOp<label>());
    Info << "Number of inverted boundary faces is " << nInvertedTria << endl;

    if( nInvertedTria == 0 )
        return 0;

    //- add additional layers around inverted points
    for(label i=0;i<2;++i)
    {
        boolList originallySelected = smoothVertex;
        forAll(smoothVertex, bpI)
            if( originallySelected[bpI] )
                forAllRow(pPoints, bpI, ppI)
                    smoothVertex[pPoints(bpI, ppI)] = true;

        if( Pstream::parRun() )
        {
            //- exchange global labels of inverted points
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

            std::map<label, labelListPMG> shareData;
            forAll(neiProcs, procI)
                shareData.insert
                (
                    std::make_pair(neiProcs[procI], labelListPMG())
                );

            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                if( !smoothVertex[bpI] )
                    continue;

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    shareData[neiProc].append(globalPointLabel[bpI]);
                }
            }

            //- exchange data with other processors
            labelListPMG receivedData;
            help::exchangeMap(shareData, receivedData);

            forAll(receivedData, j)
            {
                const label bpI = globalToLocal[receivedData[j]];

                smoothVertex[bpI] = true;
            }
        }
    }

    return nInvertedTria;
}

bool meshSurfaceOptimizer::preOptimizeSurface()
{
    Info << "Optimizing positions of surface nodes" << endl;

    bool changed(false);

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();
    this->triangles();
    this->pointTriangles();

    if( Pstream::parRun() )
    {
        surfaceEngine_.bpAtProcs();
        surfaceEngine_.globalToLocalBndPointAddressing();
        surfaceEngine_.globalBoundaryPointLabel();
        surfaceEngine_.bpNeiProcs();
    }

    boolList smoothVertex;

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    meshSurfaceMapper mapper(surfaceEngine_, meshOctree_);

    bool remapVertex(true);
    label nInvertedTria;
    label nGlobalIter(0);

    labelListPMG procBndNodes, movedPoints;

    do
    {
        label nIter(0);

        do
        {
            nInvertedTria = findInvertedVertices(smoothVertex);

            if( nInvertedTria == 0 ) break;

            changed = true;

            procBndNodes.clear();
            movedPoints.clear();
            forAll(bPoints, bpI)
            {
                if( smoothVertex[bpI] && (vertexType_[bpI] & PARTITION) )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                    {
                        procBndNodes.append(bpI);
                        continue;
                    }
                }
            }

            //- use laplacian smoothing
            # ifdef USE_OMP
            # pragma omp parallel
            # endif
            {
                LongList<labelledPoint> newPos;

                # ifdef USE_OMP
                # pragma omp for schedule(guided)
                # endif
                forAll(movedPoints, i)
                {
                    const label bpI = movedPoints[i];

                    if( vertexType_[bpI] & PROCBND )
                        continue;

                    newPos.append
                    (
                        labelledPoint(bpI, newPositionLaplacianFC(bpI))
                    );
                }

                forAll(newPos, i)
                    surfaceModifier.moveBoundaryVertexNoUpdate
                    (
                        newPos[i].pointLabel(),
                        newPos[i].coordinates()
                    );
            }

            if( Pstream::parRun() )
                nodeDisplacementLaplacianFCParallel(procBndNodes, true);

            surfaceModifier.updateGeometry(movedPoints);

            //- use surface optimizer
            # ifdef USE_OMP
            # pragma omp parallel
            # endif
            {
                LongList<labelledPoint> newPos;

                # ifdef USE_OMP
                # pragma omp for schedule(guided)
                # endif
                forAll(movedPoints, i)
                {
                    const label bpI = movedPoints[i];

                    if( vertexType_[bpI] & PROCBND )
                        continue;

                    newPos.append
                    (
                        labelledPoint(bpI, newPositionSurfaceOptimizer(bpI))
                    );
                }

                forAll(newPos, i)
                    surfaceModifier.moveBoundaryVertexNoUpdate
                    (
                        newPos[i].pointLabel(),
                        newPos[i].coordinates()
                    );
            }

            if( Pstream::parRun() )
                nodeDisplacementSurfaceOptimizerParallel(procBndNodes);

            surfaceModifier.updateGeometry(movedPoints);

            if( remapVertex )
                mapper.mapVerticesOntoSurface(movedPoints);

        } while( nInvertedTria && (++nIter < 20) );

        if( nInvertedTria )
        {
            Info << "Smoothing remaining inverted vertices " << endl;

            movedPoints.clear();
            procBndNodes.clear();
            forAll(smoothVertex, bpI)
                if( smoothVertex[bpI] )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                    {
                        procBndNodes.append(bpI);
                        continue;
                    }

                    nodeDisplacementLaplacianFC(bpI, false);
                }

            if( Pstream::parRun() )
            {
                nodeDisplacementLaplacianFCParallel(procBndNodes, false);
            }

            if( remapVertex )
                mapper.mapVerticesOntoSurface(movedPoints);

            if( nGlobalIter > 3 )
                remapVertex = false;
        }

    } while( nInvertedTria && (++nGlobalIter < 10) );

    Info << "Finished optimizing positions of surface nodes" << endl;

    return changed;
}

void meshSurfaceOptimizer::optimizeSurface(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();
    this->triangles();

    labelListPMG procBndNodes, edgePoints;
    forAll(bPoints, bpI)
    {
        if( vertexType_[bpI] & EDGE )
        {
            edgePoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndNodes.append(bpI);
        }
    }

    meshSurfaceMapper mapper(surfaceEngine_, meshOctree_);

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);
        # ifdef USE_OMP
        # pragma omp parallel if( edgePoints.size() > 1000 )
        # endif
        {
            LongList<labelledPoint> newPos;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 20)
            # endif
            forAll(edgePoints, epI)
            {
                const label bpI = edgePoints[epI];

                if( vertexType_[bpI] & PROCBND )
                    continue;

                const point newP = newEdgePositionLaplacian(bpI);
                newPos.append(labelledPoint(bpI, newP));
            }

            forAll(newPos, i)
                bMod.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
        }

        bMod.updateGeometry(edgePoints);

        if( Pstream::parRun() )
        {
            edgeNodeDisplacementParallel(procBndNodes);
        }

        mapper.mapEdgeNodes(edgePoints);
    }
    Info << endl;

    //- optimize nodes of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        procBndNodes.clear();

        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);
        # ifdef USE_OMP
        # pragma omp parallel if( vertexType_.size() > 100 )
        # endif
        {
            LongList<labelledPoint> newPos;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(bPoints, bpI)
                if( vertexType_[bpI] & PARTITION )
                {
                    if( vertexType_[bpI] & PROCBND )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            procBndNodes.append(bpI);
                        }

                        continue;
                    }

                    const point newP = newPositionLaplacianFC(bpI);
                    newPos.append(labelledPoint(bpI, newP));
                }

            forAll(newPos, i)
                bMod.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
        }

        bMod.updateGeometry();

        if( Pstream::parRun() )
        {
            nodeDisplacementLaplacianFCParallel(procBndNodes,true);
        }
    }

    Info << endl;
}

void meshSurfaceOptimizer::optimizeSurface2D(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine
    (
        const_cast<polyMeshGen&>(surfaceEngine_.mesh())
    );
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    this->triangles();

    labelListPMG procBndNodes, edgePoints, activeEdges, updatePoints;
    forAll(edges, beI)
    {
        const edge& e = edges[beI];

        if( zMinPoint[e.start()] ^ zMinPoint[e.end()] )
        {
            label bpI = bp[e.start()];
            if( !zMinPoint[e.start()] )
                bpI = bp[e.end()];

            if( vertexType_[bpI] & EDGE )
            {
                activeEdges.append(beI);

                updatePoints.append(bp[e.start()]);
                updatePoints.append(bp[e.end()]);

                edgePoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndNodes.append(bpI);
            }
        }
    }

    meshSurfaceMapper2D mapper(surfaceEngine_, meshOctree_);

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);
        # ifdef USE_OMP
        # pragma omp parallel if( edgePoints.size() > 1000 )
        # endif
        {
            LongList<labelledPoint> newPos;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 20)
            # endif
            forAll(edgePoints, epI)
            {
                const label bpI = edgePoints[epI];

                if( vertexType_[bpI] & PROCBND )
                    continue;

                const point newP = newEdgePositionLaplacian(bpI);
                newPos.append(labelledPoint(bpI, newP));
            }

            forAll(newPos, i)
                bMod.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
        }

        if( Pstream::parRun() )
        {
            edgeNodeDisplacementParallel(procBndNodes);
        }

        //- move points with maximum z coordinate
        mesh2DEngine.correctPoints();

        //- update normal, centres, etc, after the surface has been modified
        bMod.updateGeometry(updatePoints);

        //- map boundary edges to the surface
        mapper.mapVerticesOntoSurfacePatches(activeEdges);
    }
    Info << endl;

    //- optimize nodes of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        procBndNodes.clear();

        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);
        # ifdef USE_OMP
        # pragma omp parallel if( vertexType_.size() > 100 )
        # endif
        {
            LongList<labelledPoint> newPos;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(bPoints, bpI)
                if( zMinPoint[bPoints[bpI]] && (vertexType_[bpI] & PARTITION) )
                {
                    if( vertexType_[bpI] & PROCBND )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            procBndNodes.append(bpI);
                        }

                        continue;
                    }

                    const point newP = newPositionLaplacianFC(bpI);

                    newPos.append(labelledPoint(bpI, newP));
                }

            forAll(newPos, i)
                bMod.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
        }

        if( Pstream::parRun() )
        {
            nodeDisplacementLaplacianFCParallel(procBndNodes, true);
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry();
    }

    Info << endl;
}

void meshSurfaceOptimizer::untangleSurface2D()
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    this->triangles();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces = findBadFaces(badFaces, changedFace);

        Info << "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if( nBadFaces == 0 )
            break;

        //- update active points and faces affected by the movement
        //- of active points
        activeBoundaryPoint = false;
        changedFace = false;
        forAllConstIter(labelHashSet, badFaces, it)
        {
            const face& f = faces[it.key()];

            forAll(f, pI)
            {
                if( zMinPoint[f[pI]] )
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces(f[pI], pfI)] = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
            const VRWGraph& bpNeiProcs = surfaceEngine_.bpAtProcs();

            std::map<label, labelListPMG> exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- collect active points at inter-processor boundaries
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( activeBoundaryPoint[bpI] )
                {
                    forAllRow(bpNeiProcs, bpI, i)
                    {
                        const label neiProc = bpNeiProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append(it.key());
                    }
                }
            }

            //- exchange active points among the processors
            labelListPMG receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- ensure that all processors have the same nodes active
            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i]];

                //- activate this boundary point
                activeBoundaryPoint[bpI] = true;

                //- set the changeFaces for the faces attached to this point
                forAllRow(pointFaces, bPoints[bpI], pfI)
                    changedFace[pointFaces(bPoints[bpI], pfI)] = true;
            }
        }

        //- apply smoothing to the activated points
        meshSurfaceEngineModifier bMod(surfaceEngine_);

        for(label i=0;i<5;++i)
        {
            labelListPMG procBndNodes, procEdgeNodes;

            # ifdef USE_OMP
            # pragma omp parallel if( vertexType_.size() > 100 )
            # endif
            {
                LongList<labelledPoint> newPos;

                # ifdef USE_OMP
                # pragma omp for schedule(dynamic, 10)
                # endif
                forAll(bPoints, bpI)
                {
                    if( !activeBoundaryPoint[bpI] )
                        continue;

                    if( vertexType_[bpI] & PARTITION )
                    {
                        if( vertexType_[bpI] & PROCBND )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical
                            # endif
                            {
                                procBndNodes.append(bpI);
                            }

                            continue;
                        }

                        const point newP = newPositionSurfaceOptimizer(bpI);
                        newPos.append(labelledPoint(bpI, newP));
                    }
                    else if( (vertexType_[bpI] & EDGE) && (iterationI > 4) )
                    {
                        if( vertexType_[bpI] & PROCBND )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical
                            # endif
                            {
                                procEdgeNodes.append(bpI);
                            }

                            continue;
                        }

                        const point newP = newEdgePositionLaplacian(bpI);
                        newPos.append(labelledPoint(bpI, newP));
                    }
                }

                forAll(newPos, i)
                    bMod.moveBoundaryVertexNoUpdate
                    (
                        newPos[i].pointLabel(),
                        newPos[i].coordinates()
                    );
            }

            if( Pstream::parRun() )
            {
                nodeDisplacementSurfaceOptimizerParallel(procBndNodes);
                edgeNodeDisplacementParallel(procEdgeNodes);
            }
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry();

        //- update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while( ++iterationI < 10 );

    //- delete invalid data
    mesh.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
