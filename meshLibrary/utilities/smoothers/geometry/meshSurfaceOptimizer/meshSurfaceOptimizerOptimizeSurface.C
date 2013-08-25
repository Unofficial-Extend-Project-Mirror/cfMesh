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
#include "labelledPoint.H"

#include <map>
#include <omp.h>

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
    this->triangles();
    this->pointTriangles();

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
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
