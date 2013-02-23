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

#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "labelledScalar.H"

#include "helperFunctionsPar.H"
#include <omp.h>

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::findMappingDistance
(
    const labelListPMG& nodesToMap,
    std::map<label, scalar>& mappingDistance
) const
{
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    
    //- generate search distance for corner nodes
    mappingDistance.clear();
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        
        mappingDistance.insert(std::make_pair(bpI, 0.0));
        std::map<label, scalar>::iterator mIter = mappingDistance.find(bpI);
        
        const point& p = points[bPoints[bpI]];
        forAllRow(pFaces, bpI, pfI)
        {
            const scalar d = magSqr(faceCentres[pFaces(bpI, pfI)] - p);
            mIter->second = Foam::max(mIter->second, d);
        }
        
        //- safety factor
        mIter->second *= 4.0;
    }
    
    if( Pstream::parRun() )
    {
        //- make sure that corner nodesd at parallel boundaries
        //- have the same range in which they accept the corners
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const labelList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        
        //- create the map for exchanging data
        std::map<label, DynList<labelledScalar> > exchangeData;
        const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], DynList<labelledScalar>()));
        
        forAll(nodesToMap, nI)
        {
            const label bpI = nodesToMap[nI];
            
            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;
                
                exchangeData[neiProc].append
                (
                    labelledScalar(globalPointLabel[bpI], mappingDistance[bpI])
                );
            }
        }
        
        //- exchange data between processors
        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);
        
        //- select the maximum mapping distance for processor points
        forAll(receivedData, i)
        {
            const labelledScalar& ls = receivedData[i];
            
            const label bpI = globalToLocal[ls.scalarLabel()];
            
            //- choose the maximum value for the mapping distance
            std::map<label, scalar>::iterator mIter = mappingDistance.find(bpI);
            mIter->second = Foam::max(mIter->second, ls.value());
        }
    }
}
    
void meshSurfaceMapper::mapCorners(const labelListPMG& nodesToMap)
{
    const triSurfacePartitioner& sPartitioner = surfacePartitioner();
    const labelList& surfCorners = sPartitioner.corners();
    const pointField& sPoints = meshOctree_.surface().localPoints();
    const List<DynList<label> >& cornerPatches = sPartitioner.cornerPatches();
    
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& corners = mPart.corners();
    
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPatches = surfaceEngine_.pointPatches();
    
    std::map<label, scalar> mappingDistance;
    findMappingDistance(nodesToMap, mappingDistance);
    
    //- for every corner in the mesh surface find the nearest corner in the
    //- triSurface
    meshSurfaceEngineModifier sMod(surfaceEngine_);
    
    const label size = nodesToMap.size();
    # pragma omp parallel for if( size > 10 )
    for(label cornerI=0;cornerI<size;++cornerI)
    {
        const label bpI = nodesToMap[cornerI];
        if( !corners.found(bpI) )
            FatalErrorIn
            (
                "meshSurfaceMapper::mapCorners(const labelListPMG&)"
            ) << "Trying to map a point that is not a corner"
                << abort(FatalError);
        
        const point& p = points[bPoints[bpI]];
        const scalar maxDist = mappingDistance[bpI];
        
        //- find the nearest position to the given point patches
        point mapPointApprox(p);
        scalar distSqApprox;
        label iter(0);
        while( iter++ < 20 )
        {
            point newP(vector::zero);
            forAllRow(pPatches, bpI, patchI)
            {
                point np;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    pPatches(bpI, patchI),
                    mapPointApprox
                );
                
                newP += np;
            }
            
            newP /= pPatches.sizeOfRow(bpI);
            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDist )
                break;
            
            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);
        
        //- find the nearest triSurface corner for the given corner
        scalar distSq(mappingDistance[bpI]);
        point mapPoint(p);
        forAll(surfCorners, scI)
        {
            const label cornerID = surfCorners[scI];
            const point& sp = sPoints[cornerID];
            
            if( Foam::magSqr(sp - p) < distSq )
            {
                bool store(true);
                const DynList<label>& cPatches = cornerPatches[scI];
                forAllRow(pPatches, bpI, i)
                {
                    if( !cPatches.contains(pPatches(bpI, i)) )
                    {
                        store = false;
                        break;
                    }
                }
                
                if( store )
                {
                    mapPoint = sp;
                    distSq = Foam::magSqr(sp - p);
                }
            }
        }
        
        if( distSq > 1.2 * distSqApprox )
        {
            mapPoint = mapPointApprox;
        }
        
        //- move the point to the nearest corner
        sMod.moveBoundaryVertexNoUpdate(bpI, mapPoint);
    }
    
    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapEdgeNodes(const labelListPMG& nodesToMap)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPatches = surfaceEngine_.pointPatches();
    
    //- find mapping distance for selected vertices
    std::map<label, scalar> mappingDistance;
    findMappingDistance(nodesToMap, mappingDistance);
    
    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();
    
    LongList<parMapperHelper> parallelBndNodes;
    
    meshSurfaceEngineModifier sMod(surfaceEngine_);
    
    //- map point to the nearest vertex on the triSurface
    const label size = nodesToMap.size();
    # pragma omp parallel for shared(parallelBndNodes) if( size > 100 )
    for(label i=0;i<size;++i)
    {
        const label bpI = nodesToMap[i];
        const point& p = points[bPoints[bpI]];
        
        //- find patches at this edge vertex
        DynList<label> patches;
        forAllRow(pPatches, bpI, j)
            patches.append(pPatches(bpI, j));
        
        const scalar maxDist = mappingDistance[bpI];
        
        //- find approximate position of the vertex on the edge
        point mapPointApprox(p);
        scalar distSqApprox;
        label iter(0);
        while( iter++ < 20 )
        {
            point newP(vector::zero);
            forAll(patches, patchI)
            {
                point np;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    patches[patchI],
                    mapPointApprox
                );
                
                newP += np;
            }
            
            newP /= patches.size();
            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDist )
                break;
            
            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);
        
        //- find the nearest vertex on the triSurface feature edge
        point mapPoint;
        scalar distSq;
        meshOctree_.findNearestEdgePoint(p, patches, mapPoint, distSq);
        
        //- use the vertex with the smallest mapping distance
        if( distSq > 1.2 * distSqApprox )
        {
            mapPoint = mapPointApprox;
            distSq = distSqApprox;
        }
        
        //- check if the mapping distance is within the given tolerances
        if( distSq > maxDist )
        {
            //- this indicates possible problems
            //- reduce the mapping distance
            const scalar f = Foam::sqrt(maxDist / distSq);
            distSq = mappingDistance[bpI];
            mapPoint = f * (mapPoint - p) + p;
        }
        
        //- move the point to the nearest edge vertex
        sMod.moveBoundaryVertexNoUpdate(bpI, mapPoint);
        
        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # pragma omp critical
            parallelBndNodes.append
            (
                parMapperHelper
                (
                    mapPoint,
                    distSq,
                    bpI,
                    -1
                )
            );
        }
    }
    
    sMod.updateGeometry(nodesToMap);
    
    mapToSmallestDistance(parallelBndNodes);
}
        
void meshSurfaceMapper::mapCornersAndEdges()
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    labelListPMG selectedPoints;
    forAllConstIter(labelHashSet, cornerPoints, it)
        selectedPoints.append(it.key());
    
    mapCorners(selectedPoints);
    
    selectedPoints.clear();
    const labelHashSet& edgePoints = mPart.edgeNodes();
    forAllConstIter(labelHashSet, edgePoints, it)
        selectedPoints.append(it.key());
    
    mapEdgeNodes(selectedPoints);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
