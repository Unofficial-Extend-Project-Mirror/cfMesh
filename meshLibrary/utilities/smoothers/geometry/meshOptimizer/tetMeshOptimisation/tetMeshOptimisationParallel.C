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
#include "tetMeshOptimisation.H"
#include "partTetMesh.H"
#include "VRWGraph.H"
#include "helperFunctionsPar.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshOptimisation::unifyNegativePoints(boolList& negativeNode) const
{
    //- make sure that processor nodes are selected on all processors
    //- where they exist
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const DynList<label>& neiProcs = tetMesh_.neiProcs();
    const labelListPMG& globalPointLabel = tetMesh_.globalPointLabel();
    const VRWGraph& pProcs = tetMesh_.pointAtProcs();
    const Map<label>& globalToLocal = tetMesh_.globalToLocalPointAddressing();
    const labelListPMG& pAtParallelBoundaries =
        tetMesh_.pointsAtProcessorBoundaries();
    
    std::map<label, labelListPMG> selectedNegativeNodes;
    forAll(neiProcs, procI)
        selectedNegativeNodes.insert
        (
            std::make_pair(neiProcs[procI], labelListPMG())
        );
    
    forAll(pAtParallelBoundaries, i)
    {
        const label nI = pAtParallelBoundaries[i];
        
        if( !negativeNode[nI] )
            continue;
        if( !(smoothVertex[nI] & partTetMesh::PARALLELBOUNDARY) )
            continue;
        
        forAllRow(pProcs, nI, procI)
        {
            const label neiProc = pProcs(nI, procI);
            
            if( neiProc == Pstream::myProcNo() )
                continue;
            
            selectedNegativeNodes[neiProc].append(globalPointLabel[nI]);
        }
    }
    
    labelListPMG receivedNodes;
    help::exchangeMap(selectedNegativeNodes, receivedNodes);
    
    forAll(receivedNodes, i)
        negativeNode[globalToLocal[receivedNodes[i]]] = true;
}

void tetMeshOptimisation::exchangeData
(
    std::map<label, DynList<parPartTet> >& m,
    boolList* negativeNodePtr
)
{
    const DynList<label>& neiProcs = tetMesh_.neiProcs();
    const labelListPMG& globalPointLabel = tetMesh_.globalPointLabel();
    const VRWGraph& pProcs = tetMesh_.pointAtProcs();
    const Map<label>& globalToLocal = tetMesh_.globalToLocalPointAddressing();
    const labelListPMG& pAtParallelBoundaries =
        tetMesh_.pointsAtProcessorBoundaries();
    
    const LongList<point>& points = tetMesh_.points();
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const LongList<partTet>& tets = tetMesh_.tets();
    
    //- create the map holding data which will be sent to other procs
    std::map<label, LongList<parPartTet> > exchangeData;
    forAll(neiProcs, procI)
        exchangeData.insert
        (
            std::make_pair(neiProcs[procI], LongList<parPartTet>())
        );
    
    //- create storage in the m map
    m.clear();
    forAll(pAtParallelBoundaries, i)
    {
        const label pI = pAtParallelBoundaries[i];
        if( !(smoothVertex[pI] & partTetMesh::SMOOTH) )
            continue;
        if( negativeNodePtr && !(*negativeNodePtr)[pI] )
            continue;
        
        m.insert(std::make_pair(pI, DynList<parPartTet>(12)));
    }
    
    //- store local data into the maps
    forAll(tets, tetI)
    {
        const partTet& tet = tets[tetI];
        
        DynList<label> sendToProcs(3);
        for(label i=0;i<4;++i)
        {
            const label vI = tet[i];
            
            if( pProcs.sizeOfRow(vI) == 0 )
                continue;
            if( !(smoothVertex[vI] & partTetMesh::SMOOTH) )
                continue;
            if( negativeNodePtr && !(*negativeNodePtr)[vI] )
                continue;
            
            //- add processor labels
            forAllRow(pProcs, vI, procI)
                sendToProcs.appendIfNotIn(pProcs(vI, procI));
            
            //- add data into the map of proc bnd points
            DynList<parPartTet>& data = m[vI];
            data.append
            (
                parPartTet
                (
                    labelledPoint(globalPointLabel[tet[0]], points[tet[0]]),
                    labelledPoint(globalPointLabel[tet[1]], points[tet[1]]),
                    labelledPoint(globalPointLabel[tet[2]], points[tet[2]]),
                    labelledPoint(globalPointLabel[tet[3]], points[tet[3]])
                )
            );
        }
        
        if( sendToProcs.size() != 0 )
        {
            //- add data into the map 
            forAll(sendToProcs, procI)
            {
                const label neiProc = sendToProcs[procI];
                if( neiProc == Pstream::myProcNo() )
                    continue;
                
                exchangeData[neiProc].append
                (
                    parPartTet
                    (
                        labelledPoint(globalPointLabel[tet[0]], points[tet[0]]),
                        labelledPoint(globalPointLabel[tet[1]], points[tet[1]]),
                        labelledPoint(globalPointLabel[tet[2]], points[tet[2]]),
                        labelledPoint(globalPointLabel[tet[3]], points[tet[3]])
                    )
                );
            }
        }
    }
    
    //- exchange data with other processors
    LongList<parPartTet> receivedData;
    help::exchangeMap(exchangeData, receivedData);
        
    forAll(receivedData, i)
    {
        const parPartTet& tet = receivedData[i];
        
        for(label i=0;i<4;++i)
        {
            const label gpI = tet[i].pointLabel();
            
            if(
                globalToLocal.found(gpI) &&
                (smoothVertex[globalToLocal[gpI]] &partTetMesh::SMOOTH)
            )
            {
                DynList<parPartTet>& data = m[globalToLocal[gpI]];
                data.append(tet);
            }
        }
    }
}

void tetMeshOptimisation::updateBufferLayerPoints()
{
    const LongList<point>& points = tetMesh_.points();
    const labelListPMG& bufferLayerPoints = tetMesh_.bufferLayerPoints();
    const VRWGraph& pProcs = tetMesh_.pointAtProcs();
    const labelListPMG& globalPointLabel = tetMesh_.globalPointLabel();
    const Map<label>& globalToLocal = tetMesh_.globalToLocalPointAddressing();
    const DynList<label>& neiProcs = tetMesh_.neiProcs();
    
    //- create the map
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<labelledPoint>())
        );
    
    //- add points into the map
    forAll(bufferLayerPoints, pI)
    {
        const label pointI = bufferLayerPoints[pI];
        
        forAllRow(pProcs, pointI, i)
        {
            const label neiProc = pProcs(pointI, i);
            
            if( neiProc == Pstream::myProcNo() )
                continue;
            
            exchangeData[neiProc].append
            (
                labelledPoint(globalPointLabel[pointI], points[pointI])
            );
        }
    }
    
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);
    
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        
        tetMesh_.updateVertex
        (
            globalToLocal[lp.pointLabel()],
            lp.coordinates()
        );
    }
}

void tetMeshOptimisation::unifyCoordinatesParallel
(
    const boolList* negativeNodePtr
)
{
    const LongList<point>& points = tetMesh_.points();
    const DynList<label>& neiProcs = tetMesh_.neiProcs();
    const VRWGraph& pProcs = tetMesh_.pointAtProcs();
    const Map<label>& globalToLocal = tetMesh_.globalToLocalPointAddressing();
    const labelListPMG& globalPointLabel = tetMesh_.globalPointLabel();
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const labelListPMG& pAtParallelBoundaries =
        tetMesh_.pointsAtProcessorBoundaries();
    
    //- create the map
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, procI)
        exchangeData.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledPoint>())
        );
    
    //- fill in the data
    std::map<label, labelledPoint> parallelBndPoints;
    forAll(pAtParallelBoundaries, i)
    {
        const label pI = pAtParallelBoundaries[i];
        
        if( !(smoothVertex[pI] & partTetMesh::PARALLELBOUNDARY) )
            continue;
        
        parallelBndPoints[pI] = labelledPoint(1, points[pI]);
        
        if( negativeNodePtr && !(*negativeNodePtr)[pI] )
            continue;
        
        forAllRow(pProcs, pI, procI)
        {
            const label neiProc = pProcs(pI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;
            
            exchangeData[neiProc].append
            (
                labelledPoint(globalPointLabel[pI], points[pI])
            );
        }
    }
    
    //- send points to other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);
    
    //- gather the data
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        
        std::map<label, labelledPoint>::iterator iter =
            parallelBndPoints.find(globalToLocal[lp.pointLabel()]);
        
        if( iter == parallelBndPoints.end() )
            continue;
        
        ++iter->second.pointLabel();
        iter->second.coordinates() += lp.coordinates();
    }
    
    //- move the point to the averaged position
    for
    (
        std::map<label, labelledPoint>::iterator it=parallelBndPoints.begin();
        it!=parallelBndPoints.end();
        ++it
    )
    {
        const label pI = it->first;
        
        const point newP = it->second.coordinates() / it->second.pointLabel();
        tetMesh_.updateVertex(pI, newP);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
