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

#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "labelListPMG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfacePartitioner::calculatePartitionAddressing()
{
    calculateCornersAndAddressing();
    
    calculatePartitionPartitions();
    
    calculateEdgePartitions();
        
    calculatePartitionsToEdgePartitions();
    
    calculateEdgePartitionsToCorners();

}
    
void triSurfacePartitioner::calculateCornersAndAddressing()
{
    const labelListList& pointFaces = surface_.pointFaces();
    const edgeList& edges = surface_.edges();
    const labelListList& edgeFaces = surface_.edgeFaces();
    
    //- find the number of feature edges connected to each surface node
    List<direction> nEdgesAtNode(surface_.points().size(), direction(0));
    forAll(edgeFaces, eI)
    {
        if( edgeFaces[eI].size() != 2 )
            continue;
        
        const label sPatch = surface_[edgeFaces[eI][0]].region();
        const label ePatch = surface_[edgeFaces[eI][1]].region();
        
        if( sPatch != ePatch )
        {
            const edge& e = edges[eI];
            ++nEdgesAtNode[e.start()];
            ++nEdgesAtNode[e.end()];
        }
    }
    
    //- count the number of feature edges connected to each surface point
    //- corners must have 3 or more edges attached to them
    label nCorners(0);
    forAll(nEdgesAtNode, pI)
    {
        if( nEdgesAtNode[pI] < direction(3) )
            continue;
        
        ++nCorners;
    }
    
    corners_.setSize(nCorners);
    cornerPatches_.setSize(nCorners);
    nCorners = 0;
    
    //- store corner data
    DynList<label> patches(10);
    forAll(pointFaces, pointI)
    {
        if( nEdgesAtNode[pointI] < direction(3) )
            continue;
        
        const labelList& pf = pointFaces[pointI];
        
        patches.clear();
        forAll(pf, pfI)
            patches.appendIfNotIn(surface_[pf[pfI]].region());
        
        corners_[nCorners] = pointI;
        cornerPatches_[nCorners] = patches;
        ++nCorners;
    }
}
    
void triSurfacePartitioner::calculatePartitionPartitions()
{
    const labelListList& edgeFaces = surface_.edgeFaces();
    
    forAll(edgeFaces, eI)
    {
        if( edgeFaces[eI].size() != 2 )
        {
            Warning << "Surface is not a manifold!!" << endl;
            continue;
        }
        
        const label sPatch = surface_[edgeFaces[eI][0]].region();
        const label ePatch = surface_[edgeFaces[eI][1]].region();
        
        if( sPatch != ePatch )
        {
            partitionPartitions_[sPatch].insert(ePatch);
            partitionPartitions_[ePatch].insert(sPatch);
        }
    }
}

void triSurfacePartitioner::calculateEdgePartitions()
{
    const edgeList& edges = surface_.edges();
    const labelListList& pointEdges = surface_.pointEdges();
    const labelListList& edgeFaces = surface_.edgeFaces();
    
    //- make all feature edges
    boolList featureEdge(edgeFaces.size(), false);
    DynList<label> parts(10);
    forAll(edgeFaces, eI)
    {
        const labelList& eFaces = edgeFaces[eI];
        
        parts.clear();
        forAll(eFaces, efI)
            parts.appendIfNotIn(surface_[eFaces[efI]].region());
        
        if( parts.size() > 1 )
            featureEdge[eI] = true;
    }
    
    //- create a set containing corners for fast searching
    labelHashSet corners;
    forAll(corners_, i)
        corners.insert(corners_[i]);
    
    edgePartitions_.setSize(edgeFaces. size());
    edgePartitions_ = -1;
    
    label nPartitions(0);
    forAll(featureEdge, eI)
    {
        if( !featureEdge[eI] )
            continue;
        
        labelListPMG front;
        front.append(eI);
        edgePartitions_[eI] = nPartitions;
        featureEdge[eI] = false;
        
        while( front.size() )
        {
            const label eLabel = front.removeLastElement();
            const edge& e = edges[eLabel];
            
            for(label pI=0;pI<2;++pI)
            {
                const label pointI = e[pI];
                
                if( corners.found(pointI) )
                    continue;
                
                const labelList& pEdges = pointEdges[pointI];
                forAll(pEdges, peI)
                {
                    const label eJ = pEdges[peI];
                    if( featureEdge[eJ] )
                    {
                        edgePartitions_[eJ] = nPartitions;
                        featureEdge[eJ] = false;
                        front.append(eJ);
                    }
                }
            }
        }
        
        ++nPartitions;
    }
    
    Info << nPartitions << " edge partitions found!" << endl;
    
    edgePartitionEdgePartitions_.clear();
    edgePartitionEdgePartitions_.setSize(nPartitions);
}

void triSurfacePartitioner::calculatePartitionsToEdgePartitions()
{
    const labelListList& edgeFaces = surface_.edgeFaces();
    
    DynList<label> partitions(10);
    forAll(edgeFaces, eI)
    {
        if( edgePartitions_[eI] < 0 )
            continue;
        
        const labelList& eFaces = edgeFaces[eI];
        
        partitions.clear();
        forAll(eFaces, efI)
            partitions.appendIfNotIn(surface_[eFaces[efI]].region());
        
        forAll(partitions, i)
        {
            const label partI = partitions[i];
            for(label j=i+1;j<partitions.size();++j)
            {
                const label partJ = partitions[j];
                
                const std::pair<label, label> pp
                (
                    Foam::min(partI, partJ),
                    Foam::max(partI, partJ)
                );
                
                partitionsEdgeParts_[pp].insert(edgePartitions_[eI]);
            }
        }
    }
}

void triSurfacePartitioner::calculateEdgePartitionsToCorners()
{
    const labelListList& pointEdges = surface_.pointEdges();
    
    forAll(corners_, cornerI)
    {
        DynList<label> edgePartitionsAtCorner(10);
        const labelList& pEdges = pointEdges[corners_[cornerI]];
        
        forAll(pEdges, peI)
            edgePartitionsAtCorner.appendIfNotIn(edgePartitions_[pEdges[peI]]);
        
        forAll(edgePartitionsAtCorner, i)
        {
            const label epI = edgePartitionsAtCorner[i];
            if( epI < 0 )
                continue;
            for(label j=i+1;j<edgePartitionsAtCorner.size();++j)
            {
                const label epJ = edgePartitionsAtCorner[j];
                if( epJ < 0 )
                    continue;
                
                std::pair<label, label> ep
                (
                    Foam::min(epI, epJ),
                    Foam::max(epI, epJ)
                );
                
                //- create edgepartition - edge partitions addressing
                edgePartitionEdgePartitions_[ep.first].insert(ep.second);
                edgePartitionEdgePartitions_[ep.second].insert(ep.first);
                
                edgePartitionsCorners_[ep].insert(cornerI);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
