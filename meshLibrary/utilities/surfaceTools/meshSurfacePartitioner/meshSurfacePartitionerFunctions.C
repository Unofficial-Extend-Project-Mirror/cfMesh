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

#include "meshSurfacePartitioner.H"
#include "helperFunctionsPar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfacePartitioner::calculateCornersEdgesAndAddressing()
{
    const labelList& bPoints = meshSurface_.boundaryPoints();
    const labelList& bp = meshSurface_.bp();
    const edgeList& edges = meshSurface_.edges();
    const VRWGraph& edgeFaces = meshSurface_.edgeFaces();
    const labelList& facePatches = meshSurface_.boundaryFacePatches();
    
    label nPatches(0);
    forAll(facePatches, patchI)
        nPatches = Foam::max(nPatches, facePatches[patchI] + 1);
    
    corners_.clear();
    edgeNodes_.clear();
    partitionPartitions_.setSize(nPatches);
    
    nEdgesAtPoint_.clear();
    nEdgesAtPoint_.setSize(bPoints.size(), 0);
    
    forAll(edgeFaces, edgeI)
    {
        if( edgeFaces.sizeOfRow(edgeI) != 2 )
            continue;
        
        const label patch0 = facePatches[edgeFaces(edgeI, 0)];
        const label patch1 = facePatches[edgeFaces(edgeI, 1)];
        
        if( patch0 != patch1 )
        {
            const edge& e = edges[edgeI];
            ++nEdgesAtPoint_[bp[e.start()]];
            ++nEdgesAtPoint_[bp[e.end()]];
            
            partitionPartitions_[patch0].insert(patch1);
            partitionPartitions_[patch1].insert(patch0);
        }
    }
    
    if( Pstream::parRun() )
    {
        const Map<label>& otherFaceAtProc = meshSurface_.otherEdgeFaceAtProc();
        const Map<label>& otherFacePatch = meshSurface_.otherEdgeFacePatch();
        
        //- take into account feature edges at processor boundaries
        forAllConstIter(Map<label>, otherFaceAtProc, it)
        {
            const label beI = it.key();
            
            if( it() <= Pstream::myProcNo() )
                continue;
            if( otherFacePatch[beI] != facePatches[edgeFaces(beI, 0)] )
            {
                const edge& e = edges[beI];
                ++nEdgesAtPoint_[bp[e.start()]];
                ++nEdgesAtPoint_[bp[e.end()]];
            }
        }
        
        //- gather data on all processors
        std::map<label, labelListPMG> exchangeData;
        const DynList<label>& bpNeiProcs = meshSurface_.bpNeiProcs();
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], labelListPMG())
            );
        
        const Map<label>& globalToLocal =
            meshSurface_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = meshSurface_.bpAtProcs();
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();
            
            forAllRow(bpAtProcs, bpI, i)
            {
                const label procI = bpAtProcs(bpI, i);
                
                if( procI == Pstream::myProcNo() )
                    continue;
                
                labelListPMG& dts = exchangeData[procI];
                
                //- exchange data as follows:
                //- 1. global point label
                //- 2. number of feature edges connected to the vertex
                dts.append(it.key());
                dts.append(nEdgesAtPoint_[bpI]);
            }
        }
        
        //- exchange information
        labelListPMG receivedData;
        help::exchangeMap(exchangeData, receivedData);
        
        //- add the edges from other processors to the points
        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nEdges = receivedData[counter++];
            
            nEdgesAtPoint_[bpI] += nEdges;
        }
    }
    
    forAll(nEdgesAtPoint_, bpI)
    {
        if( nEdgesAtPoint_[bpI] > 2 )
        {
            corners_.insert(bpI);
        }
        else if( nEdgesAtPoint_[bpI] == 2 )
        {
            edgeNodes_.insert(bpI);
        }
    }
    
    label counter = corners_.size();
    reduce(counter, sumOp<label>());
    Info << "Found " << counter
        << " corners at the surface of the volume mesh" << endl;
    counter = edgeNodes_.size();
    reduce(counter, sumOp<label>());
    Info << "Found " << counter
        << " edge points at the surface of the volume mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
