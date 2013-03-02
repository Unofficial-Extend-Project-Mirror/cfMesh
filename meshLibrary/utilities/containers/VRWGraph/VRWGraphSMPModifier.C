/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Franjo Juretic
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "VRWGraphSMPModifier.H"
#include "labelPair.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

VRWGraphSMPModifier::VRWGraphSMPModifier(VRWGraph& graph)
:
    graph_(graph)
{}

VRWGraphSMPModifier::~VRWGraphSMPModifier()
{}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void VRWGraphSMPModifier::mergeGraphs(const List<VRWGraph>& graphParts)
{
    const label nGraphs = graphParts.size();
    const label nRows = graphParts[0].size();
    forAll(graphParts, i)
    {
        if( nRows != graphParts[i].size() )
            FatalErrorIn
            (
                "inline void Foam::VRWGraph::mergeGraphs(const List<VRWGraph>&)"
            ) << "Cannot merge graphs" << abort(FatalError);
    }
    
    //- find the number of elements in each row
    labelListPMG nElmtsInRow(nRows);
    
    # pragma omp parallel for schedule(static, 1)
    for(label rowI=0;rowI<nRows;++rowI)
    {
        label sum(0);
        for(label i=0;i<nGraphs;++i)
            sum += graphParts[i].sizeOfRow(rowI);
        
        nElmtsInRow[rowI] = sum;
    }
 
    //- set the size of graph
    setSizeAndRowSize(nElmtsInRow);
    
    //- Finally, assemble the merged graph
    # pragma omp parallel for schedule(static, 1)
    for(label rowI=0;rowI<nRows;++rowI)
    {
        forAll(graphParts, i)
        {
            const VRWGraph& gp = graphParts[i];
            for(label j=0;j<gp.sizeOfRow(rowI);++j)
                graph_(rowI, --nElmtsInRow[rowI]) = gp(rowI, j);
        }
    }
}
        
void VRWGraphSMPModifier::reverseAddressing(const VRWGraph& origGraph)
{
    graph_.setSize(0);
    labelListPMG nAppearances;
    
    label nThreads = 3 * omp_get_num_procs();
    if( origGraph.size() < 1000 )
        nThreads = 1;
    
    label minRow(INT_MAX), maxRow(-1);
    List<List<LongList<labelPair> > > dataForOtherThreads(nThreads);
    
    # pragma omp parallel num_threads(nThreads)
    {
        const label threadI = omp_get_thread_num();
        
        List<LongList<labelPair> >& dot = dataForOtherThreads[threadI];
        dot.setSize(nThreads);
        
        //- find min and max entry in the graph
        //- they are used for assigning ranges of values local for each process
        label localMinRow(INT_MAX), localMaxRow(-1);
        # pragma omp for schedule(static)
        forAll(origGraph, rowI)
        {
            forAllRow(origGraph, rowI, i)
            {
                const label entryI = origGraph(rowI, i);
                localMaxRow = Foam::max(localMaxRow, entryI);
                localMinRow = Foam::min(localMinRow, entryI);
            }
        }
        
        ++localMaxRow;
        
        # pragma omp critical
        {
            minRow = Foam::min(minRow, localMinRow);
            maxRow = Foam::max(maxRow, localMaxRow);

            nAppearances.setSize(maxRow);
        }
        
        # pragma omp barrier
        
        //- initialise appearances
        # pragma omp for schedule(static)
        for(label i=0;i<maxRow;++i)
            nAppearances[i] = 0;
        
        # pragma omp barrier

        const label range = (maxRow - minRow) / nThreads + 1;
        const label localMin = minRow + threadI * range;
        const label localMax = Foam::min(localMin + range, maxRow);
        
        //- find the number of appearances of each element in the original graph
        # pragma omp for schedule(static)
        forAll(origGraph, rowI)
        {
            forAllRow(origGraph, rowI, j)
            {
                const label entryI = origGraph(rowI, j);
        
                const label threadNo = (entryI - minRow) / range;
                
                if( threadNo == threadI )
                {
                    ++nAppearances[entryI];
                }
                else
                {
                    dot[threadNo].append(labelPair(entryI, rowI));
                }
            }
        }
        
        # pragma omp barrier

        //- count the appearances which are not local to the processor
        for(label i=0;i<nThreads;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];
          
            forAll(data, j)
                ++nAppearances[data[j].first()];
        }
        
        # pragma omp barrier

        //- allocate graph
        # pragma omp master
        setSizeAndRowSize(nAppearances);
        
        # pragma omp barrier
        
        for(label i=localMin;i<localMax;++i)
        {
            nAppearances[i] = 0;
        }
        
        //- start filling reverse addressing graph
        //- update data from processors with smaller labels
        for(label i=0;i<threadI;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];
      
            forAll(data, j)
            {
                const label entryI = data[j].first();
                graph_(entryI, nAppearances[entryI]++) = data[j].second();
            }
        }

        //- update data local to the processor
        # pragma omp for schedule(static)
        forAll(origGraph, rowI)
        {
            forAllRow(origGraph, rowI, j)
            {
                const label entryI = origGraph(rowI, j);

                if( (entryI >= localMin) && (entryI < localMax) )
                    graph_(entryI, nAppearances[entryI]++) = rowI;
            }
        }
        
        //- update data from the processors with higher labels
        for(label i=threadI+1;i<nThreads;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];
      
            forAll(data, j)
            {
                const label entryI = data[j].first();
                graph_(entryI, nAppearances[entryI]++) = data[j].second();
            }
        }
    }
}

void VRWGraphSMPModifier::optimizeMemoryUsage()
{
    label nThreads = 3 * omp_get_num_procs();
    if( graph_.size() < 1000 )
        nThreads = 1;
    
    DynList<label> nRows, nEntries;
    nRows.setSize(nThreads);
    nEntries.setSize(nThreads);
    
    LongList<rowElement> newRows;
    labelListPMG newData;

    # pragma omp parallel num_threads(nThreads)
    {
        const label threadI = omp_get_thread_num();
        nRows[threadI] = 0;
        nEntries[threadI] = 0;
        
        # pragma omp for schedule(static)
        forAll(graph_.rows_, rowI)
        {
            if( graph_.rows_[rowI].start() == VRWGraph::INVALIDROW )
                continue;
            
            ++nRows[threadI];
            nEntries[threadI] += graph_.rows_[rowI].size();
        }
        
        # pragma omp barrier
        
        # pragma omp master
        {
            //- find the number of rows
            label counter(0);
            forAll(nRows, i)
                counter += nRows[i];
            
            newRows.setSize(counter);
            
            //- find the number of data entries
            counter = 0;
            forAll(nEntries, i)
                counter += nEntries[i];
            
            newData.setSize(counter);
        }
        
        # pragma omp barrier
        
        //- find the starting position for each thread
        label rowStart(0), entryStart(0);
        for(label i=0;i<threadI;++i)
        {
            rowStart += nRows[i];
            entryStart += nEntries[i];
        }
        
        //- copy the data into the location
        # pragma omp for schedule(static)
        forAll(graph_, rowI)
        {
            rowElement& el = newRows[rowStart];
            el.start() += entryStart;
            ++rowStart;
            
            el.size() = graph_.sizeOfRow(rowI);
            forAllRow(graph_, rowI, i)
            {
                newData[entryStart] = graph_(rowI, i);
                ++entryStart;
            }
        }
    }
    
    //- replace the original data with the compressed data
    graph_.rows_.transfer(newRows);
    graph_.data_.transfer(newData);
}

void VRWGraphSMPModifier::operator=(const VRWGraph& og)
{
    graph_.data_.setSize(og.data_.size());
    graph_.rows_.setSize(og.rows_.size());
    
    # pragma omp parallel
    {
        # pragma omp for schedule(static, 1)
        forAll(graph_.data_, i)
            graph_.data_[i] = og.data_[i];
        
        # pragma omp for schedule(static, 1)
        forAll(graph_.rows_, rowI)
            graph_.rows_[rowI] = og.rows_[rowI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
