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

template<class ListType>
void VRWGraphSMPModifier::setSizeAndRowSize(const ListType& s)
{
    graph_.rows_.setSize(s.size());
    
    label nThreads = 3 * omp_get_num_procs();
    if( s.size() < 1000 )
        nThreads = 1;
    
    label nEntries(0);
    DynList<label> procEntries;
    procEntries.setSize(nThreads);
    
    # pragma omp parallel num_threads(nThreads)
    {
        label& nLocalEntries = procEntries[omp_get_thread_num()];
        nLocalEntries = 0;
        
        # pragma omp for schedule(static)
        forAll(s, i)
            nLocalEntries += s[i];
        
        # pragma omp critical
        nEntries += nLocalEntries;
        
        # pragma omp barrier
        
        # pragma omp master
        {
            graph_.data_.setSize(nEntries);
        }
        
        # pragma omp barrier
        
        label start(0);
        for(label i=0;i<omp_get_thread_num();++i)
            start += procEntries[i];
        
        # pragma omp for schedule(static)
        forAll(s, i)
        {
            graph_.rows_[i].start() = start;
            graph_.rows_[i].size() = s[i];
            start += s[i];
        }
    }
}

template<class GraphType>
void VRWGraphSMPModifier::reverseAddressing(const GraphType& origGraph)
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
            forAll(origGraph[rowI], i)
            {
                const label entryI = origGraph[rowI][i];
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
            forAll(origGraph[rowI], j)
            {
                const label entryI = origGraph[rowI][j];
        
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
            forAll(origGraph[rowI], j)
            {
                const label entryI = origGraph[rowI][j];

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

template<class ListType, class GraphType>
void VRWGraphSMPModifier::reverseAddressing
(
    const ListType& mapper,
    const GraphType& origGraph
)
{
    ListType nAppearances;
    
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
            forAll(origGraph[rowI], i)
            {
                const label entryI = mapper[origGraph[rowI][i]];
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
            forAll(origGraph[rowI], i)
            {
                const label entryI = mapper[origGraph[rowI][i]];
        
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
        {
            setSizeAndRowSize(nAppearances);
        }
        
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
            forAll(origGraph[rowI], j)
            {
                const label entryI = mapper[origGraph[rowI][j]];

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

template<class ListType>
void VRWGraphSMPModifier::reverseAddressing
(
    const ListType& mapper,
    const VRWGraph& origGraph
)
{
    ListType nAppearances;
    
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
                const label entryI = mapper[origGraph(rowI, i)];
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
            forAllRow(origGraph, rowI, i)
            {
                const label entryI = mapper[origGraph(rowI, i)];
        
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
        {
            setSizeAndRowSize(nAppearances);
        }
        
        # pragma omp barrier
        
        for(label i=localMin;i<localMax;++i)
            nAppearances[i] = 0;
        
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
                const label entryI = mapper[origGraph(rowI, j)];

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
