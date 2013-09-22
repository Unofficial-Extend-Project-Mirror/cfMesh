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

#include "error.H"
#include "helperFunctionsFrontalMarking.H"
#include "DynList.H"
#include "labelPair.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#define DEBUGFrontalMarking

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

template<class labelListType, class neiOp, class filterOp>
void frontalMarking
(
    labelListType& result,
    const label startingIndex,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    //- add the starting element
    result.clear();
    result.append(startingIndex);

    //- add the starting element to the front
    labelListPMG front;
    front.append(startingIndex);

    //- store information which element were already visited
    boolList alreadySelected(neighbourCalculator.size(), false);

    //- start with frontal marking
    while( front.size() )
    {
        const label eLabel = front.removeLastElement();

        //- find neighbours of the current element
        DynList<label> neighbours;
        neighbourCalculator(eLabel, neighbours);

        forAll(neighbours, neiI)
        {
            const label nei = neighbours[neiI];

            if( nei < 0 )
                continue;
            if( alreadySelected[nei] )
                continue;

            if( selector(nei) )
            {
                alreadySelected[nei] = true;
                front.append(nei);
                result.append(nei);
            }
        }
    }
}

template<class labelListType, class neiOp, class filterOp>
label groupMarking
(
    labelListType& elementInGroup,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    label nGroups(0);

    elementInGroup.setSize(neighbourCalculator.size());
    elementInGroup = -1;

    VRWGraph neighbouringGroups;

    label nThreads(1);

    # ifdef USE_OMP
    nThreads = 3 * omp_get_num_procs();
    # endif

    # ifdef USE_OMP
    # pragma omp parallel if( neighbourCalculator.size() > 1000 ) \
    num_threads(nThreads)
    # endif
    {
        const label chunkSize = neighbourCalculator.size() / nThreads + 1;

        const label threadI = omp_get_thread_num();

        LongList<std::pair<label, label> > threadCommPairs;

        const label minEl = threadI * chunkSize;

        const label maxEl =
            Foam::min
            (
                neighbourCalculator.size(),
                minEl + chunkSize
            );

        for(label elI=minEl;elI<maxEl;++elI)
        {
            if( elementInGroup[elI] != -1 )
                continue;
            if( !selector(elI) )
                continue;

            label groupI;
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            groupI = nGroups++;

            elementInGroup[elI] = groupI;
            labelListPMG front;
            front.append(elI);

            while( front.size() )
            {
                const label eLabel = front.removeLastElement();

                DynList<label> neighbours;
                neighbourCalculator(eLabel, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( (nei < 0) || (elementInGroup[nei] != -1) )
                        continue;

                    if( (nei < minEl) || (nei >= maxEl) )
                    {
                        //- this is a communication interface between
                        //- two threads
                        threadCommPairs.append(std::make_pair(elI, nei));
                    }
                    else if( selector(nei) )
                    {
                        //- this element is part of the same group
                        elementInGroup[nei] = groupI;
                        front.append(nei);
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- find group to neighbouring groups addressing
        List<DynList<label> > localNeiGroups(nGroups);
        forAll(threadCommPairs, cfI)
        {
            const std::pair<label, label>& lp = threadCommPairs[cfI];
            const label groupI = elementInGroup[lp.first];
            const label neiGroup = elementInGroup[lp.second];

            if( (neiGroup >= nGroups) || (groupI >= nGroups) )
                FatalError << "neiGroup " << neiGroup
                    << " groupI " << groupI << " are >= than "
                    << "nGroups " << nGroups << abort(FatalError);

            if( neiGroup != -1 )
            {
                localNeiGroups[groupI].appendIfNotIn(neiGroup);
                localNeiGroups[neiGroup].appendIfNotIn(groupI);
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            neighbouringGroups.setSize(nGroups);

            forAll(localNeiGroups, groupI)
            {
                const DynList<label>& lGroups = localNeiGroups[groupI];

                neighbouringGroups.appendIfNotIn(groupI, groupI);

                forAll(lGroups, i)
                    neighbouringGroups.appendIfNotIn(groupI, lGroups[i]);
            }
        }
    }

    DynList<label> globalGroupLabel;
    globalGroupLabel.setSize(nGroups);

    bool changed;

    do
    {
        changed = false;

        DynList<label> newGroupLabel;
        newGroupLabel.setSize(globalGroupLabel.size());
        newGroupLabel = -1;

        //- reduce the information about the groups
        forAllReverse(neighbouringGroups, groupI)
        {
            forAllRow(neighbouringGroups, groupI, ngI)
            {
                const label neiGroup = neighbouringGroups(groupI, ngI);

                if( neiGroup >= groupI )
                    continue;

                forAllRow(neighbouringGroups, groupI, i)
                {
                    const label grI = neighbouringGroups(groupI, i);
                    neighbouringGroups.appendIfNotIn(neiGroup, grI);
                }
            }
        }

        label counter(0);
        forAll(neighbouringGroups, groupI)
        {
            if( newGroupLabel[groupI] != -1 )
                continue;

            forAllRow(neighbouringGroups, groupI, ngI)
                newGroupLabel[neighbouringGroups(groupI, ngI)] = counter;

            ++counter;
        }

        if( counter < nGroups )
        {
            changed = true;
            nGroups = counter;
        }

        globalGroupLabel = newGroupLabel;

        if( Pstream::parRun() )
        {
            //- reduce the groups over processors
            labelList nGroupsAtProc(Pstream::nProcs());
            nGroupsAtProc[Pstream::myProcNo()] = nGroups;

            Pstream::gatherList(nGroupsAtProc);
            Pstream::scatterList(nGroupsAtProc);

            label startGroup(0);
            for(label procI=0;procI<Pstream::myProcNo();++procI)
                startGroup += nGroupsAtProc[procI];

            //- find the neighbouring groups
            std::map<label, DynList<label> > neiGroups;

            //- collect groups on other processors
            //- this operator implements the algorithm for exchanging data
            //- over processors and collects information which groups
            //- are connected ovr inter-processor boundaries
            neighbourCalculator.collectGroups
            (
                neiGroups,
                elementInGroup,
                newGroupLabel,
                startGroup
            );
        }

    } while( changed );

    //- set the global group label
    forAll(elementInGroup, elI)
    {
        if( elementInGroup[elI] < 0 )
            continue;

        elementInGroup[elI] = globalGroupLabel[elementInGroup[elI]];
    }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
