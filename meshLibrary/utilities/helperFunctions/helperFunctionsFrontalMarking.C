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

    forAll(neighbourCalculator, elI)
    {
        if( elementInGroup[elI] != -1 )
            continue;
        if( !selector(elI) )
            continue;

        elementInGroup[elI] = nGroups;
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

                if( selector(nei) )
                {
                    elementInGroup[nei] = nGroups;
                    front.append(nei);
                }
            }
        }

        //- increment the number of groups
        ++nGroups;
    }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
