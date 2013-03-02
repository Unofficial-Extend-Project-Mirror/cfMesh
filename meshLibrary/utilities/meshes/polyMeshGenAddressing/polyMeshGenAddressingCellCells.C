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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCells() const
{
    if( ccPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellCells() const")
            << "cellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellListPMG& cells = mesh_.cells();
        
        const labelList& own = mesh_.owner();
        const labelList& nei = mesh_.neighbour();
        
        //- create the storage
        ccPtr_ = new VRWGraph();
        VRWGraph& cellCellAddr = *ccPtr_;
        
        labelList nNei(cells.size());
        
        const label nThreads = 3 * omp_get_num_procs();
        
        # pragma omp parallel num_threads(nThreads)
        {
            # pragma omp for schedule(static)
            forAll(nNei, i)
                nNei[i] = 0;
            
            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];
                
                DynList<label> neiCells;
                
                forAll(c, fI)
                {
                    label neiCell = own[c[fI]];
                    if( (neiCell == cellI) && (nei[c[fI]] != -1) ) 
                        neiCell = nei[c[fI]];
    
                    if( neiCell != cellI )
                        neiCells.appendIfNotIn(neiCell);
                }
                
                nNei[cellI] = neiCells.size();
            }
            
            # pragma omp barrier
            
            # pragma omp master
            VRWGraphSMPModifier(cellCellAddr).setSizeAndRowSize(nNei);
            
            # pragma omp barrier
            
            //- fill the graph with data
            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];
                
                DynList<label> neiCells;
                
                forAll(c, fI)
                {
                    label neiCell = own[c[fI]];
                    if( (neiCell == cellI) && (nei[c[fI]] != -1) ) 
                        neiCell = nei[c[fI]];
    
                    if( neiCell != cellI )
                        neiCells.appendIfNotIn(neiCell);
                }
                
                cellCellAddr.setRow(cellI, neiCells);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::cellCells() const
{
    if( !ccPtr_ )
        calcCellCells();

    return *ccPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
