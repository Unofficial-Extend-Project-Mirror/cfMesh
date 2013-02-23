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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "DynList.H"
#include "VRWGraphSMPModifier.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellEdges() const
{
    if( cePtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellEdges() const")
            << "cellEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellListPMG& cells = mesh_.cells();
        const VRWGraph& fe = faceEdges();

        cePtr_ = new VRWGraph();
        VRWGraph& cellEdgeAddr = *cePtr_;
        
        labelList nEdges(cells.size());
        
        const label nThreads = 3 * omp_get_num_procs();

        # pragma omp parallel num_threads(nThreads) if( cells.size() > 10000 )
        {
            # pragma omp for schedule(static)
            forAll(nEdges, i)
                nEdges[i] = 0;
            
            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];
                
                DynList<label, 32> cEdges;
                
                forAll(c, fI)
                {
                    const label faceI = c[fI];
                    
                    forAllRow(fe, faceI, eI)
                        cEdges.appendIfNotIn(fe(faceI, eI));
                }
                
                nEdges[cellI] = cEdges.size();
            }
            
            # pragma omp barrier
            
            # pragma omp master
            VRWGraphSMPModifier(cellEdgeAddr).setSizeAndRowSize(nEdges);
            
            # pragma omp barrier
            
            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];
                
                DynList<label, 32> cEdges;
                
                forAll(c, fI)
                {
                    const label faceI = c[fI];
                    
                    forAllRow(fe, faceI, eI)
                        cEdges.appendIfNotIn(fe(faceI, eI));
                }
                
                cellEdgeAddr.setRow(cellI, cEdges);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::cellEdges() const
{
    if( !cePtr_ )
        calcCellEdges();

    return *cePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
