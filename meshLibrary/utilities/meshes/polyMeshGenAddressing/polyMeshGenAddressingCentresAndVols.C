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

Class
    polyMeshGenAddressing

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCentresAndVols() const
{
    if( cellCentresPtr_ || cellVolumesPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }
    
    const cellListPMG& cells = mesh_.cells();

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(cells.size());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(cells.size());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);
}

void polyMeshGenAddressing::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    const labelList& own = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    const label nCells = cells.size();
    
    # pragma omp parallel for if( nCells > 1000 )
    for(label cellI=0;cellI<nCells;++cellI)
    {
        const cell& c = cells[cellI];
        
        //- approximate the centre first
        vector cEst(vector::zero);
        forAll(c, fI)
            cEst += fCtrs[c[fI]];
        
        cEst /= c.size();
        
        //- start evaluating the volume and the cell centre
        vector cellCentre(vector::zero);
        scalar cellVol(0.0);
        
        forAll(c, fI)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = (fAreas[c[fI]] & (fCtrs[c[fI]] - cEst));
            
            if( own[c[fI]] != cellI )
                pyr3Vol *= -1.0;
            
            pyr3Vol = Foam::max(pyr3Vol, VSMALL);

            // Calculate face-pyramid centre
            const vector pc = (3.0/4.0)*fCtrs[c[fI]] + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            cellCentre += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            cellVol += pyr3Vol;
        }
        
        cellCtrs[cellI] = cellCentre / cellVol;
        cellVols[cellI] = cellVol / 3.0;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::cellCentres() const
{
    if( !cellCentresPtr_ )
        calcCellCentresAndVols();

    return *cellCentresPtr_;
}

const scalarField& polyMeshGenAddressing::cellVolumes() const
{
    if( !cellVolumesPtr_ )
        calcCellCentresAndVols();

    return *cellVolumesPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
