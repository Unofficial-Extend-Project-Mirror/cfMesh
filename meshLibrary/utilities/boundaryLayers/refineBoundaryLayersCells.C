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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::generateNewCellsPrism
(
    const label cellI,
    DynList<DynList<DynList<label, 8>, 10> >& cellsFromCell
)
{

}

//- generate new cells from a hex at a feature edge
void refineBoundaryLayers::generateNewCellsEdgeHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 64>& cellsFromCell
)
{

}

//- generate new cells from a hex at a corner
void refineBoundaryLayers::generateNewCellsCornerHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
)
{

}

void refineBoundaryLayers::generateNewCells()
{
    labelList nCellsFromCell(mesh_.cells().size(), 1);
    labelList refType(mesh_.cells().size(), 0);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwners = mse.faceOwners();

    //- calculate the number new cells generated from a cell
    forAll(faceOwners, bfI)
    {
        const label cellI = faceOwners[bfI];

        nCellsFromCell[cellI] *= nLayersAtBndFace_[bfI];

        ++refType[cellI];
    }

    label nNewCells(0);
    forAll(nCellsFromCell, cellI)
        nNewCells += (nCellsFromCell[cellI] - 1);

    forAll(nCellsFromCell, cellI)
    {
        Info << "\nCell " << cellI << endl;
        Info << "nCellsFromCell " << nCellsFromCell[cellI] << endl;
        Info << "Ref type " << refType[cellI] << endl;
    }
    Info << "Number of newly generated cells " << nNewCells << endl;

    //- generate new cells
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    cellListPMG& cells = meshModifier.cellsAccess();
    label nCells = cells.size();
    cells.setSize(nCells+nNewCells);

    VRWGraph newCellsFromCell(refType.size());

    VRWGraph pointNewFaces;
    pointNewFaces.reverseAddressing(newFaces_);

    forAll(nCellsFromCell, cellI)
    {
        if( refType[cellI] == 0 )
        {
            newCellsFromCell.append(cellI, cellI);
        }
        else if( refType[cellI] == 1 )
        {
            //- generate new cells from this prism
            DynList<DynList<DynList<label, 8>, 10> > cellsFromCell;
            generateNewCellsPrism(cellI, cellsFromCell);

            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 8>, 10>& nc = cellsFromCell[cI];

                cell& c = cells[cI==0?cellI:nCells++];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 8>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else if( refType[cellI] == 2 )
        {
            //- generate new cell from a hex cell where two layers intersect
            //- generate mostly hex cells
            DynList<DynList<DynList<label, 4>, 6>, 64> cellsFromCell;
            generateNewCellsEdgeHex(cellI, cellsFromCell);

            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 4>, 6>& nc = cellsFromCell[cI];

                cell& c = cells[cI==0?cellI:nCells++];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 4>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else if( refType[cellI] == 3 )
        {
            //- generate new cells from a hex at a corner where three
            //- layers intersect
            //- generate mostly hex cells
            DynList<DynList<DynList<label, 4>, 6>, 256> cellsFromCell;
            generateNewCellsCornerHex(cellI, cellsFromCell);

            //- new points have been generated
            pointNewFaces.setSize(mesh_.points().size());

            //- recognise face cells in the graph of new faces
            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 4>, 6>& nc = cellsFromCell[cI];

                cell& c = cells[cI==0?cellI:nCells++];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 4>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void refineBoundaryLayers::generateNewCells()"
            ) << "Cannot refine boundary layer for cell "
              << cellI << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

