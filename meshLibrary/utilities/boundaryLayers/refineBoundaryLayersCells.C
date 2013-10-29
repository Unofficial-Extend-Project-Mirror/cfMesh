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
    const cell& c = mesh_.cells()[cellI];

    cellsFromCell.setSize(1);

    DynList<DynList<label, 8>, 10>& cellFaces = cellsFromCell[0];

    forAll(c, fI)
    {
        forAllRow(facesFromFace_, c[fI], cfI)
        {
            const label nfI = facesFromFace_(c[fI], cfI);

            DynList<label, 8> cf;
            cf.setSize(newFaces_.sizeOfRow(nfI));

            forAllRow(newFaces_, nfI, pI)
                cf[pI] = newFaces_(nfI, pI);

            cellFaces.append(cf);
        }
    }
}

//- generate new cells from a hex at a feature edge
void refineBoundaryLayers::generateNewCellsEdgeHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 64>& cellsFromCell
)
{
    const cell& c = mesh_.cells()[cellI];

    cellsFromCell.setSize(1);

    DynList<DynList<label, 4>, 6>& cellFaces = cellsFromCell[0];

    forAll(c, fI)
    {
        forAllRow(facesFromFace_, c[fI], cfI)
        {
            const label nfI = facesFromFace_(c[fI], cfI);

            DynList<label, 4> cf;
            cf.setSize(newFaces_.sizeOfRow(nfI));

            forAllRow(newFaces_, nfI, pI)
                cf[pI] = newFaces_(nfI, pI);

            cellFaces.append(cf);
        }
    }
}

//- generate new cells from a hex at a corner
void refineBoundaryLayers::generateNewCellsCornerHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
)
{
    const cell& c = mesh_.cells()[cellI];

    cellsFromCell.setSize(1);

    DynList<DynList<label, 4>, 6>& cellFaces = cellsFromCell[0];

    forAll(c, fI)
    {
        forAllRow(facesFromFace_, c[fI], cfI)
        {
            const label nfI = facesFromFace_(c[fI], cfI);

            DynList<label, 4> cf;
            cf.setSize(newFaces_.sizeOfRow(nfI));

            forAllRow(newFaces_, nfI, pI)
                cf[pI] = newFaces_(nfI, pI);

            cellFaces.append(cf);
        }
    }
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
    /*
    forAll(nCellsFromCell, cellI)
    {
        Info << "\nCell " << cellI << endl;
        Info << "nCellsFromCell " << nCellsFromCell[cellI] << endl;
        Info << "Ref type " << refType[cellI] << endl;
    }
    Info << "Number of newly generated cells " << nNewCells << endl;
*/
    //- create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    const label numFacesBefore = newFaces_.size();

    //- set the number of cells to the new value
    cellListPMG& cells = meshModifier.cellsAccess();
    label nCells = cells.size();
    cells.setSize(nCells+nNewCells);

    //- start creating new cells
    //- store the information which new cells were generated from
    //- an existing cell
    VRWGraph newCellsFromCell(refType.size());

    VRWGraph pointNewFaces;
    pointNewFaces.reverseAddressing(newFaces_);

    forAll(nCellsFromCell, cellI)
    {
        if( refType[cellI] == 0 )
        {
            newCellsFromCell.append(cellI, cellI);

            cell& c = cells[cellI];

            forAll(c, fI)
            {
                if( facesFromFace_.sizeOfRow(c[fI]) != 1 )
                    FatalError << "Crap!" << abort(FatalError);
                c[fI] = facesFromFace_(c[fI], 0);
            }
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
                        FatalError << "1.Sranje" << abort(FatalError);
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
                        FatalError << "2.Sranje" << abort(FatalError);
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
                        FatalError << "3.Sranje" << abort(FatalError);
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

    //- update cell sets
    mesh_.updateCellSubsets(newCellsFromCell);
    newCellsFromCell.setSize(0);

    //- point-faces addressing is not needed any more
    pointNewFaces.setSize(0);

    //- copy the newFaces to the mesh
    const label nOrigInternalFaces = mesh_.nInternalFaces();
    const label nNewInternalFaces =
        facesFromFace_(mesh_.boundaries()[0].patchStart(), 0);

    //- store internal faces originating from existing faces
    Info << "Copying internal faces " << endl;
    Info << "Original number of internal faces " << nOrigInternalFaces << endl;
    labelListPMG newFaceLabel(newFaces_.size());
    faces.setSize(newFaces_.size());

    label currFace = 0;
    for(label faceI=0;faceI<nOrigInternalFaces;++faceI)
    {
        forAllRow(facesFromFace_, faceI, ffI)
        {
            face& f = faces[currFace];
            newFaceLabel[currFace] = currFace;
            ++currFace;

            const label newFaceI = facesFromFace_(faceI, ffI);

            f.setSize(newFaces_.sizeOfRow(newFaceI));

            forAll(f, pI)
                f[pI] = newFaces_(newFaceI, pI);
        }
    }

    //- store newly-generated internal faces
    Info << "Copying newly generated internal faces" << endl;
    Info << "nNewInternalFaces " << currFace << endl;
    Info << "numFacesBefore " << numFacesBefore << endl;
    Info << "Total number of faces " << newFaces_.size() << endl;

    for(label faceI=numFacesBefore;faceI<newFaces_.size();++faceI)
    {
        newFaceLabel[faceI] = currFace;
        face& f = faces[currFace];
        ++currFace;

        f.setSize(newFaces_.sizeOfRow(faceI));

        forAll(f, pI)
            f[pI] = newFaces_(faceI, pI);
    }

    //- store new boundary faces
    Info << "Copying boundary faces " << endl;
    Info << "currFace " << currFace << endl;
    Info << "Faces size " << faces.size() << endl;
    Info << "Initial number of faces " << facesFromFace_.size() << endl;
    PtrList<writePatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        const label newStart = currFace;
        label nNewFacesInPatch(0);
        for(label fI=0;fI<size;++fI)
        {
            const label faceI = start + fI;

            forAllRow(facesFromFace_, faceI, nfI)
            {
                face& f = faces[currFace];

                //- update the new label
                const label origFaceI = facesFromFace_(faceI, nfI);
                newFaceLabel[origFaceI] = currFace;
                facesFromFace_(faceI, nfI) = currFace;
                ++currFace;

                //- copy the face into the mesh
                f.setSize(newFaces_.sizeOfRow(origFaceI));
                forAll(f, pI)
                    f[pI] = newFaces_(origFaceI, pI);

                ++nNewFacesInPatch;
            }
        }

        //- update patch
        boundaries[patchI].patchStart() = newStart;
        boundaries[patchI].patchSize() = nNewFacesInPatch;
    }

    if( Pstream::parRun() )
    {
        Info << "Copying processor faces" << endl;
        //- copy faces at inter-processor boundaries
        PtrList<writeProcessorPatch>& procBoundaries =
            meshModifier.procBoundariesAccess();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            const label newStart = currFace;
            label nNewFacesInPatch(0);
            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;
                forAllRow(facesFromFace_, faceI, nfI)
                {
                    face& f = faces[currFace];
                    newFaceLabel[faceI] = currFace;

                    //- update the new label
                    const label origFaceI = facesFromFace_(faceI, nfI);
                    facesFromFace_(faceI, nfI) = currFace;
                    ++currFace;

                    //- copy the face into the mesh
                    f.setSize(newFaces_.sizeOfRow(origFaceI));
                    forAll(f, pI)
                        f[pI] = newFaces_(origFaceI, pI);

                    ++nNewFacesInPatch;
                }
            }

            //- update patch
            procBoundaries[patchI].patchStart() = newStart;
            procBoundaries[patchI].patchSize() = nNewFacesInPatch;
        }
    }

    Info << "Faces after refinement " << faces << endl;
    Info << "newFaceLabel " << newFaceLabel << endl;

    //- update face subsets
    Info << "Updating subsets" << endl;
    mesh_.updateFaceSubsets(facesFromFace_);
    facesFromFace_.setSize(0);
    newFaces_.setSize(0);

    //- update cells to match the faces
    Info << "Updating cells to match new faces" << endl;
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            c[fI] = newFaceLabel[c[fI]];
    }

    Info << "Cleaning mesh " << endl;
    meshModifier.clearAll();

    Info << "Finished generating new cells " << endl;
    //::exit(1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

