/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "triangulateBoundaryLayerBaseFaces.H"
#include "decomposeFaces.H"
#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateBoundaryLayerBaseFaces::classifyMeshFaces()
{
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    faceType_.setSize(faces.size());
    ownerInColumn_.setSize(faces.size());
    neiInColumn_.setSize(faces.size());

    //- find cells in column
    labelList cellInColumn(mesh_.cells().size(), -1);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(layerCellsInColumn_, scI)
    {
        forAllRow(layerCellsInColumn_, scI, ccI)
            cellInColumn[layerCellsInColumn_(scI, ccI)] = scI;
    }

    //- analyse internal faces
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    for(label faceI=0;faceI<mesh_.nInternalFaces();++faceI)
    {
        const label colOwn = cellInColumn[owner[faceI]];
        const label colNei = cellInColumn[neighbour[faceI]];

        ownerInColumn_[faceI] = colOwn;
        neiInColumn_[faceI] = colNei;

        direction fType(NONE);

        if( (colOwn != -1) && (colOwn == colNei) )
        {
            fType |= BASEFACE;
            fType |= FACEINCOLUMN;
        }
        else if( (colOwn != -1) && (colNei != -1) && (colOwn != colNei) )
        {
            fType |= FACEBETWEENCOLUMNS;
        }
        else if
        (
            ((colOwn != -1) && (colNei == -1)) ||
            ((colOwn == -1) && (colNei != -1))
        )
        {
            fType |= BASEFACE;
            fType |= BASEFACEINSIDEBND;
        }
    }

    //- analyse boundary faces
    label startFace = mesh_.nInternalFaces();
    label endFace = faces.size();
    if( Pstream::parRun() )
        endFace = mesh_.procBoundaries()[0].patchStart();

    for(label faceI=startFace;faceI<endFace;++faceI)
    {
        const label colOwn = cellInColumn[owner[faceI]];
        ownerInColumn_[faceI] = colOwn;
        neiInColumn_[faceI] = -1;

        direction fType(NONE);

        if( colOwn != -1 )
        {
            fType |= BASEFACE;
            fType |= BASEFACEATBND;
        }

        faceType_[faceI] = fType;
    }

    //- analse faces at inter-processor boundaries
    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            labelList patchNeiColumn(procBoundaries[patchI].patchSize());

            const label start = procBoundaries[patchI].patchStart();

            forAll(patchNeiColumn, fI)
                patchNeiColumn[fI] = cellInColumn[owner[start+fI]];

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                patchNeiColumn.byteSize()
            );

            toOtherProc << patchNeiColumn;
        }

        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList neiColumn;
            fromOtherProc >> neiColumn;

            const label start = procBoundaries[patchI].patchStart();
            forAll(neiColumn, fI)
            {
                const label faceI = start + fI;

                const label colOwn = cellInColumn[owner[faceI]];
                const label colNei = neiColumn[fI];

                ownerInColumn_[faceI] = colOwn;
                neiInColumn_[faceI] = colNei;

                direction fType(NONE);

                if( (colOwn != -1) && (colOwn == colNei) )
                {
                    fType |= BASEFACE;
                    fType |= FACEINCOLUMN;
                }
                else if
                (
                    (colOwn != -1) && (colNei != -1) && (colOwn != colNei)
                )
                {
                    fType |= FACEBETWEENCOLUMNS;
                }
                else if
                (
                    ((colOwn != -1) && (colNei == -1)) ||
                    ((colOwn == -1) && (colNei != -1))
                )
                {
                    fType |= BASEFACE;
                    fType |= BASEFACEINSIDEBND;
                }
            }
        }
    }
}

bool triangulateBoundaryLayerBaseFaces::createNewPoints()
{
    boolList splitColumn(layerCellsInColumn_.size(), false);
    label nInvalidColumns(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50) \
    reduction(+:nInvalidColumns)
    # endif
    forAll(layerCellsInColumn_, colI)
    {
        forAllRow(layerCellsInColumn_, colI, i)
        {
            if( invertedCell_[layerCellsInColumn_(colI, i)] )
            {
                ++nInvalidColumns;
                splitColumn[colI] = true;
                break;
            }
        }
    }

    //- check if there exist invalid cells in the boundary layer
    if( nInvalidColumns == 0 )
        return false;

    //- find faces requiring new points
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    //- mark layer
    boolList decomposeFace(owner.size(), false);
    List<DynList<label> > baseFacesInColumns(layerCellsInColumn_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(splitColumn, scI)
    {
        if( splitColumn[scI] )
        {
            //- find faces in this column which shall be triangulated
            labelHashSet cellsInColumn(layerCellsInColumn_.sizeOfRow(scI));

            forAllRow(layerCellsInColumn_, scI, ccI)
                cellsInColumn.insert(layerCellsInColumn_(scI, ccI));

            //- select faces shared between two column cells
            forAllConstIter(labelHashSet, cellsInColumn, it)
            {
                const cell& c = cells[it.key()];

                forAll(c, fI)
                {
                    if( faceType_[c[fI]] & BASEFACE )
                    {
                        decomposeFace[c[fI]] = true;
                    }
                }
            }
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(splitColumn, scI)
    {
        if( splitColumn[scI] )
        {
            //- find faces in this column which shall be triangulated
            labelHashSet cellsInColumn(layerCellsInColumn_.sizeOfRow(scI));

            DynList<label>& basesInColumn = baseFacesInColumns[scI];
            basesInColumn.clear();

            forAllRow(layerCellsInColumn_, scI, ccI)
                cellsInColumn.insert(layerCellsInColumn_(scI, ccI));

            //- select faces shared between two column cells
            bool validDecomposition(true);
            forAllConstIter(labelHashSet, cellsInColumn, it)
            {
                label nSelectedFaces(0);

                const cell& c = cells[it.key()];

                forAll(c, fI)
                {
                    if( decomposeFace[c[fI]] )
                    {
                        ++nSelectedFaces;

                        if( faceType_[c[fI]] & BASEFACEATBND )
                            basesInColumn.append(c[fI]);
                    }
                }

                if( nSelectedFaces > 2 )
                {
                    //- the cell shall be split in more than one direction
                    //- this situation is not allowed
                    validDecomposition = false;
                }
            }

            if( validDecomposition )
            {
                //- find and sort all base faces
                labelHashSet processedCell(cellsInColumn.size());

                bool found;
                do
                {
                    found = false;

                    const label currBaseFace = basesInColumn.lastElement();

                    //- find the cell in the boundary layer attached to
                    //- the current base face. Do not select previously
                    //- processed cells
                    label cLabel(-1);
                    if( !processedCell.found(owner[currBaseFace]) )
                    {
                        cLabel = owner[currBaseFace];
                    }
                    else if
                    (
                        neighbour[currBaseFace] >= 0 &&
                        !processedCell.found(neighbour[currBaseFace])
                    )
                    {
                        cLabel = neighbour[currBaseFace];
                    }

                    if( cellsInColumn.found(cLabel) )
                    {
                        //- find other face marked for decomposition
                        const cell& c = cells[cLabel];

                        forAll(c, fI)
                        {
                            if( c[fI] != currBaseFace )
                                continue;

                            if( faceType_[c[fI]] & BASEFACE  )
                            {
                                basesInColumn.append(c[fI]);
                                processedCell.insert(cLabel);
                                found = true;
                                break;
                            }
                        }
                    }
                } while( found );
            }
            else
            {
                //- do not decompose the cells assigned to this column
                //- it is not a column after all. It is a matrix of cells
                //- at exitting faces and/or corners
                forAllConstIter(labelHashSet, cellsInColumn, it)
                {
                    const cell& c = cells[it.key()];

                    forAll(c, fI)
                        decomposeFace[c[fI]] = false;
                }

                splitColumn[scI] = false;
            }
        }
    }

    //- create points located at face centres
    pointFieldPMG& points = mesh_.points();
    faceCentreLabel_.setSize(faces.size());
    faceCentreLabel_ = -1;

    forAll(decomposeFace, faceI)
    {
        if( decomposeFace[faceI] )
        {
            const point newP = faces[faceI].centre(points);

            faceCentreLabel_[faceI] = points.size();
            points.append(newP);
        }
    }

    //- smooth variation of face-centre points


    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(baseFacesInColumns, scI)
    {
        const DynList<label>& basesInColumn = baseFacesInColumns[scI];

        if( basesInColumn.size() < 2 )
            continue;

        DynList<point> faceCentres(basesInColumn.size());

        //- get face centres
        forAll(basesInColumn, baseI)
        {
            faceCentres[baseI] = points[faceCentreLabel_[basesInColumn[baseI]]];
        }

        //- find the thickness of each layer
        DynList<scalar> layerThickness(layerCellsInColumn_.sizeOfRow(scI));
        for(label layerI=1;layerI<layerThickness.size();++layerI)
        {
            layerThickness[layerI] =
                mag(faceCentres[layerI] - faceCentres[0]);
        }

        //- position all new vertices on a line
        const vector vec = faceCentres.lastElement() - faceCentres[0];
        const scalar magVec = mag(vec) + VSMALL;

        for(label layerI=1;layerI<layerThickness.size();++layerI)
        {
            const point newP =
                faceCentres[0] + layerThickness[layerI]/magVec * vec;

//            const face& bf = faces[facesFromFace(basesInColumn[layerI], 0)];
//            points[bf[2]] = newP;
        }
    }

    //- create missing quad faces as they are internal faces

    //- start by counting the number of new faces in each cell
//    label nNewFaces(0);
//    label nNewCells(0);

//    # ifdef USE_OMP
//    # pragma omp parallel for schedule(dynamic, 50) \
    reduction(+:nNewCells) reduction(+:nNewFacess)
//    # endif
//    forAll(splitColumn, scI)
//    {
//        if( splitColumn[scI] )
//        {

//        }
//    }

    return true;
}

void triangulateBoundaryLayerBaseFaces::createNewFacesAndCells()
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
