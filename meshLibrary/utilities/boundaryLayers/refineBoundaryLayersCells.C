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

//#define DEBUGLayer

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
    const labelList& owner = mesh_.owner();

    # ifdef DEBUGLayer
    Pout << "New cells from cell " << cellI << endl;
    # endif

    const label startBoundary = mesh_.boundaries()[0].patchStart();

    //- find the number of lyers for this cell
    label nLayers(1), baseFace(-1);
    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace_.size()) )
            continue;

        if( nLayersAtBndFace_[bfI] < 2 )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        nLayers = nLayersAtBndFace_[bfI];
        baseFace = fI;
    }

    # ifdef DEBUGLayer
    Pout << "Number of layers " << nLayers << endl;
    Pout << "Base face " << baseFace << " has points "
         << mesh_.faces()[c[baseFace]] << endl;
    forAll(c, fI)
    {
        Pout << "Faces from face " << fI << " are "
             << facesFromFace_[c[fI]] << endl;

        forAllRow(facesFromFace_, c[fI], i)
            Pout << "Face " << facesFromFace_(c[fI], i)
                 << " is " << newFaces_[facesFromFace_(c[fI], i)] << endl;
    }
    # endif

    //- set the number of layers
    cellsFromCell.setSize(nLayers);

    //- distribute existing faces into new cells
    label otherBaseFace(-1);
    forAll(c, fI)
    {
        if( fI == baseFace )
        {
            const label faceI = facesFromFace_(c[fI], 0);
            DynList<label, 8> f;
            f = newFaces_[faceI];
            cellsFromCell[nLayers-1].append(f);
        }
        else if( facesFromFace_.sizeOfRow(c[fI]) == 1 )
        {
            const label faceI = facesFromFace_(c[fI], 0);
            otherBaseFace = fI;
            DynList<label, 8> f;
            f = newFaces_[faceI];
            cellsFromCell[0].append(f);
        }
        else
        {
            forAllRow(facesFromFace_, c[fI], cfI)
            {
                const label nfI = facesFromFace_(c[fI], cfI);

                DynList<label, 8> cf;
                cf = newFaces_[nfI];

                if( owner[c[fI]] != cellI )
                    cf = help::reverseFace(cf);

                cellsFromCell[Foam::max(nLayers-1-cfI, 0)].append(cf);
            }
        }
    }

    //- generate missing faces
    const faceListPMG& faces = mesh_.faces();
    const face& bf = faces[c[baseFace]];
    const face& obf = faces[c[otherBaseFace]];
    for(label layerI=1;layerI<nLayers;++layerI)
    {
        //- create new face from points at the same height
        DynList<label, 8> cf;
        forAll(bf, pI)
        {
            const label pointI = bf[pI];

            # ifdef DEBUGLayer
            Pout << "Split edges at point " << pointI << " are "
                 << splitEdgesAtPoint_[pointI] << endl;
            # endif

            label seI(-1);
            if( splitEdgesAtPoint_.sizeOfRow(pointI) == 1 )
            {
                seI = splitEdgesAtPoint_(pointI, 0);
            }
            else
            {
                forAllRow(splitEdgesAtPoint_, pointI, sepI)
                {
                    const label seJ = splitEdgesAtPoint_(pointI, sepI);
                    const edge& se = splitEdges_[seJ];

                    if( obf.which(se.end()) >= 0 || obf.which(se.start()) >= 0 )
                    {
                        seI = seJ;
                        break;
                    }
                }
            }

            cf.append(newVerticesForSplitEdge_(seI, layerI));
        }

        //- add faces to cells
        cellsFromCell[nLayers-layerI].append(cf);
        cellsFromCell[nLayers-1-layerI].append(cf);
    }

    # ifdef DEBUGLayer
    Pout << "New cells from cell " << cellI << " are " << cellsFromCell << endl;
    //::exit(1);

    Pout << "1. Newly generated cells " << cellsFromCell << endl;

    //- check if all generated cells are topologically closed
    forAll(cellsFromCell, cI)
    {
        const DynList<DynList<label, 8>, 10>& cellFaces = cellsFromCell[cI];

        DynList<edge, 12> edges;
        DynList<label, 12> nAppearances;

        forAll(cellFaces, fI)
        {
            const DynList<label, 8>& f = cellFaces[fI];

            forAll(f, eI)
            {
                const edge e(f[eI], f.fcElement(eI));

                const label pos = edges.containsAtPosition(e);

                if( pos < 0 )
                {
                    edges.append(e);
                    nAppearances.append(1);
                }
                else
                {
                    ++nAppearances[pos];
                }
            }
        }

        forAll(nAppearances, eI)
            if( nAppearances[eI] != 2 )
            {
                Pout << "Prism cell " << cI << " edge " << edges[eI]
                    << " is present " << nAppearances[eI] << " times!" << endl;
                abort(FatalError);
            }
    }
    # endif
}

//- generate new cells from a hex at a feature edge
void refineBoundaryLayers::generateNewCellsEdgeHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 64>& cellsFromCell
)
{
    const faceListPMG& faces = mesh_.faces();
    const cell& c = mesh_.cells()[cellI];

    # ifdef DEBUGLayer
    Pout << "Generating new cells from edge cell " << cellI << endl;
    # endif

    const label startBoundary = mesh_.boundaries()[0].patchStart();
    const PtrList<writeProcessorPatch>& procBnd = mesh_.procBoundaries();

    //- find the number of layers for this cell
    FixedList<label, 2> layersInDirection(-1), dirFace;
    label currDir(0);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace_.size()) )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        if( nLayersAtBndFace_[bfI] < 2 )
            continue;

        layersInDirection[currDir] = nLayersAtBndFace_[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    # ifdef DEBUGLayer
    Pout << "Directions " << dirFace << endl;
    Pout << "Number of layers in direction " << layersInDirection << endl;
    # endif

    if( layersInDirection[0] < 0 || layersInDirection[1] < 0 )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::generateNewCellsEdgeHex("
            "const label, DynList<DynList<DynList<label, 4>, 6>, 64>&)"
        ) << "Cannot refine edge hex " << cellI << abort(FatalError);
    }

    //- set the number of layers
    cellsFromCell.setSize(layersInDirection[0] * layersInDirection[1]);

    //- find the shared edge between the boundary faces
    const edge commonEdge =
        help::sharedEdge(faces[c[dirFace[0]]], faces[c[dirFace[1]]]);
    const label donorFace = dirFace[0];

    # ifdef DEBUGLayer
    Pout << "Common edge " << commonEdge << endl;
    Pout << "Donor face " << donorFace << endl;
    Pout << "Donor face points " << faces[c[donorFace]] << endl;
    # endif

    //- find the face attached to the starting point of the edge and
    //- the face attached to the end point of the edge
    label faceAtStart(-1), faceAtEnd(-1);

    forAll(c, fI)
    {
        if(
            (faces[c[fI]].which(commonEdge.start()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
            faceAtStart = fI;

        if(
            (faces[c[fI]].which(commonEdge.end()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
            faceAtEnd = fI;
    }

    # ifdef DEBUGLayer
    Pout << "Face at start " << faces[c[faceAtStart]] << endl;
    Pout << "Face at end " << faces[c[faceAtEnd]] << endl;
    # endif

    //- check the orientation of cross-split faces
    //- checking face at starting point
    edge e = help::sharedEdge(faces[c[donorFace]], faces[c[faceAtStart]]);
    label pos = help::positionOfEdgeInFace(e, faces[c[faceAtStart]]);
    bool startFaceOrientation(true);
    if( e.start() != faces[c[faceAtStart]][pos] )
        startFaceOrientation = false;

    pos = mesh_.faceIsInProcPatch(c[faceAtStart]);
    bool neiProcFaceStart(false);
    if( (pos >= 0) && !procBnd[pos].owner() )
        neiProcFaceStart = true;

    # ifdef DEBUGLayer
    Pout << "Start face shared edge " << e << " is at position " << pos << endl;
    Pout << "Starting face orientation " << startFaceOrientation << endl;
    # endif

    //- checking edge at end point
    e = help::sharedEdge(faces[c[donorFace]], faces[c[faceAtEnd]]);
    pos = help::positionOfEdgeInFace(e, faces[c[faceAtEnd]]);
    bool endFaceOrientation(true);
    if( e.start() == faces[c[faceAtEnd]][pos] )
        endFaceOrientation = false;

    pos = mesh_.faceIsInProcPatch(c[faceAtEnd]);
    bool neiProcFaceEnd(false);
    if( (pos >= 0) && !procBnd[pos].owner() )
        neiProcFaceEnd = true;

    # ifdef DEBUGLayer
    Pout << "End face shared edge " << e << " is at position " << pos << endl;
    Pout << "End face orientation " << endFaceOrientation << endl;
    # endif

    //- fill up the matrix of points for this cell
    //- the matrix is used for generation of new cells
    FixedList<DynList<DynList<label> >, 2> facePoints;

    //- fill in the data for a face at the starting point
    if( startFaceOrientation )
    {
        DynList<DynList<label> >& fp = facePoints[0];
        fp.setSize(layersInDirection[0]+1);
        forAll(fp, i)
            fp[i].setSize(layersInDirection[1]+1);

        const label fLabel = c[faceAtStart];
        forAllRow(facesFromFace_, fLabel, fI)
        {
            const label nfI = facesFromFace_(fLabel, fI);

            label i = (fI % layersInDirection[0]);
            label j = (fI / layersInDirection[0]);

            if( neiProcFaceStart )
            {
                i = (fI / layersInDirection[1]);
                j = (fI % layersInDirection[1]);
            }

            # ifdef DEBUGLayer
            Pout << "1. i = " << i << " j = " << j << endl;
            Pout << "New face " << newFaces_[nfI] << endl;
            # endif

            facePoints[0][i][j] = newFaces_(nfI, 0);
            facePoints[0][i+1][j] = newFaces_(nfI, 1);
            facePoints[0][i+1][j+1] = newFaces_(nfI, 2);
            facePoints[0][i][j+1] = newFaces_(nfI, 3);
        }
    }
    else
    {
        DynList<DynList<label> >& fp = facePoints[0];
        fp.setSize(layersInDirection[0]+1);
        forAll(fp, i)
            fp[i].setSize(layersInDirection[1]+1);

        const label fLabel = c[faceAtStart];
        forAllRow(facesFromFace_, fLabel, fI)
        {
            const label nfI = facesFromFace_(fLabel, fI);

            label i = (fI / layersInDirection[1]);
            label j = (fI % layersInDirection[1]);

            if( neiProcFaceStart )
            {
                i = (fI % layersInDirection[0]);
                j = (fI / layersInDirection[0]);
            }

            # ifdef DEBUGLayer
            Pout << "2. i = " << i << " j = " << j << endl;
            Pout << "New face " << newFaces_[nfI] << endl;
            # endif

            facePoints[0][i][j] = newFaces_(nfI, 0);
            facePoints[0][i+1][j] = newFaces_(nfI, 3);
            facePoints[0][i+1][j+1] = newFaces_(nfI, 2);
            facePoints[0][i][j+1] = newFaces_(nfI, 1);
        }
    }

    //- fill in the data for a face at the end point
    if( endFaceOrientation )
    {
        DynList<DynList<label> >& fp = facePoints[1];
        fp.setSize(layersInDirection[0]+1);
        forAll(fp, i)
            fp[i].setSize(layersInDirection[1]+1);

        const label fLabel = c[faceAtEnd];
        forAllRow(facesFromFace_, fLabel, fI)
        {
            const label nfI = facesFromFace_(fLabel, fI);

            label i = (fI % layersInDirection[0]);
            label j = (fI / layersInDirection[0]);

            if( neiProcFaceEnd )
            {
                i = (fI / layersInDirection[1]);
                j = (fI % layersInDirection[1]);
            }

            # ifdef DEBUGLayer
            Pout << "3. i = " << i << " j = " << j << endl;
            Pout << "New face " << newFaces_[nfI] << endl;
            # endif

            facePoints[1][i][j] = newFaces_(nfI, 0);
            facePoints[1][i+1][j] = newFaces_(nfI, 1);
            facePoints[1][i+1][j+1] = newFaces_(nfI, 2);
            facePoints[1][i][j+1] = newFaces_(nfI, 3);
        }
    }
    else
    {
        DynList<DynList<label> >& fp = facePoints[1];
        fp.setSize(layersInDirection[0]+1);
        forAll(fp, i)
            fp[i].setSize(layersInDirection[1]+1);

        const label fLabel = c[faceAtEnd];
        forAllRow(facesFromFace_, fLabel, fI)
        {
            const label nfI = facesFromFace_(fLabel, fI);

            label i = (fI / layersInDirection[1]);
            label j = (fI % layersInDirection[1]);

            if( neiProcFaceEnd )
            {
                i = (fI % layersInDirection[0]);
                j = (fI / layersInDirection[0]);
            }

            # ifdef DEBUGLayer
            Pout << "4. i = " << i << " j = " << j << endl;
            Pout << "New face " << newFaces_[nfI] << endl;
            # endif

            facePoints[1][i][j] = newFaces_(nfI, 0);
            facePoints[1][i+1][j] = newFaces_(nfI, 3);
            facePoints[1][i+1][j+1] = newFaces_(nfI, 2);
            facePoints[1][i][j+1] = newFaces_(nfI, 1);
        }
    }

    # ifdef DEBUGLayer
    Pout << "Face points " << facePoints << endl;
    # endif

    //- generate new cells from the current cells
    label counter(0);
    for(label i=layersInDirection[0]-1;i>=0;--i)
    {
        for(label j=layersInDirection[1]-1;j>=0;--j)
        {
            DynList<DynList<label, 4>, 6>& cellFaces = cellsFromCell[counter];
            cellFaces.clear();

            if(
                (i == layersInDirection[0] - 1) &&
                (j == layersInDirection[1] - 1)
            )
            {
                # ifdef DEBUGLayer
                Pout << "1. Generating cell i = " << i << " j = " << j << endl;
                # endif

                //- this cell might not be a quad and
                //- it therefore requires special treatment
                forAll(c, fI)
                {
                    if( (fI == dirFace[0]) || (fI == dirFace[1]) )
                        continue;

                    const label cfI = facesFromFace_.sizeOfRow(c[fI]) - 1;

                    const label nfI = facesFromFace_(c[fI], cfI);

                    DynList<label, 4> cf;
                    cf.setSize(newFaces_.sizeOfRow(nfI));

                    forAllRow(newFaces_, nfI, pI)
                        cf[pI] = newFaces_(nfI, pI);

                    cellFaces.append(cf);
                }

                //- generate two missing faces
                DynList<label, 4> mf;
                mf.setSize(4);

                //- first missing face
                mf[0] = facePoints[0][i][j];
                mf[1] = facePoints[1][i][j];
                mf[2] = facePoints[1][i][j+1];
                mf[3] = facePoints[0][i][j+1];
                cellFaces.append(mf);

                //- second missing face
                mf[1] = facePoints[0][i+1][j];
                mf[2] = facePoints[1][i+1][j];
                mf[3] = facePoints[1][i][j];
                cellFaces.append(mf);
            }
            else
            {
                # ifdef DEBUGLayer
                Pout << "2. Generating cell i = " << i << " j = " << j << endl;
                # endif

                //- generate a hex cell from the matrix of cell points
                cellFaces.setSize(6);
                forAll(cellFaces, cfI)
                    cellFaces[cfI].setSize(4);

                //- first face
                cellFaces[0][0] = facePoints[0][i][j];
                cellFaces[0][1] = facePoints[0][i][j+1];
                cellFaces[0][2] = facePoints[0][i+1][j+1];
                cellFaces[0][3] = facePoints[0][i+1][j];

                //- second face
                cellFaces[1][0] = facePoints[1][i][j];
                cellFaces[1][1] = facePoints[1][i+1][j];
                cellFaces[1][2] = facePoints[1][i+1][j+1];
                cellFaces[1][3] = facePoints[1][i][j+1];

                //- third face
                cellFaces[2][0] = facePoints[0][i][j];
                cellFaces[2][1] = facePoints[0][i+1][j];
                cellFaces[2][2] = facePoints[1][i+1][j];
                cellFaces[2][3] = facePoints[1][i][j];

                //- fourth face
                cellFaces[3][0] = facePoints[0][i][j+1];
                cellFaces[3][1] = facePoints[1][i][j+1];
                cellFaces[3][2] = facePoints[1][i+1][j+1];
                cellFaces[3][3] = facePoints[0][i+1][j+1];

                //- fifth face
                cellFaces[4][0] = facePoints[0][i][j];
                cellFaces[4][1] = facePoints[1][i][j];
                cellFaces[4][2] = facePoints[1][i][j+1];
                cellFaces[4][3] = facePoints[0][i][j+1];

                //- sixth face
                cellFaces[5][0] = facePoints[0][i+1][j];
                cellFaces[5][1] = facePoints[0][i+1][j+1];
                cellFaces[5][2] = facePoints[1][i+1][j+1];
                cellFaces[5][3] = facePoints[1][i+1][j];
            }

            ++counter;
        }
    }

    # ifdef DEBUGLayer
    Pout << "Newly generated cells " << cellsFromCell << endl;

    //- check if all generated cells are topologically closed
    forAll(cellsFromCell, cI)
    {
        const DynList<DynList<label, 4>, 6>& cellFaces = cellsFromCell[cI];

        DynList<edge, 12> edges;
        DynList<label, 12> nAppearances;

        forAll(cellFaces, fI)
        {
            const DynList<label, 4>& f = cellFaces[fI];

            forAll(f, eI)
            {
                const edge e(f[eI], f.fcElement(eI));

                const label pos = edges.containsAtPosition(e);

                if( pos < 0 )
                {
                    edges.append(e);
                    nAppearances.append(1);
                }
                else
                {
                    ++nAppearances[pos];
                }
            }
        }

        forAll(nAppearances, eI)
            if( nAppearances[eI] != 2 )
            {
                Pout << "Edge hex cell " << cI << " edge " << edges[eI]
                    << " is present " << nAppearances[eI] << " times!" << endl;
                abort(FatalError);
            }
    }

    //::exit(1);
    # endif
}

//- generate new cells from a hex at a corner
void refineBoundaryLayers::generateNewCellsCornerHex
(
    const label cellI,
    DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
)
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();

    # ifdef DEBUGLayer
    Pout << "Generating new cells from cell " << cellI << endl;
    Pout << "Cell faces " << c << endl;
    # endif

    const label startBoundary = mesh_.boundaries()[0].patchStart();

    //- find the number of layers for this cell
    FixedList<label, 3> layersInDirection(-1), dirFace;
    FixedList<bool, 6> usedDirection(false);
    label currDir(0);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace_.size()) )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        if( nLayersAtBndFace_[bfI] < 2 )
            continue;

        usedDirection[fI] = true;
        layersInDirection[currDir] = nLayersAtBndFace_[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    //- find a common point for all three boundary faces
    FixedList<DynList<label, 4>, 3> bndFaces;
    forAll(dirFace, i)
    {
        bndFaces[i].setSize(4);
        forAll(faces[c[dirFace[i]]], pI)
            bndFaces[i][pI] = faces[c[dirFace[i]]][pI];
    }

    const label commonPoint = help::sharedVertex(bndFaces);

    # ifdef DEBUGLayer
    Pout << "Used directions " << usedDirection << endl;
    Pout << "Layers in direction " << layersInDirection << endl;
    Pout << "dirFace " << dirFace << endl;
    Pout << "Common point " << commonPoint << endl;

    forAll(dirFace, i)
        Pout << "bnd face " << i << " is " << faces[c[dirFace[i]]] << endl;
    # endif

    //- find the position of the common point in each boundary face
    const face& baseFace = faces[c[dirFace[0]]];
    const label posInBndFace = baseFace.which(commonPoint);

    FixedList<label, 3> splitEdgeDirection;

    forAllRow(splitEdgesAtPoint_, commonPoint, i)
    {
        const edge& se = splitEdges_[splitEdgesAtPoint_(commonPoint, i)];

        if( se == baseFace.faceEdge(posInBndFace) )
        {
            //- this edge is in j direction
            splitEdgeDirection[1] = splitEdgesAtPoint_(commonPoint, i);
        }
        else if( se == baseFace.faceEdge(baseFace.rcIndex(posInBndFace)) )
        {
            //- this edge is in i diretion
            splitEdgeDirection[0] = splitEdgesAtPoint_(commonPoint, i);
        }
        else if( splitEdgesAtPoint_.sizeOfRow(commonPoint) == 3 )
        {
            //- this point is in k direction
            splitEdgeDirection[2] = splitEdgesAtPoint_(commonPoint, i);
        }
        else
        {
            //- this situation is not allowed
            FatalErrorIn
            (
                "void refineBoundaryLayers::generateNewCellsCornerHex("
                "const label, DynList<DynList<DynList<label, 4>, 6>, 256>&)"
            ) << "Cannot refine layer for cell " << cellI << abort(FatalError);
        }
    }

    # ifdef DEBUGLayer
    forAll(splitEdgeDirection, i)
        Pout << "Split edge in direction " << i << " has nodes "
             << splitEdges_[splitEdgeDirection[i]]
             << " number of points on split edge "
             << newVerticesForSplitEdge_.sizeOfRow(splitEdgeDirection[i])
             << endl;
    # endif

    //- create a matrix of points forming new cells
    DynList<DynList<DynList<label> > > cellPoints;
    cellPoints.setSize
    (
        newVerticesForSplitEdge_.sizeOfRow(splitEdgeDirection[0])
    );

    forAll(cellPoints, i)
    {
        const label nJ =
            newVerticesForSplitEdge_.sizeOfRow(splitEdgeDirection[1]);

        cellPoints[i].setSize(nJ);

        forAll(cellPoints[i], j)
        {
            const label nK =
                newVerticesForSplitEdge_.sizeOfRow(splitEdgeDirection[2]);

            cellPoints[i][j].setSize(nK);
            cellPoints[i][j] = -1;
        }
    }

    //- find the direction od other boundary faces
    //- in the local coordinate system
    FixedList<label, 3> permutation;
    permutation[0] = 0;

    label helper = help::positionOfEdgeInFace
    (
        baseFace.faceEdge(baseFace.rcIndex(posInBndFace)),
        faces[c[dirFace[1]]]
    );

    if( helper >= 0 )
    {
        permutation[1] = 1;
        permutation[2] = 2;
    }
    else
    {
        permutation[1] = 2;
        permutation[2] = 1;
    }

    const edge seDir0 = splitEdges_[splitEdgeDirection[0]];
    const label nLayersI = layersInDirection[permutation[2]];
    const edge seDir1 = splitEdges_[splitEdgeDirection[1]];
    const label nLayersJ = layersInDirection[permutation[1]];
    const edge seDir2 = splitEdges_[splitEdgeDirection[2]];
    const label nLayersK = layersInDirection[permutation[0]];

    //- start filling in the data
    //- base face is located at k = 0 and has the opposite rotation
    forAllRow(facesFromFace_, c[dirFace[0]], fI)
    {
        const label nfI = facesFromFace_(c[dirFace[0]], fI);

        const label i = (fI / nLayersJ);
        const label j = (fI % nLayersJ);

        # ifdef DEBUGLayer
        Pout << "1. i = " << i << " j = " << j << endl;
        Pout << "New face " << newFaces_[nfI] << endl;
        # endif

        cellPoints[i][j][0] = newFaces_(nfI, 0);
        cellPoints[i+1][j][0] = newFaces_(nfI, 3);
        cellPoints[i+1][j+1][0] = newFaces_(nfI, 2);
        cellPoints[i][j+1][0] = newFaces_(nfI, 1);
    }

    //- permutation[1] is at j = 0
    forAllRow(facesFromFace_, c[dirFace[permutation[1]]], fI)
    {
        const label nfI = facesFromFace_(c[dirFace[permutation[1]]], fI);
        const label i = (fI % nLayersI);
        const label k = (fI / nLayersI);

        # ifdef DEBUGLayer
        Pout << "2. i = " << i << " k = " << k << endl;
        Pout << "New face " << newFaces_[nfI] << endl;
        # endif

        cellPoints[i][0][k] = newFaces_(nfI, 0);
        cellPoints[i+1][0][k] = newFaces_(nfI, 1);
        cellPoints[i+1][0][k+1] = newFaces_(nfI, 2);
        cellPoints[i][0][k+1] = newFaces_(nfI, 3);
    }

    //- permutation[2] is at i = 0
    forAllRow(facesFromFace_, c[dirFace[permutation[2]]], fI)
    {
        const label nfI = facesFromFace_(c[dirFace[permutation[2]]], fI);
        const label j = (fI / nLayersK);
        const label k = (fI % nLayersK);

        # ifdef DEBUGLayer
        Pout << "3. j = " << j << " k = " << k << endl;
        Pout << "New face " << newFaces_[nfI] << endl;
        # endif

        cellPoints[0][j][k] = newFaces_(nfI, 0);
        cellPoints[0][j][k+1] = newFaces_(nfI, 1);
        cellPoints[0][j+1][k+1] = newFaces_(nfI, 2);
        cellPoints[0][j+1][k] = newFaces_(nfI, 3);
    }

    //- find the face attached to the end of each split edge
    const PtrList<writeProcessorPatch>& procBnd = mesh_.procBoundaries();
    forAll(usedDirection, dirI)
    {
        if( usedDirection[dirI] )
            continue;

        label procPatch = mesh_.faceIsInProcPatch(c[dirI]);
        bool neiProcFace(false);
        if(
            procPatch >= 0 &&
            !procBnd[procPatch].owner()
        )
            neiProcFace = true;

        if( faces[c[dirI]].which(seDir0.end()) >= 0 )
        {
            const edge shEdge =
                help::sharedEdge
                (
                    faces[c[dirI]],
                    faces[c[dirFace[permutation[0]]]]
                );

            //- this face is attached to the end of the edge in direction 0
            //- it has its i coordinate equal to nLayersI
            if( shEdge.end() == seDir0.end() )
            {
                //- desired orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label j = (fI / nLayersK);
                    label k = (fI % nLayersK);

                    if( neiProcFace )
                    {
                        j = (fI % nLayersJ);
                        k = (fI / nLayersJ);
                    }

                    # ifdef DEBUGLayer
                    Pout << "4. j = " << j << " k = " << k << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[nLayersI][j][k] = newFaces_(nfI, 0);
                    cellPoints[nLayersI][j][k+1] = newFaces_(nfI, 1);
                    cellPoints[nLayersI][j+1][k+1] = newFaces_(nfI, 2);
                    cellPoints[nLayersI][j+1][k] = newFaces_(nfI, 3);
                }
            }
            else
            {
                //- opposite orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label j = (fI % nLayersJ);
                    label k = (fI / nLayersJ);

                    if( neiProcFace )
                    {
                        j = (fI / nLayersK);
                        k = (fI % nLayersK);
                    }

                    # ifdef DEBUGLayer
                    Pout << "5. j = " << j << " k = " << k << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[nLayersI][j][k] = newFaces_(nfI, 0);
                    cellPoints[nLayersI][j][k+1] = newFaces_(nfI, 3);
                    cellPoints[nLayersI][j+1][k+1] = newFaces_(nfI, 2);
                    cellPoints[nLayersI][j+1][k] = newFaces_(nfI, 1);
                }
            }
        }
        else if( faces[c[dirI]].which(seDir1.end()) >= 0 )
        {
            const edge shEdge =
                help::sharedEdge
                (
                    faces[c[dirI]],
                    faces[c[dirFace[permutation[0]]]]
                );

            if( shEdge.end() == seDir1.end() )
            {
                //- desired orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label i = (fI / nLayersK);
                    label k = (fI % nLayersK);

                    if( neiProcFace )
                    {
                        i = (fI % nLayersI);
                        k = (fI / nLayersI);
                    }
                    # ifdef DEBUGLayer
                    Pout << "6. i = " << i << " k = " << k << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[i][nLayersJ][k] = newFaces_(nfI, 0);
                    cellPoints[i][nLayersJ][k+1] = newFaces_(nfI, 1);
                    cellPoints[i+1][nLayersJ][k+1] = newFaces_(nfI, 2);
                    cellPoints[i+1][nLayersJ][k] = newFaces_(nfI, 3);
                }
            }
            else
            {
                //- opposite orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label i = (fI % nLayersI);
                    label k = (fI / nLayersI);

                    if( neiProcFace )
                    {
                        i = (fI / nLayersK);
                        k = (fI % nLayersK);
                    }

                    # ifdef DEBUGLayer
                    Pout << "7. i = " << i << " k = " << k << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[i][nLayersJ][k] = newFaces_(nfI, 0);
                    cellPoints[i+1][nLayersJ][k] = newFaces_(nfI, 1);
                    cellPoints[i+1][nLayersJ][k+1] = newFaces_(nfI, 2);
                    cellPoints[i][nLayersJ][k+1] = newFaces_(nfI, 3);
                }
            }
        }
        else if( faces[c[dirI]].which(seDir2.end()) >= 0 )
        {
            const edge shEdge =
                help::sharedEdge
                (
                    faces[c[dirI]],
                    faces[c[dirFace[permutation[2]]]]
                );

            if( shEdge.end() == seDir2.end() )
            {
                //- desired orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label i = (fI % nLayersI);
                    label j = (fI / nLayersI);

                    if( neiProcFace )
                    {
                        i = (fI / nLayersJ);
                        j = (fI % nLayersJ);
                    }

                    # ifdef DEBUGLayer
                    Pout << "8. i = " << i << " j = " << j << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[i][j][nLayersK] = newFaces_(nfI, 0);
                    cellPoints[i+1][j][nLayersK] = newFaces_(nfI, 1);
                    cellPoints[i+1][j+1][nLayersK] = newFaces_(nfI, 2);
                    cellPoints[i][j+1][nLayersK] = newFaces_(nfI, 3);
                }
            }
            else
            {
                //- opposite orientation
                forAllRow(facesFromFace_, c[dirI], fI)
                {
                    const label nfI = facesFromFace_(c[dirI], fI);

                    label i = (fI / nLayersJ);
                    label j = (fI % nLayersJ);

                    if( neiProcFace )
                    {
                        i = (fI % nLayersI);
                        j = (fI / nLayersI);
                    }

                    # ifdef DEBUGLayer
                    Pout << "9. i = " << i << " j = " << j << endl;
                    Pout << "New face " << newFaces_[nfI] << endl;
                    # endif

                    cellPoints[i][j][nLayersK] = newFaces_(nfI, 0);
                    cellPoints[i][j+1][nLayersK] = newFaces_(nfI, 1);
                    cellPoints[i+1][j+1][nLayersK] = newFaces_(nfI, 2);
                    cellPoints[i+1][j][nLayersK] = newFaces_(nfI, 3);
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void refineBoundaryLayers::generateNewCellsCornerHex("
                "const label, DynList<DynList<DynList<label, 4>, 6>, 256>&)"
            ) << "Corner cell " << cellI << " may not be a hex"
              << abort(FatalError);
        }
    }

    # ifdef DEBUGLayer
    Pout << "cellPoints before inner points " << cellPoints << endl;
    # endif

    //- generate points inside the cell via transfinite interpolation
    const pointFieldPMG& points = mesh_.points();
    forAll(cellPoints, i)
    {
        forAll(cellPoints[i], j)
        {
            forAll(cellPoints[i][j], k)
            {
                if( cellPoints[i][j][k] < 0 )
                {
                    //- generate a vertex via transfinite interpolation of
                    //- vertices at faces
                    const scalar u
                    (
                        Foam::mag
                        (
                            points[cellPoints[i][0][0]] -
                            points[cellPoints[0][0][0]]
                        ) /
                        seDir0.mag(points)
                    );

                    const scalar v
                    (
                        Foam::mag
                        (
                            points[cellPoints[0][j][0]] -
                            points[cellPoints[0][0][0]]
                        ) /
                        seDir1.mag(points)
                    );

                    const scalar w
                    (
                        Foam::mag
                        (
                            points[cellPoints[0][0][k]] -
                            points[cellPoints[0][0][0]]
                        ) /
                        seDir2.mag(points)
                    );

                    //- corner points
                    const point& v000 = points[cellPoints[0][0][0]];
                    const point& v100 = points[cellPoints[nLayersI][0][0]];
                    const point& v110 =
                        points[cellPoints[nLayersI][nLayersJ][0]];
                    const point& v010 = points[cellPoints[0][nLayersJ][0]];
                    const point& v001 = points[cellPoints[0][0][nLayersK]];
                    const point& v101 =
                        points[cellPoints[nLayersI][0][nLayersK]];
                    const point& v111 =
                        points[cellPoints[nLayersI][nLayersJ][nLayersK]];
                    const point& v011 =
                        points[cellPoints[0][nLayersJ][nLayersK]];

                    //- calculate coordinates of the new vertex
                    const point newP =
                        (1.0 - u) * (1.0 - v) * (1.0 - w) * v000 +
                        u * (1.0 - v) * (1.0 - w) * v100 +
                        u * v * (1.0 - w) * v110 +
                        (1.0 - u) * v * (1.0 - w) * v010 +
                        (1.0 - u) * (1.0 - v) * w * v001 +
                        u * (1.0 - v) * w * v101 +
                        u * v * w * v111 +
                        (1.0 - u) * v * w * v011;

                    //- add the point to the mesh
                    cellPoints[i][j][k] = points.size();
                    mesh_.appendVertex(newP);
                }
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "cellPoints " << cellPoints << endl;
    # endif

    //- create new cells from
    cellsFromCell.setSize(nLayersI * nLayersJ * nLayersK);

    # ifdef DEBUGLayer
    Pout << "Number of new cells " << cellsFromCell.size() << endl;
    # endif

    label counter(0);
    for(label i=nLayersI-1;i>=0;--i)
    {
        for(label j=nLayersJ-1;j>=0;--j)
        {
            for(label k=nLayersK-1;k>=0;--k)
            {
                DynList<DynList<label, 4>, 6>& cellFaces =
                    cellsFromCell[counter++];

                if( (i == nLayersI-1) && (j == nLayersJ-1) && (k == nLayersK-1) )
                {
                    cellFaces.clear();

                    //- this cell may not be a hex
                    # ifdef DEBUGLayer
                    Pout << "1. Generating cell i = " << i
                         << " j = " << j
                         << " k = " << k << endl;
                    # endif

                    //- add the faces originating from internal faces
                    //- these faces may not be quads
                    forAll(c, fI)
                    {
                        if( usedDirection[fI] )
                            continue;

                        const label cfI = facesFromFace_.sizeOfRow(c[fI]) - 1;

                        const label nfI = facesFromFace_(c[fI], cfI);

                        DynList<label, 4> cf;
                        cf.setSize(newFaces_.sizeOfRow(nfI));

                        forAllRow(newFaces_, nfI, pI)
                            cf[pI] = newFaces_(nfI, pI);

                        cellFaces.append(cf);
                    }

                    //- generate two missing faces
                    DynList<label, 4> mf;
                    mf.setSize(4);

                    //- first missing face
                    mf[0] = cellPoints[i][j][k];
                    mf[1] = cellPoints[i][j][k+1];
                    mf[2] = cellPoints[i][j+1][k+1];
                    mf[3] = cellPoints[i][j+1][k];
                    cellFaces.append(mf);

                    //- second missing face
                    mf[1] = cellPoints[i+1][j][k];
                    mf[2] = cellPoints[i+1][j][k+1];
                    mf[3] = cellPoints[i][j][k+1];
                    cellFaces.append(mf);

                    //- third missing face
                    mf[1] = cellPoints[i][j+1][k];
                    mf[2] = cellPoints[i+1][j+1][k];
                    mf[3] = cellPoints[i+1][j][k];
                    cellFaces.append(mf);
                }
                else
                {
                    //- generate a hex cell
                    # ifdef DEBUGLayer
                    Pout << "2. Generating cell i = " << i
                         << " j = " << j << " k = " << k << endl;
                    # endif

                    //- generate a hex cell from the matrix of cell points
                    cellFaces.setSize(6);
                    forAll(cellFaces, cfI)
                        cellFaces[cfI].setSize(4);

                    //- first face
                    cellFaces[0][0] = cellPoints[i][j][k];
                    cellFaces[0][1] = cellPoints[i][j][k+1];
                    cellFaces[0][2] = cellPoints[i][j+1][k+1];
                    cellFaces[0][3] = cellPoints[i][j+1][k];

                    //- second face
                    cellFaces[1][0] = cellPoints[i+1][j][k];
                    cellFaces[1][1] = cellPoints[i+1][j+1][k];
                    cellFaces[1][2] = cellPoints[i+1][j+1][k+1];
                    cellFaces[1][3] = cellPoints[i+1][j][k+1];

                    //- third face
                    cellFaces[2][0] = cellPoints[i][j][k];
                    cellFaces[2][1] = cellPoints[i][j+1][k];
                    cellFaces[2][2] = cellPoints[i+1][j+1][k];
                    cellFaces[2][3] = cellPoints[i+1][j][k];

                    //- fourth face
                    cellFaces[3][0] = cellPoints[i][j][k+1];
                    cellFaces[3][1] = cellPoints[i+1][j][k+1];
                    cellFaces[3][2] = cellPoints[i+1][j+1][k+1];
                    cellFaces[3][3] = cellPoints[i][j+1][k+1];

                    //- fifth face
                    cellFaces[4][0] = cellPoints[i][j][k];
                    cellFaces[4][1] = cellPoints[i+1][j][k];
                    cellFaces[4][2] = cellPoints[i+1][j][k+1];
                    cellFaces[4][3] = cellPoints[i][j][k+1];

                    //- sixth face
                    cellFaces[5][0] = cellPoints[i][j+1][k];
                    cellFaces[5][1] = cellPoints[i][j+1][k+1];
                    cellFaces[5][2] = cellPoints[i+1][j+1][k+1];
                    cellFaces[5][3] = cellPoints[i+1][j+1][k];
                }
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "New cells from corner cell " << cellI
         << " are " << cellsFromCell << endl;
    //::exit(1);
    # endif
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

        if( nLayersAtBndFace_[bfI] > 1 )
            ++refType[cellI];
    }

    label nNewCells(0);
    forAll(nCellsFromCell, cellI)
        nNewCells += (nCellsFromCell[cellI] - 1);

    # ifdef DEBUGLayer
    forAll(nCellsFromCell, cellI)
    {
        Info << "\nCell " << cellI << endl;
        Info << "nCellsFromCell " << nCellsFromCell[cellI] << endl;
        Info << "Ref type " << refType[cellI] << endl;
    }
    #  endif

    const label totalNumNewCells = returnReduce(nNewCells, sumOp<label>());
    Info << "Number of newly generated cells " << totalNumNewCells << endl;

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
            //- this cell is not refined
            //- update face labels
            newCellsFromCell.append(cellI, cellI);

            cell& c = cells[cellI];

            forAll(c, fI)
                c[fI] = facesFromFace_(c[fI], 0);
        }
        else if( refType[cellI] == 1 )
        {
            //- generate new cells from this prism refined in one direction
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

                # ifdef DEBUGLayer
                Info << "Adding cell " << (cI==0?cellI:nCells)
                     << " originating from cell " << cellI << endl;
                # endif

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

    //- update cell sets
    mesh_.updateCellSubsets(newCellsFromCell);
    newCellsFromCell.setSize(0);

    //- point-faces addressing is not needed any more
    pointNewFaces.setSize(0);

    //- copy the newFaces to the mesh
    const label nOrigInternalFaces = mesh_.nInternalFaces();

    # ifdef DEBUGLayer
    Info << "Copying internal faces " << endl;
    Info << "Original number of internal faces " << nOrigInternalFaces << endl;
    # endif

    //- store internal faces originating from existing faces
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
    # ifdef DEBUGLayer
    Info << "Copying newly generated internal faces" << endl;
    Info << "nNewInternalFaces " << currFace << endl;
    Info << "numFacesBefore " << numFacesBefore << endl;
    Info << "Total number of faces " << newFaces_.size() << endl;
    # endif

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
    # ifdef DEBUGLayer
    Info << "Copying boundary faces " << endl;
    Info << "currFace " << currFace << endl;
    Info << "Faces size " << faces.size() << endl;
    Info << "Initial number of faces " << facesFromFace_.size() << endl;
    # endif

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
        # ifdef DEBUGLayer
        Info << "Copying processor faces" << endl;
        # endif

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
            procBoundaries[patchI].patchStart() = newStart;
            procBoundaries[patchI].patchSize() = nNewFacesInPatch;
        }
    }

    # ifdef DEBUGLayer
    Info << "Faces after refinement " << faces << endl;
    Info << "newFaceLabel " << newFaceLabel << endl;
    # endif

    //- update face subsets
    Info << "Updating subsets" << endl;
    mesh_.updateFaceSubsets(facesFromFace_);
    facesFromFace_.setSize(0);
    newFaces_.setSize(0);

    //- update cells to match the faces
    # ifdef DEBUGLayer
    Info << "Updating cells to match new faces" << endl;
    # endif

    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            c[fI] = newFaceLabel[c[fI]];
    }

    # ifdef DEBUGLayer
    Info << "Cleaning mesh " << endl;
    # endif
    //- delete all adressing which is no longer up-to-date
    meshModifier.clearAll();

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<edge> edges;
                DynList<label> nAppearances;
                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, eI)
                    {
                        const edge e = f.faceEdge(eI);

                        const label pos = edges.containsAtPosition(e);

                        if( pos < 0 )
                        {
                            edges.append(e);
                            nAppearances.append(1);
                        }
                        else
                        {
                            ++nAppearances[pos];
                        }
                    }
                }

                bool badCell(false);
                forAll(nAppearances, i)
                    if( nAppearances[i] != 2 )
                    {
                        badCell = true;
                        break;

                    }

                if( badCell )
                {
                    Pout << "Cell " << cellI
                         << " is not topologically closed" << endl;

                    forAll(c, fI)
                        Pout << "Face " << c[fI] << " with points "
                             << faces[c[fI]] << endl;
                    Pout << "Cell edges " << edges << endl;
                    Pout << "nAppearances " << nAppearances << endl;
                }
            }
        }

        returnReduce(1, sumOp<label>());
    }

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(faces, faceI)
            {
                Pout << "Face " << faceI << " owner " << owner[faceI]
                     << " neighbour " << neighbour[faceI]
                     << " face points " << faces[faceI] << endl;
            }

            forAll(cells, cellI)
                Pout << "Cell " << cellI << " has faces "
                     << cells[cellI] << endl;
        }

        returnReduce(procI, maxOp<label>());
    }
    # endif

    Info << "Finished generating new cells " << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

