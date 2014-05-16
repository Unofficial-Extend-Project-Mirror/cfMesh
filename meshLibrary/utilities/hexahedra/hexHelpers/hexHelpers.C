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

#include "hexHelpers.H"
#include "helperFunctions.H"
#include "FRWGraph.H"
#include "VRWGraphSMPModifier.H"
#include "hexMatcher.H"
#include "polyMeshGenAddressing.H"
#include "labelPair.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace hexHelpers
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool allHexMesh(const polyMeshGen& mesh)
{
    const labelList& owner = mesh.owner();
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();

    hexMatcher matcher;

    bool allHex(true);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        if( !matcher.matchShape(true, faces, owner, cellI, c) )
            allHex = false;

        if( !allHex )
            cellI = cells.size();
    }

    returnReduce(allHex, maxOp<bool>());

    return allHex;
}

label findColumnCells
(
    const polyMeshGen& mesh,
    const label faceI,
    labelLongList& columnCells
)
{
    columnCells.setSize(0);

    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    boolList addedCell(mesh.cells().size(), false);

    labelLongList front;
    front.append(faceI);

    while( front.size() )
    {
        const label fLabel = front.removeLastElement();

        if( !addedCell[owner[fLabel]] )
        {
            const cell& own = cells[owner[fLabel]];
            const label ownF = own.opposingFaceLabel(fLabel, faces);

            addedCell[owner[fLabel]] = true;
            columnCells.append(owner[fLabel]);
            front.append(ownF);
        }
        if( (neighbour[fLabel] >= 0) && !addedCell[neighbour[fLabel]] )
        {
            const cell& nei = cells[neighbour[fLabel]];
            const label neiF = nei.opposingFaceLabel(fLabel, faces);

            addedCell[neighbour[fLabel]] = true;
            columnCells.append(neighbour[fLabel]);
            front.append(neiF);
        }
    }

    return returnReduce(columnCells.size(), sumOp<label>());
}

label findColumnCells
(
    polyMeshGen& mesh,
    const label faceI,
    const word& columnCellSet
)
{
    labelLongList columnCells;
    findColumnCells(mesh, faceI, columnCells);

    const label cID = mesh.addCellSubset(columnCellSet);
    forAll(columnCells, i)
        mesh.addCellToSubset(cID, columnCells[i]);

    return columnCells.size();
}

bool selfIntersectingColumn(const polyMeshGen& mesh, const boolList& columnCells)
{
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();

    bool selfIntersecting(false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(cells, cellI)
    {
        if( !columnCells[cellI] || selfIntersecting )
            continue;

        const cell& c = cells[cellI];

        DynList<label> columnFaces;

        forAll(c, fI)
        {
            label nei = neighbour[c[fI]];
            if( nei < 0 )
                continue;
            if( nei == cellI )
                nei = owner[c[fI]];

            if( columnCells[nei] )
                columnFaces.append(c[fI]);
        }

        if( columnFaces.size() > 1 )
        {
            forAllReverse(columnFaces, i)
            {
                const face& f = faces[columnFaces[i]];

                for(label j=i-1;j>=0;j--)
                    if( help::shareAnEdge(f, faces[columnFaces[j]]) )
                        selfIntersecting = true;
            }
        }
    }

    return selfIntersecting;
}

bool selfIntersectingColumn(const polyMeshGen& mesh, const labelLongList& columnCells)
{
    boolList cellsInColumn(mesh.cells().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(columnCells, i)
        cellsInColumn[columnCells[i]] = true;

    return selfIntersectingColumn(mesh, cellsInColumn);
}

bool selfIntersectingColumn(const polyMeshGen& mesh, const word& columnCellSet)
{
    labelLongList columnCells;

    const label setID = mesh.cellSubsetIndex(columnCellSet);
    mesh.cellsInSubset(setID, columnCells);

    return selfIntersectingColumn(mesh, columnCells);
}

//- find all columns in hex mesh
label findAllColumns(polyMeshGen& mesh, const word& cellSetPrefix)
{
    FatalErrorIn
    (
        "label findAllColumns(polyMeshGen&, const word&)"
    ) << "Not implemented!" << exit(FatalError);

    FRWGraph<label, 3> cellInSheets(mesh.cells().size(), -1);

    label columnI(0);

    return columnI;
}

//- find all cells in a sheet containing a given face
label findSheetCells
(
    const polyMeshGen& mesh,
    const edge& sheetEdge,
    labelLongList& cellsInSheet
)
{
    cellsInSheet.setSize(0);

    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();

    boolList addedCell(cells.size(), false);

    const edgeList& edges = mesh.addressingData().edges();
    const VRWGraph& pointEdges = mesh.addressingData().pointEdges();
    const VRWGraph& edgeCells = mesh.addressingData().edgeCells();
    const VRWGraph& cellEdges = mesh.addressingData().cellEdges();

    labelLongList front;

    forAllRow(pointEdges, sheetEdge.start(), peI)
    {
        if( edges[pointEdges(sheetEdge.start(), peI)] == sheetEdge )
        {
            front.append(pointEdges(sheetEdge.start(), peI));
            break;
        }
    }

    while( front.size() )
    {
        const label eLabel = front.removeLastElement();
        const edge& e = edges[eLabel];

        forAllRow(edgeCells, eLabel, ecI)
        {
            const label cellI = edgeCells(eLabel, ecI);

            if( addedCell[cellI] )
                continue;

            cellsInSheet.append(cellI);
            addedCell[cellI] = true;

            const cell& c = cells[cellI];

            DynList<label, 2> sheetBndFaces;
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                if( (f.which(e.start()) >= 0) ^ (f.which(e.end()) >= 0) )
                    sheetBndFaces.append(c[fI]);
            }

            if( sheetBndFaces.size() != 2 )
            {
                Warning << "It seems that cell " << cellI
                    << " is not a hex!" << endl;
                return -1;
            }

            const face& sbf1 = faces[sheetBndFaces[0]];
            const face& sbf2 = faces[sheetBndFaces[1]];

            forAllRow(cellEdges, cellI, ceI)
            {
                const label cEdge = cellEdges(cellI, ceI);

                if( cEdge == eLabel )
                    continue;

                const edge& ce = edges[cEdge];

                if(
                    ((sbf1.which(ce[0]) >= 0) && (sbf2.which(ce[1]) >= 0)) ||
                    ((sbf1.which(ce[1]) >= 0) && (sbf2.which(ce[0]) >= 0))
                )
                {
                    front.append(cEdge);
                }
            }
        }
    }

    return cellsInSheet.size();
}

label findSheetCells
(
    polyMeshGen& mesh,
    const edge& sheetEdge,
    const word& sheetCellSet
)
{
    labelLongList sheetCells;
    findSheetCells(mesh, sheetEdge, sheetCells);

    const label cID = mesh.addCellSubset(sheetCellSet);
    forAll(sheetCells, i)
        mesh.addCellToSubset(cID, sheetCells[i]);

    return sheetCells.size();
}

//- check if any of the cells in the sheet share two or more faces
//- with another cell in the sheet
bool hasSheetDigons(const polyMeshGen& mesh, const boolList& sheetCells)
{
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const cellListPMG& cells = mesh.cells();

    bool hasDigons(false);

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # pragma omp parallel for num_threads(nThreads) schedule(static)
    # endif
    forAll(sheetCells, cellI)
    {
        if( !sheetCells[cellI] )
            continue;

        const cell& c = cells[cellI];

        labelHashSet neis;
        forAll(c, fI)
        {
            label nei = owner[c[fI]];
            if( nei == cellI )
                nei = neighbour[c[fI]];

            if( nei < 0 )
                continue;

            if( neis.found(nei) )
            {
                hasDigons = true;
                break;
            }

            neis.insert(nei);
        }

        //- break the loop in case any digons are found
        if( hasDigons )
            cellI = sheetCells.size();
    }

    return hasDigons;
}

bool hasSheetDigons
(
    const polyMeshGen& mesh,
    const labelLongList& sheetCells
)
{
    boolList sCells(mesh.cells().size());

    # ifdef USE_OMP
    # pragma omp parallel if( mesh.cells().size() > 1000 )
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(sCells, i)
            sCells[i] = false;

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for schedule(static, 1)
        # endif
        forAll(sheetCells, i)
            sCells[sheetCells[i]] = true;
    }

    return hasSheetDigons(mesh, sCells);
}

bool hasSheetDigons(const polyMeshGen& mesh, const word& sheetCellSet)
{
    const label sheetID = mesh.cellSubsetIndex(sheetCellSet);

    labelLongList sheetCells;
    mesh.cellsInSubset(sheetID, sheetCells);

    return hasSheetDigons(mesh, sheetCells);
}

//- find all sheets in hex mesh
void findAllSheets(polyMeshGen& mesh, const word& cellSetPrefix)
{
    FatalErrorIn
    (
        "void findAllSheets(polyMeshGen& mesh, const word& cellSetPrefix)"
    ) << "Not implemented!" << exit(FatalError);
}

// Modification tools
//- find all cells creating a single column and colapse the column
//- the column is collapsed by merging the point at the given position
//- in the selected face with the opposite point
bool collapseColumn
(
    polyMeshGen& mesh,
    const label faceI,
    const label positionInFace
)
{
    labelLongList columnCells;
    findColumnCells(mesh, faceI, columnCells);

    const label pointI = mesh.faces()[faceI][positionInFace];

    return collapseColumn(mesh, columnCells, pointI);
}

//- collapse the column by merging the vertices in the layer of edges
//- containing the given point with the opposite layer of edges
bool collapseColumn
(
    polyMeshGen& mesh,
    const labelLongList& columnCells,
    const label pointI
)
{
    polyMeshGenModifier meshModifier(mesh);

    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    //- check whether the collumn can be collapsed
    boolList removeCell(cells.size(), false), removeFace(faces.size(), false);

    forAll(columnCells, i)
        removeCell[columnCells[i]] = true;

    //- it is not possible to collapse self-intersecting columns
    //- check the number of neighbour cells in the column
    if( selfIntersectingColumn(mesh, removeCell) )
        return false;

    //- find edges of the column
    LongList<edge> edges;
    VRWGraph pEdges;
    pEdges.setSize(points.size());

    forAll(cells, cellI)
    {
        if ( !removeCell[cellI] )
            continue;

        const cell& c = cells[cellI];

        //- mark faces which are part of this column
        label faceI(-1);

        forAll(c, fI)
        {
            const label cfI = c[fI];
            label nei = neighbour[cfI];
            if( nei < 0 )
                continue;
            if( nei == cellI )
                nei = owner[cfI];

            if( removeCell[nei] )
            {
                faceI = cfI;
                break;
            }
        }

        if( faceI == -1 )
            continue;

        //- create 4 edges connecting vertices of opposite faces
        const label oppositeFace =
            c.opposingFaceLabel(faceI, faces);
        const face& f = faces[faceI];
        const face& of = faces[oppositeFace];

        //- mark faces for removal
        removeFace[faceI] = true;
        removeFace[oppositeFace] = true;

        forAll(c, fI)
        {
            if( (c[fI] == faceI) || (c[fI] == oppositeFace) )
                continue;

            const face& cf = faces[c[fI]];
            forAll(cf, eI)
            {
                const edge e = cf.faceEdge(eI);

                bool edgeExists(false);
                forAllRow(pEdges, e.start(), peI)
                    if( edges[pEdges(e.start(), peI)] == e )
                    {
                        edgeExists = true;
                        break;
                    }

                if( edgeExists )
                    continue;

                if( (
                        (f.which(e.start()) != -1) &&
                        (of.which(e.end()) != -1)
                    ) ||
                    (
                        (f.which(e.end()) != -1) &&
                        (of.which(e.start()) != -1)
                    )
                )
                {
                    pEdges.append(e.start(), edges.size());
                    pEdges.append(e.end(), edges.size());
                    edges.append(e);
                }
            }
        }
    }

    //- find all points connected to the selected point
    boolList removePoint(mesh.points().size(), false);

    labelLongList front;
    front.append(pointI);
    while( front.size() )
    {
        const label pLabel = front.removeLastElement();

        removePoint[pLabel] = true;

        forAllRow(pEdges, pLabel, peI)
        {
            const label nei = edges[pEdges(pLabel, peI)].otherVertex(pLabel);

            if( removePoint[nei] )
                continue;

            front.append(nei);
        }
    }

    //- pass through the faces marked for removal
    //- merge the opposite vertices
    //- update nodes of faces connected to merged nodes
    VRWGraph pFaces;
    pFaces.reverseAddressing(points.size(), faces);

    forAll(removeFace, faceI)
    {
        if( !removeFace[faceI] )
            continue;

        const face& f = faces[faceI];
        label pos(-1);
        forAll(f, pI)
        {
            if( removePoint[f[pI]] )
            {
                pos = pI;
                break;
            }
        }

        if( pos < 0 )
            FatalErrorIn
            (
                "bool collapseColumn(polyMeshGen&,"
                "const labelLongList&, const label"
            ) << "Cannot find position in face " << faceI << exit(FatalError);

        const label pPos = f[pos];
        const label pNextNext = f[(pos+2) % 4];
        const point newP = 0.5 * (points[pPos] + points[pNextNext]);

        points[pNextNext] = newP;

        forAllRow(pFaces, pPos, pfI)
        {
            face& cf = faces[pFaces(pPos, pfI)];

            const label pJ = cf.which(pPos);
            cf[pJ] = pNextNext;
        }
    }

    //- update neighbour cells of the collapsed cells
    forAll(cells, cellI)
    {
        if( !removeCell[cellI] )
            continue;

        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            if( removeFace[c[fI]] )
                continue;

            const face& cf = faces[c[fI]];

            for(label fJ=fI+1;fJ<c.size();++fJ)
            {
                if( removeFace[c[fJ]] )
                    continue;

                if( cf == faces[c[fJ]] )
                {
                    //- found duplicate face
                    label neiF = neighbour[c[fI]];
                    if( neiF == cellI )
                        neiF = owner[c[fI]];

                    label neiS = neighbour[c[fJ]];
                    if( neiS == cellI )
                        neiS = owner[c[fJ]];

                    if( (neiF == owner[c[fI]]) && (neiS > neiF) )
                    {
                        label pos(-1);
                        forAll(cells[neiS], i)
                        {
                            if( cells[neiS][i] == c[fJ] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        removeFace[c[fJ]] = true;
                        cells[neiS][pos] = c[fI];
                    }
                    else if( (neiS == owner[c[fJ]]) && (neiF > neiS) )
                    {
                        label pos(-1);
                        forAll(cells[neiF], i)
                        {
                            if( cells[neiF][i] == c[fI] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        removeFace[c[fI]] = true;
                        cells[neiF][pos] = c[fJ];
                    }
                    else if( (neiF == neighbour[c[fI]]) && (neiS > neiF) )
                    {
                        label pos(-1);
                        forAll(cells[neiS], i)
                        {
                            if( cells[neiS][i] == c[fJ] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        removeFace[c[fJ]] = true;
                        faces[c[fI]] = cf.reverseFace();
                        cells[neiS][pos] = c[fI];
                    }
                    else if( (neiS == neighbour[c[fI]]) && (neiF > neiS) )
                    {
                        label pos(-1);
                        forAll(cells[neiF], i)
                        {
                            if( cells[neiF][i] == c[fJ] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        removeFace[c[fI]] = true;
                        faces[c[fJ]] = faces[c[fJ]].reverseFace();
                        cells[neiF][pos] = c[fJ];
                    }
                    else if( (neiF == -1) && (neiS != -1) )
                    {
                        removeFace[c[fJ]] = true;
                        label pos(-1);

                        forAll(cells[neiS], i)
                        {
                            if( cells[neiS][i] == c[fJ] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        cells[neiS][pos] = c[fI];
                    }
                    else if( (neiS == -1) && (neiF != -1) )
                    {
                        removeFace[c[fI]] = true;
                        label pos(-1);

                        forAll(cells[neiF], i)
                        {
                            if( cells[neiF][i] == c[fI] )
                            {
                                pos = i;
                                break;
                            }
                        }

                        cells[neiF][pos] = c[fJ];
                    }
                    else
                    {
                        removeFace[c[fI]] = true;
                        removeFace[c[fJ]] = true;
                    }
               }
            }
        }

        //- set the cell size to zero
        cells[cellI].clear();
    }

    //- remove cells from the mesh and the unused vertices
    meshModifier.removeFaces(removeFace);

    label nCells(0);
    forAll(cells, cellI)
    {
        if( removeCell[cellI] )
            continue;

        if( cellI > nCells )
            cells[nCells].transfer(cells[cellI]);

        ++nCells;
    }
    cells.setSize(nCells);

    meshModifier.removeUnusedVertices();

    return true;
}

bool collapseColumn
(
    polyMeshGen& mesh,
    const word& columnCellSet,
    const label pointI
)
{
    labelLongList columnCells;
    const label cID = mesh.cellSubsetIndex(columnCellSet);
    mesh.cellsInSubset(cID, columnCells);

    return collapseColumn(mesh, columnCells, pointI);
}

//- find all cells creating a sheet containing a given face and extract
//- the sheet into a set of faces
bool extractSheet(polyMeshGen& mesh, const edge& sheetEdge)
{
    labelLongList sheetCells;
    findSheetCells(mesh, sheetEdge, sheetCells);

    return extractSheet(mesh, sheetCells);
}

//- extract the sheet of cells
bool extractSheet(polyMeshGen& mesh, const boolList& sheetCells)
{
    if( hasSheetDigons(mesh, sheetCells) )
    {
        Warning << "Sheet cannot be removed! It contains digons!!" << endl;
        return false;
    }

    polyMeshGenModifier meshModifier(mesh);

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const label nIntFaces = mesh.nInternalFaces();

    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    //- calculate point-faces
    VRWGraph pFaces;
    pFaces.reverseAddressing(mesh.points().size(), faces);

    labelLongList newPointLabel(mesh.points().size(), -1);

    //- find faces which shall be removed from the mesh
    //- as a consequence of sheet extraction
    boolList removeFace(faces.size(), false);
    boolList protectedFace(faces.size(), false);
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        if( sheetCells[owner[faceI]] && sheetCells[neighbour[faceI]] )
        {
            removeFace[faceI] = true;
        }
        else if
        (
            (sheetCells[owner[faceI]] && !sheetCells[neighbour[faceI]]) ||
            (!sheetCells[owner[faceI]] && sheetCells[neighbour[faceI]])
        )
        {
            protectedFace[faceI] = true;
        }
    }

    forAll(mesh.boundaries(), patchI)
    {
        const label start = mesh.boundaries()[patchI].patchStart();
        const label size = mesh.boundaries()[patchI].patchSize();

        for(label fI=0;fI<size;++fI)
        {
            const label faceI = start + fI;

            if( !sheetCells[owner[faceI]] )
                continue;

            const cell& c = cells[owner[faceI]];
            const label ofI = c.opposingFaceLabel(faceI, faces);

            if( removeFace[ofI] )
                removeFace[faceI] = true;
        }
    }

    //- make sure that the opposite faces to the protected faces
    //- are not selected for removal
    forAll(protectedFace, faceI)
    {
        if( protectedFace[faceI] )
        {
            if( sheetCells[owner[faceI]] )
            {
                const cell& c = cells[owner[faceI]];
                const label ofI = c.opposingFaceLabel(faceI, faces);

                removeFace[ofI] = false;
                removeFace[faceI] = false;
            }
            else if( sheetCells[neighbour[faceI]] )
            {
                const cell& c = cells[neighbour[faceI]];
                const label ofI = c.opposingFaceLabel(faceI, faces);

                removeFace[ofI] = false;
                removeFace[faceI] = false;
            }
        }
    }

    //- remove cells from the mesh and connect neighbouring sheets
    forAll(sheetCells, cellI)
    {
        if( !sheetCells[cellI] )
            continue;

        cell& c = cells[cellI];

        //- this cell will be collapsed
        //- find the first protected face
        //- merges it with the opposite face
        label posF(-1);

        forAll(c, fI)
        {
            if( protectedFace[c[fI]] )
            {
                posF = fI;
                break;
            }
        }

        if( posF == -1 )
        {
            //- this cell is removed without any re-connection
            //- of neighbouring cells
            forAll(c, fI)
                removeFace[c[fI]] = true;

            continue;
        }
        else
        {
            //- merge the protected face with the neighbouring face
            const label faceI = c[posF];
            const label ofI = c.opposingFaceLabel(faceI, faces);

            //- find neighbour cells over this face pair
            label neiF = owner[faceI];
            if( neiF == cellI )
                neiF = neighbour[faceI];

            label neiS = owner[ofI];
            if( neiS == cellI )
                neiS = neighbour[ofI];

            if( neiF == -1 )
            {
                label pos(-1);
                forAll(cells[neiS], fI)
                    if( ofI == cells[neiS][fI] )
                    {
                        pos = fI;
                        break;
                    }

                removeFace[ofI] = true;
                cells[neiS][pos] = faceI;
            }
            else if( neiS == -1 )
            {
                label pos(-1);
                forAll(cells[neiF], fI)
                    if( faceI == cells[neiF][fI] )
                    {
                        pos = fI;
                        break;
                    }

                removeFace[faceI] = true;
                cells[neiF][pos] = ofI;
            }
            else if( neiS > neiF )
            {
                label pos(-1);
                forAll(cells[neiS], fI)
                    if( cells[neiS][fI] == ofI )
                    {
                        pos = fI;
                        break;
                    }

                cells[neiS][pos] = faceI;
                removeFace[ofI] = true;

                if( faceI == neighbour[neiF] )
                    faces[faceI] = faces[faceI].reverseFace();
            }
            else if( neiF > neiS )
            {
                label pos(-1);
                forAll(cells[neiF], fI)
                    if( cells[neiF][fI] == faceI )
                    {
                        pos = fI;
                        break;
                    }

                cells[neiF][pos] = ofI;
                removeFace[faceI] = true;

                if( ofI == neighbour[neiS] )
                    faces[ofI] = faces[ofI].reverseFace();
            }
        }
    }

    //- re-connect faces near the removed sheet
    //- and move the points
    forAll(removeFace, faceI)
    {
        if( !removeFace[faceI] )
            continue;
        if( protectedFace[faceI] )
            continue;

        const face& f = faces[faceI];

        //- find the sheet cell comprising of this face
        //- and find the opposite face
        label cellI = owner[faceI];
        if( !sheetCells[cellI] )
            cellI = neighbour[faceI];

        const label ofI = cells[cellI].opposingFaceLabel(faceI, faces);

        if( protectedFace[ofI] )
            continue;

        //- check which edges are shared by protected faces
        FixedList<bool, 4> protectedEdge(false);

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            forAllRow(pFaces, e.start(), pfI)
            {
                const label fJ = pFaces(e.start(), pfI);

                if( !protectedFace[fJ] )
                    continue;

                const face& of = faces[fJ];

                if( of.which(e.end()) != -1 )
                {
                    protectedEdge[eI] = true;
                    protectedEdge[(eI+2)%4] = true;
                    break;
                }
            }
        }

        //- check if none of the edges are protected
        bool noneProtected(true);
        forAll(protectedEdge, i)
            if( protectedEdge[i] )
            {
                noneProtected = false;
                break;
            }

        if( !noneProtected )
        {
            forAll(f, eI)
            {
                if( protectedEdge[eI] )
                    continue;

                const edge e = f.faceEdge(eI);

                if(
                    (newPointLabel[e.start()] != -1) &&
                    (newPointLabel[e.end()] != -1)
                )
                    continue;

                const point newP = e.centre(points);

                newPointLabel[e.start()] = Foam::min(e.start(), e.end());
                points[e.start()] = newP;

                newPointLabel[e.end()] = newPointLabel[e.start()];
                points[e.end()] = newP;
            }
        }
        else
        {
            point newP(vector::zero);
            forAll(f, pI)
                newP += points[f[pI]];
            newP /= f.size();

            points[f[0]] = newP;
            forAll(f, pI)
            {
                const label pointI = f[pI];

                points[pointI] = newP;
                newPointLabel[pointI] = f[0];
            }
        }
    }

    //- re-connect faces
    forAll(newPointLabel, pointI)
    {
        if( newPointLabel[pointI] < 0 )
            continue;
        if( newPointLabel[pointI] == pointI )
            continue;

        forAllRow(pFaces, pointI, pfI)
        {
            face& pf = faces[pFaces(pointI, pfI)];
            const label pos = pf.which(pointI);
            pf[pos] = newPointLabel[pointI];
        }
    }

    //- remove cells from the mesh and the unused vertices
    meshModifier.removeFaces(removeFace);

    label nCells(0);
    forAll(sheetCells, cellI)
    {
        if( sheetCells[cellI] )
            continue;

        if( nCells < cellI )
            cells[nCells].transfer(cells[cellI]);

        ++nCells;
    }
    cells.setSize(nCells);

    meshModifier.removeUnusedVertices();

    return true;
}

bool extractSheet(polyMeshGen& mesh, const labelLongList& cellsInSheet)
{
    boolList sheetCells(mesh.cells().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for
    # endif
    forAll(cellsInSheet, i)
        sheetCells[cellsInSheet[i]] = true;

    return extractSheet(mesh, sheetCells);
}

bool extractSheet(polyMeshGen& mesh, const word& sheetCellSet)
{
    labelLongList sheetCells;
    const label cID = mesh.cellSubsetIndex(sheetCellSet);
    mesh.cellsInSubset(cID, sheetCells);

    return extractSheet(mesh, sheetCells);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hexHelpers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
