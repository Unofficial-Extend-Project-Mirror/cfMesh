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

Class


Description
    Functions for detection sheets and columns in the hexhedral mesh,
    and the utilities whih can modify the topology of the hexahedral mesh

SourceFiles
    hexHelpers.C

\*---------------------------------------------------------------------------*/

#ifndef hexHelpers_H
#define hexHelpers_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace hexHelpers Declaration
\*---------------------------------------------------------------------------*/

namespace hexHelpers
{
    //- check if the supplied mesh contains hexahedra, only
    bool allHexMesh(const polyMeshGen& mesh);

    //- find all cells in a column containing a given face
    label findColumnCells
    (
        const polyMeshGen& mesh,
        const label faceI,
        labelLongList& columnCells
    );
    label findColumnCells
    (
        polyMeshGen& mesh,
        const label faceI,
        const word& columnCellSet
    );

    //- check if the column self-intersects itself
    //- returns true in case of existing self-intersections
    bool selfIntersectingColumn
    (
        const polyMeshGen& mesh,
        const boolList& columnCells
    );
    bool selfIntersectingColumn
    (
        const polyMeshGen& mesh,
        const labelLongList& columnCells
    );
    bool selfIntersectingColumn
    (
        const polyMeshGen& mesh,
        const word& columnCellSet
    );

    //- find all columns in hex mesh
    label findAllColumns(polyMeshGen& mesh, const word& cellSetPrefix);

    //- find all cells in a sheet containing a given face
    label findSheetCells
    (
        const polyMeshGen& mesh,
        const edge& sheetEdge,
        labelLongList& cellsInSheet
    );
    label findSheetCells
    (
        polyMeshGen& mesh,
        const edge& sheetEdge,
        const word& sheetCellSet
    );

    //- check if the sheet has any digons
    bool hasSheetDigons(const polyMeshGen& mesh, const boolList& sheetCells);
    bool hasSheetDigons
    (
        const polyMeshGen& mesh,
        const labelLongList& sheetCells
    );
    bool hasSheetDigons(const polyMeshGen& mesh, const word& sheetCellSet);

    //- find all sheets in hex mesh
    void findAllSheets(polyMeshGen& mesh, const word& cellSetPrefix);

    // Modification tools
    //- find all cells creating a single column and colapse the column
    //- the column is collapsed by merging the point at the given position
    //- in the selected face with the opposite point
    bool collapseColumn
    (
        polyMeshGen& mesh,
        const label faceI,
        const label positionInFace
    );

    //- collapse the column by merging the vertices in the layer of edges
    //- containing the given point with the opposite layer of edges
    bool collapseColumn
    (
        polyMeshGen& mesh,
        const labelLongList& columnCells,
        const label pointI
    );
    bool collapseColumn
    (
        polyMeshGen& mesh,
        const word& columnCellSet,
        const label pointI
    );

    //- find all cells creating a sheet containing a given face and extract
    //- the sheet into a set of faces
    bool extractSheet(polyMeshGen& mesh, const edge& sheetEdge);

    //- extract the sheet of cells
    bool extractSheet(polyMeshGen& mesh, const boolList& sheetCells);
    bool extractSheet(polyMeshGen& mesh, const labelLongList& cellsInSheet);
    bool extractSheet(polyMeshGen& mesh, const word& sheetCellSet);

} // End namespace hexHelpers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
