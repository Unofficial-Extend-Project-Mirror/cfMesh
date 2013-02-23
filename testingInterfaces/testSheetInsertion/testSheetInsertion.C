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

Application
    Test for smoothers

Description
    - reads the mesh and tries to untangle negative volume cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "hexHelpers.H"
#include "Time.H"
#include "polyMeshGen.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    Info << "Finished reading mesh" << endl;

    hexHelpers::findColumnCells(pmg, 0, "column_0");

    Info << "Finding sheet 1" << endl;

    edge e = pmg.faces()[0].faceEdge(0);
    hexHelpers::findSheetCells(pmg, e, "sheet_0");

    Info << "Finding sheet 2" << endl;

    e = pmg.faces()[0].faceEdge(1);
    hexHelpers::findSheetCells(pmg, e, "sheet_1");

    Info << "Collapsing column" << endl;
    forAll(pmg.cells(), cellI)
    {
        const cell& c = pmg.cells()[cellI];

        bool allInternal(true);
        forAll(c, fI)
            if( c[fI] >= pmg.nInternalFaces() )
            {
                allInternal = false;
                break;
            }

        if( allInternal )
        {
            hexHelpers::collapseColumn(pmg, c[1], 1);
            //break;
        }
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
