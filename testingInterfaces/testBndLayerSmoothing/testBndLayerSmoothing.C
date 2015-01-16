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
#include "meshOptimizer.H"
#include "boundaryLayers.H"
#include "Time.H"
#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "meshOptimizer.H"
#include "boolList.H"

#include "boundaryLayerOptimisation.H"

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

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( help::isnan(p) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    //boundaryLayers(pmg).addLayerForAllPatches();
    //pmg.clearAddressingData();

    meshSurfaceEngine mse(pmg);
    boundaryLayerOptimisation blOpt(pmg, mse);

    const labelList& faceCell = mse.faceOwners();
    const boolList& isBaseFace = blOpt.isBaseFace();

    labelLongList layerCells;
    boolList layerCell(pmg.cells().size(), false);
    const label blCellsId = pmg.addCellSubset("boundaryLayerCells");
    forAll(isBaseFace, bfI)
    {
        pmg.addCellToSubset(blCellsId, faceCell[bfI]);
        layerCell[faceCell[bfI]] = true;
        layerCells.append(faceCell[bfI]);
    }

    //- find points in boundary layer
    boolList pointInBoundaryLayer(pmg.points().size(), false);
    const cellListPMG& cells = pmg.cells();
    const faceListPMG& faces = pmg.faces();
    forAll(cells, cellI)
    {
        if( layerCell[cellI] )
        {
            const cell& c = cells[cellI];
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                    pointInBoundaryLayer[f[pI]] = true;
            }
        }
    }

    //- check if there exist faces with all vertices in the boundary layer
    const label allPointsInLayerId = pmg.addCellSubset("allPointsInLayer");
    forAll(cells, cellI)
    {
        if( !layerCell[cellI] )
        {
            bool allInBndLayer(true);

            const cell& c = cells[cellI];
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                    allInBndLayer = false;
            }

            if( allInBndLayer )
            {
                Info << "Cell " << cellI
                     << " has all points in the layer" << endl;
                pmg.addCellToSubset(allPointsInLayerId, cellI);
            }
        }
    }

    blOpt.optimiseLayer(5, 0.15, 0.4);

    pmg.clearAddressingData();

    meshOptimizer mOpt(pmg);
    mOpt.lockCellsInSubset("boundaryLayerCells");
    mOpt.untangleMeshFV(5, 50, 0);
    //meshOptimizer mOpt(pmg);
    //mOpt.optimizeBoundaryLayer();

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( help::isnan(p) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
