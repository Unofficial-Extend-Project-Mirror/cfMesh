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

Application
    Test for boundary layers

Description
    -

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMeshGen.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "polyMeshGenChecks.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //boundaryLayers bndLayers(pmg);
    //bndLayers.addLayerForPatch("inlet");
    //bndLayers.addLayerForPatch("symmetryplane");
    //bndLayers.createOTopologyLayers();
    //bndLayers.addLayerForAllPatches();

    Info << "Starting bnd layer refinement "
         << runTime.elapsedClockTime() << endl;

    //polyMeshGenChecks::checkMesh(pmg, true);

    refineBoundaryLayers refLayers(pmg);

    refineBoundaryLayers::readSettings(meshDict, refLayers);

    refLayers.refineLayers();

    Info << "Finished with bnd layer refinement "
         << runTime.elapsedClockTime() << endl;

    polyMeshGenChecks::checkMesh(pmg, true);
    return 0;
    pmg.write();
    //meshOctree* octreePtr = NULL;
    //meshOptimizer(*octreePtr, pmg).preOptimize();

    //pmg.addressingData().checkMesh(true);

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
