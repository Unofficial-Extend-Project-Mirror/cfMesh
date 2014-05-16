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
    Test of surface smoother

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngine.H"
#include "Time.H"
#include "objectRegistry.H"
#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "polyMeshGen.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    objectRegistry registry(runTime);

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            registry.time().system(),
            registry,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surf(registry.path()/surfaceFile);

    // construct the octree
    meshOctree moc(surf);
    meshOctreeCreator(moc, meshDict).createOctreeBoxes();

    polyMeshGen pmg(registry);
    pmg.read();


    mo.optimizeSurface();
    mo.untangleSurface();

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
