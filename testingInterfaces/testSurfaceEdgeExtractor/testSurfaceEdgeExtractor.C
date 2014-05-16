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
    Test for edge extraction

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEdgeExtractor.H"
#include "triSurf.H"

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
            registry.time().constant(),
            registry,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const fileName surfaceFile = meshDict.lookup("surfaceFile");

    triSurf surf(registry.path()/surfaceFile);

    // construct the octree
    meshOctree mo(surf);
    meshOctreeCreator(mo, meshDict).createOctreeWithRefinedBoundary(8);

    // construct the polyMeshGen
    polyMeshGen pmg(registry);
    pmg.read();

    // find regions for boundary vertices
    meshSurfaceEngine mse(pmg);

    labelList pointRegion(pmg.points().size(), -1);

    Info << "Finding boundary regions for boundary vertices" << endl;
    const labelList& bPoints = mse.boundaryPoints();
    forAll(bPoints, bpI)
    {
        point& p = pmg.points()[bPoints[bpI]];
        point np;
        mo.findNearestSurfacePoint(np, pointRegion[bPoints[bpI]], p);
        p = np;
    }

    Info << "Extracting edges" << endl;
    meshSurfaceEdgeExtractor(pmg, mo, pointRegion);

    pmg.addressingData().checkMesh(true);

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
