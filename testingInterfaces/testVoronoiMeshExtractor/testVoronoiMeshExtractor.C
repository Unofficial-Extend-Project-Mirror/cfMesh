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
    Test of the voronoi mesh

Description
    - creates an octree and creates voronoi mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "voronoiMeshExtractor.H"
#include "triSurf.H"
#include "writeMeshEnsight.H"

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
	
	const fileName surfaceFile = meshDict.lookup("surfaceFile");

    triSurf surf(registry.path()/surfaceFile);

	// construct the octree
    meshOctree mo(surf);
	meshOctreeCreator(mo, meshDict).createOctreeBoxes();
	
	Info<< "Execution time for octree creation = "
        << runTime.elapsedCpuTime()
        << " s\n" << endl << endl;
	
	polyMeshGen pmg(registry);
	voronoiMeshExtractor vmg(mo, meshDict, pmg);
	
	vmg.createMesh();
	
	pmg.write();
	writeMeshEnsight(pmg, "voronoiMesh");
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
