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
    Test for boundary layers

Description
    - 

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "boundaryLayers.H"
#include "writeMeshEnsight.H"
#include "meshOptimizer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"
	
	objectRegistry registry(runTime);
	
/*	labelList patchStarts(mesh.boundaryMesh().size());
	labelList patchSizes(mesh.boundaryMesh().size());
	
	forAll(mesh.boundaryMesh(), patchI)
	{
		patchStarts[patchI] = mesh.boundaryMesh()[patchI].start();
		patchSizes[patchI] = mesh.boundaryMesh()[patchI].size();
	}
	
	polyMeshGen pmg
	(
		registry,
		mesh.points(),
		mesh.faces(),
		mesh.cells(),
		mesh.boundaryMesh().names(),
		patchStarts,
		patchSizes
	);
*/
	polyMeshGen pmg(registry);
	pmg.read();
	//writeMeshEnsight(pmg, "meshWithoutBndLayers");
	
	boundaryLayers bndLayers(pmg);
	//bndLayers.addLayerForPatch("inlet");
	//bndLayers.addLayerForPatch("symmetryplane");
	//bndLayers.createOTopologyLayers();
	bndLayers.addLayerForAllPatches();
	
	//pmg.write();
	//meshOctree* octreePtr = NULL;
	//meshOptimizer(*octreePtr, pmg).preOptimize();
	
	writeMeshEnsight(pmg, "meshWithBndLayers");
	//pmg.addressingData().checkMesh(true);
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
