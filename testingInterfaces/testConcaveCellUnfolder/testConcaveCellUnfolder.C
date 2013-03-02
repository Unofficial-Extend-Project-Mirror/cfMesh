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
    Test of the octree dual

Description
    - Testing interface for correcting concave cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "dualUnfoldConcaveCells.H"
#include "meshSurfaceOptimizer.H"
#include "triSurf.H"
#include "writeMeshEnsight.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"
    
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
    meshOctreeCreator(mo, meshDict).createOctreeWithRefinedBoundary(8);
    
    // construct the polyMeshGen
    labelList starts(mesh.boundaryMesh().names().size());
    labelList nFacesInPatch(starts.size());
    
    forAll(mesh.boundaryMesh(), patchI)
    {
        starts[patchI] = mesh.boundaryMesh()[patchI].start();
        nFacesInPatch[patchI] = mesh.boundaryMesh()[patchI].size();
    }
    
    polyMeshGen pmg
    (
        registry,
        mesh.points(),
        mesh.faces(),
        mesh.cells(),
        mesh.boundaryMesh().names(),
        starts,
        nFacesInPatch
    );
    
    //- optimize surface to get rid of nearby vertices
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(pmg);
    meshSurfaceOptimizer(*msePtr, mo).optimizeSurface();
    deleteDemandDrivenData(msePtr);

    
    dualUnfoldConcaveCells(pmg, mo).unfoldInvalidCells();
    
    pmg.write();
    writeMeshEnsight(pmg, "correctedEdges");
    pmg.addressingData().checkMesh(true);
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
