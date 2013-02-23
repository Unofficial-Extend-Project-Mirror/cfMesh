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
    
    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
		surfaceFile = ".."/surfaceFile;

    triSurf surf(registry.path()/surfaceFile);
	
	if( meshDict.found("subsetFileName") )
	{
		fileName subsetFileName = meshDict.lookup("subsetFileName");
        if( Pstream::parRun() )
			subsetFileName = ".."/subsetFileName;
		surf.readFaceSubsets(registry.path()/subsetFileName);
	}

	// construct the octree
    meshOctree moc(surf);
	meshOctreeCreator(moc, meshDict).createOctreeBoxes();

	polyMeshGen pmg(registry);
	pmg.read();
    
    // make sure that mesh contains all surface patches
    if(
        Pstream::parRun() &&
        (pmg.patchNames().size() != surf.patches().size())
    )
    {
        wordList patchNames(surf.patches().size());
        labelList patchStart(patchNames.size(), pmg.nInternalFaces());
        labelList nFacesInPatch(patchNames.size(), 0);
        
        forAll(patchNames, patchI)
            patchNames[patchI] = surf.patches()[patchI].name();
        
        forAll(pmg.patchNames(), patchI)
        {
            label pos(-1);
            forAll(patchNames, nameI)
                if( patchNames[nameI] == pmg.patchNames()[patchI] )
                {
                    pos = nameI;
                    break;
                }
                
            nFacesInPatch[pos] = pmg.numFacesInPatch()[patchI];
            patchStart[pos] = pmg.patchStart()[patchI];
                
            for(label i=pos+1;i<patchNames.size();++i)
                patchStart[i] += nFacesInPatch[pos];
        }
        
        polyMeshGenModifier meshModifier(pmg);
        meshModifier.patchNamesAccess() = patchNames;
        meshModifier.patchStartAccess() = patchStart;
        meshModifier.nFacesInPatchAccess() = nFacesInPatch;
    }
	
	meshSurfaceEngine ms(pmg);
	meshSurfaceOptimizer mo(ms, moc);

    /*
	labelList pointRegion(pmg.points().size(), -1);
	for(label faceI=mesh.nInternalFaces();faceI<pmg.faces().size();faceI++)
	{
		const face& f = pmg.faces()[faceI];
		
		forAll(f, pI)
			pointRegion[f[pI]] = 1;
    }
    
	mo.preOptimizeSurface(pointRegion);
	*/
    mo.optimizeSurface();
    
    pmg.write();
	writeMeshEnsight(pmg, "smoothMesh");
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
