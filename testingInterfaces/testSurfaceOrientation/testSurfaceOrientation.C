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
    Test the triSurfaceClassifyEdges

Description
    - creates an octree and checks which edges in the surface mesh are convex
     or concave

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "meshOctreeAutomaticRefinement.H"
#include "triSurfaceClassifyEdges.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "cartesianMeshExtractor.H"
#include "triSurf.H"

#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

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

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surf(runTime.path()/surfaceFile);

    // construct the octree
    meshOctree mo(surf);
    meshOctreeCreator(mo, meshDict).createOctreeBoxes();

    meshOctreeAutomaticRefinement(mo, meshDict, false).automaticRefinement();

    Info<< "Execution time for octree creation = "
        << runTime.elapsedCpuTime()
        << " s\n" << endl << endl;

    triSurfaceClassifyEdges tce(mo);
    const List<direction>& edgeTypes = tce.edgeTypes();
    const edgeLongList& edges = surf.edges();

    const label featureNodes = surf.addPointSubset("featureNodes");
    const label convexId = surf.addPointSubset("convexEdges");
    const label concaveId = surf.addPointSubset("concaveEdges");

    forAll(edgeTypes, edgeI)
    {
        if( !(edgeTypes[edgeI] & triSurfaceClassifyEdges::FEATUREEDGE) )
            continue;

        surf.addPointToSubset(featureNodes, edges[edgeI].start());
        surf.addPointToSubset(featureNodes, edges[edgeI].end());

        if( edgeTypes[edgeI] & triSurfaceClassifyEdges::CONVEXEDGE )
        {
            surf.addPointToSubset(convexId, edges[edgeI].start());
            surf.addPointToSubset(convexId, edges[edgeI].end());
        }
        else if( edgeTypes[edgeI] & triSurfaceClassifyEdges::CONCAVEEDGE )
        {
            surf.addPointToSubset(concaveId, edges[edgeI].start());
            surf.addPointToSubset(concaveId, edges[edgeI].end());
        }
    }

    Info << "Writting surface to a file" << endl;
    surf.writeSurface("surfWithEdgeTypes.fms");

/*    if( Pstream::parRun() )
    {
        std::ostringstream ss;
        ss << Pstream::myProcNo();
        writeMeshEnsight(pmg, "cartesianMesh"+ss.str());
    }
    else
    {
        writeMeshEnsight(pmg, "cartesianMesh");
    }
*/
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
