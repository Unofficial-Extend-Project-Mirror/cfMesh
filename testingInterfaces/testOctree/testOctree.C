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
    Creates the octree

Description
    - creates an octree based on the information specified in a dictionary

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeAutomaticRefinement.H"
#include "meshOctreeCreator.H"
#include "Time.H"
#include "objectRegistry.H"
#include "triSurf.H"
#include "writeOctreeEnsight.H"

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

    Info << "Path " << runTime.path() << endl;
    Info << "Root path " << runTime.rootPath() << endl;
    Info << "Case name " << runTime.caseName() << endl;

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surf(runTime.path()/surfaceFile);

    if( meshDict.found("subsetFileName") )
    {
        fileName subsetFileName = meshDict.lookup("subsetFileName");
        if( Pstream::parRun() )
            subsetFileName = ".."/subsetFileName;
        surf.readFaceSubsets(runTime.path()/subsetFileName);
    }

    surf.edgeFaces();
    surf.faceEdges();

    const scalar startTime = runTime.elapsedClockTime();

    meshOctree* moPtr = new meshOctree(surf);
    meshOctreeCreator* mcPtr = new meshOctreeCreator(*moPtr, meshDict);
    mcPtr->createOctreeBoxes();
    deleteDemandDrivenData(mcPtr);

    for(label i=1;i<moPtr->numberOfLeaves()-1;++i)
    {
        if( (moPtr->returnLeaf(i+1).coordinates() <= moPtr->returnLeaf(i).coordinates()) )
            FatalError << "1. Wrong ordering "
            << moPtr->returnLeaf(i+1).coordinates()
            << " and " << moPtr->returnLeaf(i).coordinates()
            << abort(FatalError);
        if( !(moPtr->returnLeaf(i).coordinates() >= moPtr->returnLeaf(i-1).coordinates()) )
            FatalError << "2. Wrong ordering"
            << moPtr->returnLeaf(i).coordinates()
            << " and " << moPtr->returnLeaf(i-1).coordinates()
            << abort(FatalError);
    }
    /*  if( Pstream::parRun() )
        {
            std::ostringstream ss;

            ss << Pstream::myProcNo();

            //writeOctreeEnsight(*moPtr, "insideBoxes"+ss.str(), meshOctreeCube::INSIDE);
            //writeOctreeEnsight(*moPtr, "unknownBoxes"+ss.str(), meshOctreeCube::UNKNOWN);
            //writeOctreeEnsight(*moPtr, "outsideBoxes"+ss.str(), meshOctreeCube::OUTSIDE);
            writeOctreeEnsight(*moPtr, "dataBoxes"+ss.str(), meshOctreeCube::DATA);
        }
        else
        {
            writeOctreeEnsight(*moPtr, "insideBoxes", meshOctreeCube::INSIDE);
            //writeOctreeEnsight(*moPtr, "unknownBoxes", meshOctreeCube::UNKNOWN);
            writeOctreeEnsight(*moPtr, "outsideBoxes", meshOctreeCube::OUTSIDE);
            writeOctreeEnsight(*moPtr, "dataBoxes", meshOctreeCube::DATA);
        }
    */
    deleteDemandDrivenData(moPtr);

    Info<< "Execution time for octree creation = "
        << (runTime.elapsedClockTime() - startTime)
        << " s\n" << endl << endl;

    label sizeOf = sizeof(meshOctreeCube);
    Info << "Size of single cube is " << sizeOf << endl;
    sizeOf = sizeof(direction);
    Info << "Size of direction is " << sizeOf << endl;
    sizeOf = sizeof(labelList);
    Info << "Size of labelList is " << sizeOf << endl;
    sizeOf = sizeof(UList<direction>);
    Info << "Size of UList<direction> is " << sizeOf << endl;
    sizeOf = sizeof(FixedList<direction, 4>);
    Info << "Size of FixedList<direction, 4> is " << sizeOf << endl;
    sizeOf = sizeof(meshOctreeCube*);
    Info << "Size of meshOctreeCube* is " << sizeOf << endl;
    sizeOf = sizeof(meshOctreeCubeCoordinates);
    Info << "Size of meshOctreeCubeCoordinates " << sizeOf << endl;
    sizeOf = sizeof(meshOctreeCubeBasic);
    Info << "Size of meshOctreeCubeBasic " << sizeOf << endl;

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
