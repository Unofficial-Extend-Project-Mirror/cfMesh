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
    Test of the cartesian mesh

Description
    - creates an octree and creates cartesian mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "meshOctreeAutomaticRefinement.H"
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
    meshOctreeCreator moc(mo, meshDict);
    moc.createOctreeBoxes();

    //meshOctreeAutomaticRefinement(mo, meshDict, false).automaticRefinement();

    Info<< "Execution time for octree creation = "
        << runTime.elapsedCpuTime()
        << " s\n" << endl << endl;

    polyMeshGen pmg(runTime);
    cartesianMeshExtractor cmg(mo, meshDict, pmg);

    //cmg.decomposeSplitHexes();
    cmg.createMesh();

    const pointFieldPMG& points = pmg.points();
/*    forAll(points, pointI)
    {
        const point p = points[pointI];
        for(label pointJ=pointI+1;pointJ<points.size();++pointJ)
            if( mag(p - points[pointJ]) < 1e-10 )
                Info << "Points " << pointI << " and " << pointJ
                     << " are duplicates " << endl;
    }
*/
    boolList usedPoint(points.size());
    const faceListPMG& faces = pmg.faces();
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        forAll(f, pI)
            usedPoint[f[pI]] = true;
    }

    forAll(usedPoint, pointI)
        if( !usedPoint[pointI] )
            Info << "Point " << pointI << " is not used !!" << endl;

    pmg.write();

/*    if( Pstream::parRun() )
    {
        std::ostringstream ss;
        ss << Pstream::myProcNo();
    }
*/
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
