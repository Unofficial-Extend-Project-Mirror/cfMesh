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

Description
    Reads the specified surface and writes it in the fms format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "coordinateModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

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

    triSurf surface(runTime.path()/surfaceFile);

//    argList::noParallel();
//    argList::validArgs.clear();
//    argList::validArgs.append("output surface file");

//    argList args(argc, argv);

    //const fileName outFileName(args.args()[1]);

    triSurfModifier sMod(surface);
    pointField& pts = sMod.pointsAccess();

    coordinateModifier cMod(meshDict.subDict("geometryModification"));

    //- transform points
    forAll(pts, i)
    {
        pts[i] = cMod.modifiedPoint(pts[i]);
    }

    Info << "Writting transformed surface" << endl;
    surface.writeSurface("transformedSurf.fms");

    //- apply backward transformation
    forAll(pts, i)
        pts[i] = cMod.backwardModifiedPoint(pts[i]);

    Info << "Writting backward transformed surface" << endl;
    surface.writeSurface("backwardTransformedPoints.fms");

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
