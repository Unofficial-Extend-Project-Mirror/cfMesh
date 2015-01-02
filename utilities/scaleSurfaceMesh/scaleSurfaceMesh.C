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
#include "checkMeshDict.H"
#include "surfaceMeshGeometryModification.H"

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

    checkMeshDict cmd(meshDict);

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surface(runTime.path()/surfaceFile);

    surfaceMeshGeometryModification gMod(surface, meshDict);

    //- modify points
    const triSurf* modSurfPtr = gMod.modifyGeometry();

    Info << "Writting modified surface" << endl;
    modSurfPtr->writeSurface("modifiedSurf.stl");

    # ifdef DEBUGScaling
    //- apply backward modification
    Info << "Here" << endl;
    coordinateModifier cMod(meshDict.subDict("anisotropicSources"));
    Info << "Starting modifications" << endl;
    forAll(surface.points(), i)
    {
        Info << "\nOrig point " << i << " coordinates " << surface.points()[i]
             << " modified point " << modSurfPtr->points()[i] << endl;
        const point p = cMod.backwardModifiedPoint(modSurfPtr->points()[i]);

        if( mag(p - surface.points()[i]) > 1e-14 )
        {
            Warning << "Point " << i << " is different "
                    << p
                    << " from original " << surface.points()[i]
                    << " modified point "
                    << cMod.modifiedPoint(surface.points()[i]) << endl;
            ::exit(0);
        }
    }

    Info << "Backscaling Ok" << endl;
    ::exit(0);
    # endif

    surfaceMeshGeometryModification bgMod(*modSurfPtr, meshDict);
    const triSurf* backModSurfPtr = bgMod.revertGeometryModification();

    Info << "Writting backward transformed surface" << endl;
    backModSurfPtr->writeSurface("backwardModifiedSurf.stl");

    # ifdef DEBUGScaling
    forAll(backModSurfPtr->points(), pI)
        if( mag(backModSurfPtr->points()[pI] - surface.points()[pI]) > 1e-14 )
            Warning << "Point " << pI << " is different "
                    << backModSurfPtr->points()[pI]
                    << " from original " << surface.points()[pI] << endl;
    # endif

    deleteDemandDrivenData(modSurfPtr);
    deleteDemandDrivenData(backModSurfPtr);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
