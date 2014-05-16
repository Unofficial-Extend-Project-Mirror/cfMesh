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
    Test for smoothers

Description
    - reads the mesh and tries to untangle negative volume cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOptimizer.H"
#include "Time.H"
#include "polyMeshGen.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    Info << "Finished reading mesh" << endl;

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( isnan(p.x()) || isnan(p.y()) || isnan(p.z()) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    meshOptimizer mo(pmg);

    mo.optimizeMeshFV();
    //mo.optimizeLowQualityFaces();

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( isnan(p.x()) || isnan(p.y()) || isnan(p.z()) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
