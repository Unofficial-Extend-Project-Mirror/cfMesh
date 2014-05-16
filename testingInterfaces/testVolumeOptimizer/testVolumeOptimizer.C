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
#include "volumeOptimizer.H"
#include "Time.H"
#include "polyMeshGen.H"
#include "partTetMeshSimplex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    DynList<labelledPoint> points;

    points.append(labelledPoint(0, point(-0.5, -0.5, -1))); // point which gets moved

    points.append(labelledPoint(1, point(0, 0, 0)));
    points.append(labelledPoint(2, point(1, 0, 0)));
    points.append(labelledPoint(3, point(1, 1, 0)));
    points.append(labelledPoint(4, point(0, 1, 0)));

    points.append(labelledPoint(5, point(0, 0, 1)));
    points.append(labelledPoint(6, point(1, 0, 1)));
    points.append(labelledPoint(7, point(1, 1, 1)));
    points.append(labelledPoint(8, point(0, 1, 1)));

    points.append(labelledPoint(9, point(0.5, 0.25, 0.999)));

    Info << "Points before movement " << points << endl;

    DynList<parPartTet> tets;
    //- tets at z = 0
    tets.append(parPartTet(points[1], points[2], points[9], points[0]));
    tets.append(parPartTet(points[2], points[3], points[9], points[0]));
    tets.append(parPartTet(points[3], points[4], points[9], points[0]));
    tets.append(parPartTet(points[4], points[1], points[9], points[0]));

    //- tets at z = 1
    tets.append(parPartTet(points[6], points[5], points[8], points[0]));
    tets.append(parPartTet(points[6], points[8], points[7], points[0]));

    //- tets at y = 0
    tets.append(parPartTet(points[2], points[1], points[5], points[0]));
    tets.append(parPartTet(points[5], points[6], points[2], points[0]));

    //- tets at y = 1
    tets.append(parPartTet(points[4], points[3], points[7], points[0]));
    tets.append(parPartTet(points[7], points[8], points[4], points[0]));

    //- tets at x = 0
    tets.append(parPartTet(points[5], points[1], points[4], points[0]));
    tets.append(parPartTet(points[5], points[4], points[8], points[0]));

    //- tets at x = 1
    tets.append(parPartTet(points[2], points[6], points[7], points[0]));
    tets.append(parPartTet(points[2], points[7], points[3], points[0]));

    partTetMeshSimplex simplex(tets, 0);

    volumeOptimizer volOpt(simplex);

    volOpt.optimizeNodePosition(1e-4);

    Info << "Points after movement " << simplex.pts() << endl;

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
