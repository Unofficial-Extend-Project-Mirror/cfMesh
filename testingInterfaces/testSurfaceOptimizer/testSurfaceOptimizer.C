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
    Test for smoothers

Description
    - reads the mesh and tries to untangle negative volume cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "surfaceOptimizer.H"
#include "Time.H"
#include "polyMeshGen.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    DynList<point> points;
    DynList<triFace> trias;

    points.append(point(-0.5, 1.2, 0)); // point which gets moved

    points.append(point(0, 0, 0));
    points.append(point(1, 0, 0));
    points.append(point(1, 1, 0));
    points.append(point(0.5, -0.01, 0));
    points.append(point(0, 1, 0));

    Info << "Points before movement " << points << endl;

    trias.append(triFace(0, 1, 2));
    trias.append(triFace(0, 2, 3));
    trias.append(triFace(0, 3, 4));
    trias.append(triFace(0, 4, 5));
    trias.append(triFace(0, 5, 1));

    surfaceOptimizer sOpt(points, trias);
    sOpt.optimizePoint(1e-3);

    Info << "Points after movement " << points << endl;

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
