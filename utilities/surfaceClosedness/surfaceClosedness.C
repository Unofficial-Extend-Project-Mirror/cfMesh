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
    Checks if the surface is geometrically closed

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

#define DEBUGStitch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{

    argList::noParallel();
    argList::validArgs.clear();
    argList::validOptions.insert("noCleanup", "");
    argList::validOptions.insert("group", "");
    argList::validArgs.append("input surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);

    triSurface ts(inFileName);

    const List<labelledTri>& triangles = ts.localFaces();
    const pointField& points = ts.localPoints();
    const labelListList& edgeFaces = ts.edgeFaces();

    forAll(edgeFaces, eI)
        if( edgeFaces[eI].size() != 2 )
            Info << "Edge " << eI << " is an open edge " << endl;


    vector sum(vector::zero);

    forAll(triangles, tI)
    {
        const vector n = triangles[tI].normal(points);

        sum += n;
    }

    if( mag(sum) > 1e-12 )
    {
        Info << "Surface is not closed " << sum << endl;
    }
    else
    {
        Info << "Surface Ok" << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
