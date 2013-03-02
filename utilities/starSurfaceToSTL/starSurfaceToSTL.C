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

Description
    Converts STAR CD surfaces into stl or other formats

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

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
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if (outFileName == inFileName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite input file."
            << exit(FatalError);
    }

    fileName vFile(inFileName+".vrt");

    IFstream vertexFile(vFile.c_str());

    pointField points(100);
    label index(0);

    bool finished(false);

    do
    {
        token t(vertexFile);
        if( t.isLabel() )
        {
            point p;
            p.x() = readScalar(vertexFile);
            p.y() = readScalar(vertexFile);
            p.z() = readScalar(vertexFile);

            points.newElmt(index++) = p;
        }
        else
        {
            finished = true;
        }
    } while( !finished );
   
    points.setSize(index);

    Info << "points " << points << endl;

    index = 0;
    SLList<labelledTri> triangles;

    fileName cFile(inFileName+".cel");
    IFstream cellFile(cFile);

    finished = false;
    do
    {
        token t(cellFile);
        if( t.isLabel() )
        {
            label v0 = readLabel(cellFile) - 1;
            label v1 = readLabel(cellFile) - 1;
            label v2 = readLabel(cellFile) - 1;
            label v3 = readLabel(cellFile) - 1;
            readLabel(cellFile);
            readLabel(cellFile);
            readLabel(cellFile);
            readLabel(cellFile);

            if( v2 > 0 && v3 > 0 )
            {
                if( v2 == v3 )
                {
                    triangles.append(labelledTri(v0,v1,v2,0));
                }
                else
                {
                    triangles.append(labelledTri(v0,v1,v2,0));
                    triangles.append(labelledTri(v0,v2,v3,0));
                }
            }

            readLabel(cellFile);
            readLabel(cellFile);
        }
        else
        {
            finished = true;
        }
    } while( !finished );

    List<labelledTri> triFaces(triangles);

    triSurface surf(triFaces, points);

    Info << "Writing : " << outFileName << endl;
    surf.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
