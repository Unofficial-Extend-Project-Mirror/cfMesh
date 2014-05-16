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

    if( outFileName == inFileName )
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite input file."
            << exit(FatalError);
    }

    IFstream surfFile(inFileName);
    
    //- parse the number of vertices
    token t;
    do
    {
        surfFile >> t;
    } while( !t.isLabel() );
    pointField points(t.labelToken());
    Info << "Surface has " << points.size() << " points" << endl;
    
    //- find the number of elements
    do
    {
        surfFile >> t;
    } while( !t.isLabel() );
    List<labelledTri> triFaces(t.labelToken());
    Info << "Surface has " << triFaces.size() << " facets" << endl;
    
    //- read vertices
    do
    {
        surfFile >> t;
    } while( t.isWord() && (t.wordToken() != "end_header") );
    
    forAll(points, pointI)
    {
        point& p = points[pointI];
        
        surfFile >> p.x();
        surfFile >> p.y();
        surfFile >> p.z();
    }
    
    //- read triangles
    forAll(triFaces, triI)
    {
        label nPts;
        surfFile >> nPts;
        if( nPts != 3 )
        {
            Info << "Facet " << triI << " is not a triangle!!" << endl;
            Warning << "Cannot convert this surface!" << endl;
            return 0;
        }
        
        for(label i=0;i<nPts;++i)
            surfFile >> triFaces[triI][i];
        
        triFaces[triI].region() = 0;
    }

    triSurface surf(triFaces, points);

    Info << "Writing : " << outFileName << endl;
    surf.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
