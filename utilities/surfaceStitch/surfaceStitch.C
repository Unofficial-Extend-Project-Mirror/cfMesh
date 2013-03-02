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
    Closing holes in the triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

#include "DynList.H"
#include "LongList.H"

//#define DEBUGStitch

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

    triSurface* ts = new triSurface(inFileName);

    label counter(0);

    LongList<labelledTri> createdTriangles;
    do
    {
        Info << "Creating new triangle" << endl;
        createdTriangles.clear();

        triSurface& surf = *ts;

        const pointField& points = surf.localPoints();
        const labelListList& eFaces = surf.edgeFaces();
        const edgeList& edges = surf.edges();
        const labelListList& pEdges = surf.pointEdges();
        const labelListList& pFaces = surf.pointFaces();

        //- mark the open edges
        boolList problemEdges(eFaces.size(), false);

        forAll(eFaces, eI)
        {
            if( eFaces[eI].size() < 2 )
                problemEdges[eI] = true;
        }
        
        //- find the vertices which have exactly two open edges attached to them
        //- create triangles which fit best into their surrounding
        forAll(pEdges, pointI)
        {
            const labelList& pe = pEdges[pointI];
            
            DynList<label> openEdges;
            forAll(pe, peI)
            {
                if( problemEdges[pe[peI]] )
                    openEdges.append(pe[peI]);
            }
            
            //- skip vertices attached to more than two open edges
            if( openEdges.size() != 2 )
                continue;
            
            //- skip vertices attached to only one triangle
            if( pFaces[pointI].size() == 1 )
                continue;
            
            const edge& e0 = edges[openEdges[0]];
            const edge& e1 = edges[openEdges[1]];
            
            const labelledTri& t0 = surf[eFaces[openEdges[0]][0]];
            const labelledTri& t1 = surf[eFaces[openEdges[1]][0]];
            
            labelledTri t;
            t[0] = pointI;
            if( e0.start() == pointI )
            {
                t[1] = e1.otherVertex(pointI);
                t[2] = e0.otherVertex(pointI);
                t.region() = t0.region();
            }
            else if( e1.start() == pointI )
            {
                t[1] = e0.otherVertex(pointI);
                t[2] = e1.otherVertex(pointI);
                t.region() = t1.region();
            }
            else
            {
                continue;
                FatalError << "Strange" << exit(FatalError);
            }
            
            Info << "Creating triangle for point " << pointI << endl;
            Info << "triangle is " << t << endl;
            
            //- check dot product between the normals
            vector n0 = t0.normal(points);
            n0 /= mag(n0) + VSMALL;
            
            vector n1 = t1.normal(points);
            n1 /= mag(n1) + VSMALL;
            
            vector n = t.normal(points);
            n /= mag(n) + VSMALL;
            
            //- find the maximum dot product between the normals of
            //- of existing triangles and the newly generated one
            scalar q = (n & n0);
            q = Foam::max(q, (n & n1));
            
            if( q > 0.9 )
            {
                createdTriangles.append(t);
                ++counter;
                
                forAll(openEdges, i)
                    problemEdges[openEdges[i]] = false;
            }
        }
        
        if( createdTriangles.size() )
        {
            List<labelledTri> newTriangles = surf.localFaces();
            
            label nTriangles(newTriangles.size());
            newTriangles.setSize(nTriangles+createdTriangles.size());
            forAll(createdTriangles, i)
                newTriangles[nTriangles++] = createdTriangles[i];
            
            triSurface* newts =
                new triSurface
                (
                    newTriangles,
                    surf.patches(),
                    surf.points()
                );

            deleteDemandDrivenData(ts);
            ts = newts;
        }

    } while( createdTriangles.size() );

    ts->write(inFileName, true);
    deleteDemandDrivenData(ts);

    Info << "Added " << counter << " triangles." << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
