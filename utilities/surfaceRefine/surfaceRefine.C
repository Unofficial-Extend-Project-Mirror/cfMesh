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
	Refinement of triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

#include "faceTessalation.H"
#include "delaunayPoint.H"

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
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    triSurface* ts = new triSurface(inFileName);

    //- refine the surface where needed
    const pointField& points = ts->points();

    SLList<delaunayPoint> dp;

    forAll(points, pI)
        dp.append(delaunayPoint(points[pI],pI));
    
    faceTessalation ft(dp);

    const IDLList<elementFace>& ef = ft.faces();

    List<labelledTri> tessalationFaces(ft.faces().size());
    label faceI(0);

    for(IDLList<elementFace>::const_iterator fIter = ef.begin();
        fIter != ef.end();
        ++fIter
    )
    {
        bool store(true);
        
        for(direction i=0;i<3;i++)
            if(  fIter().pointI(i)->number() < 0 )
                store = false;
        
        if( store )
        {
            labelledTri ltri
            (
                fIter().pointI(0)->number(),
                fIter().pointI(1)->number(),
                fIter().pointI(2)->number(),
                0
            );
                
            tessalationFaces[faceI++] = ltri;
        }
    }

    tessalationFaces.setSize(faceI);

    const List<labelledTri>& faces = ts->localFaces();
    boolList foundFace(faces.size(), false);
    const labelListList& pointFaces = ts->pointFaces();
    
    forAll(tessalationFaces, tfI)
    {
        const labelledTri& tf = tessalationFaces[tfI];
        const labelList& pf = pointFaces[tessalationFaces[tfI][0]];

        forAll(pf, pfI)
        {
            const labelledTri& ltri = faces[pf[pfI]];

            bool found(true);

            forAll(tf, pI)
            {
                bool foundVrt(false);
                forAll(ltri, pJ)
                    if( tf[pI] == ltri[pJ] )
                    {
                        foundVrt = true;
                        break;
                    }

                if( !foundVrt ) found = false;
            }

            if( found )
            {
                foundFace[pf[pfI]] = true;
                break;
            }
        }
    }
    
    Info << "FoundFace " << foundFace << endl;

    List<labelledTri> newTriangles(faceI);
    faceI = 0;
    pointField newPoints(points);
    label pointI = points.size();

    forAll(foundFace, fI)
        if( foundFace[fI] )
        {
            newTriangles[faceI++] = faces[fI];
        }
        else
        {
            List<labelledTri> nf(3);
            const edgeList edg = faces[fI].edges();

            forAll(edg, eI)
            {
                const edge& e = edg[eI];
                nf[eI] =
                    labelledTri
                    (
                        e.start(),
                        e.end(),
                        pointI,
                        faces[fI].region()
                    );
            }

            forAll(nf, nfI)
                newTriangles[faceI++] = nf[nfI];

            const point p = faces[fI].centre(points);
            newPoints.newElmt(pointI++) = p;
        }

    newPoints.setSize(pointI);
    newTriangles.setSize(faceI);

    triSurface newSurf
    (
        newTriangles,
        ts->patches(),
        newPoints
    );
    
    newSurf.write(outFileName, true);

    delete ts;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
