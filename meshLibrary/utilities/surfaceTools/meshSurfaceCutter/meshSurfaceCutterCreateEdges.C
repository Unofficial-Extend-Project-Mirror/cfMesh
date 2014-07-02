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

\*---------------------------------------------------------------------------*/

#include "meshSurfaceCutter.H"

// #define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void meshSurfaceCutter::createEdges
(
    DynList<edge>& edges,
    labelListList& faceEdges
) const
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    
    //- set sizes of lists
    faceEdges.setSize(faces.size());
    //- create point faces addressing
    labelListList pointFaces(points.size());
    forAll(pointFaces, pI)
        pointFaces[pI].setSize(12);
    List<direction> nAppearances(points.size(), direction(0));

    forAll(faceEdges, fI)
    {
        const face& f = faces[fI];

        forAll(f, pI)
            pointFaces[f[pI]].newElmt(nAppearances[f[pI]]++) = fI;
    }
    
    forAll(nAppearances, pI)
    {
        if( nAppearances[pI] > 6 )
            Warning << "More than 6 faces share a point" << endl;

        pointFaces[pI].setSize(nAppearances[pI]);
    }

    //- create poly edges
    forAll(faceEdges, fI)
    {
        faceEdges[fI].setSize(faces[fI].size());
        faceEdges[fI] = -1;
    }
        
    register label count(0);

    forAll(faceEdges, fI)
    {
        const edgeList fe = faces[fI].edges();

        forAll(fe, eI)
            if( faceEdges[fI][eI] == -1 )
            {
                const edge& e = fe[eI];

                const labelList& efaces = pointFaces[e.start()];
                
                labelList edgeFaces(15);
                labelList edgeInFace(15);
                
                edgeFaces[0] = fI;
                edgeInFace[0] = eI;
                short i(1);

                forAll(efaces, efI)
                    if( efaces[efI] != fI )
                    {
                        const face& nf = faces[efaces[efI]];
                        const edgeList nedges = nf.edges();
                        
                        forAll(nedges, neI)
                            if( nedges[neI] == e )
                            {
                                edgeFaces.newElmt(i) = efaces[efI];
                                edgeInFace.newElmt(i++) = neI;
                                
                                break;
                            }
                    }

                //- store the edge
                edges.append(e);

                for(short j=0;j<i;j++)
                    faceEdges[edgeFaces[j]][edgeInFace[j]] = count;

                count++;
            }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
