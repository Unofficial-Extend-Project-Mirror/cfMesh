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
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCutter::checkFaceDirections()
{
    labelList own(polyFaces_.size(),-1);
    labelList nei(nIntFaces_,-1);

    help::calculateOwnersAndNeighbours
    (
        polyFaces_,
        polyCells_,
        own,
        nei
    );

    //- check if the boundary face normals point outside
    for(register label faceI=patchStart_[0];faceI < polyFaces_.size();faceI++)
    {
        const cell& c = polyCells_[own[faceI]];
        
        const edgeList fEdges = polyFaces_[faceI].edges();
        
        forAll(c, fI)
            if( faceI != c[fI] )
            {
                const edgeList nEdges = polyFaces_[c[fI]].edges();
                
                forAll(fEdges, eI)
                {
                    bool found(false);
                    
                    forAll(nEdges, eJ)
                        if( fEdges[eI] == nEdges[eJ] )
                        {
                            if(
                                (
                                    (own[c[fI]] == own[faceI]) &&
                                    (fEdges[eI].start() == nEdges[eJ].start())
                                ) ||
                                (
                                    (nei[c[fI]] == own[faceI]) &&
                                    (fEdges[eI].start() == nEdges[eJ].end())
                                )
                            )
                                polyFaces_[faceI] =
                                    polyFaces_[faceI].reverseFace();
                        }
                        
                    if( found )
                        break;
                }
            }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
