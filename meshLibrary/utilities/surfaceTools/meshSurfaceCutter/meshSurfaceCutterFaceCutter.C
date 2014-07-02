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
#include "polyMeshGenAddressing.H"

#define DEBUG_meshSurfaceCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool meshSurfaceCutter::faceCutter
(
    const label faceI
)
{
    const face& f = mesh_.faces()[faceI];
    bool boundaryCut(false);

    # ifdef DEBUG_meshSurfaceCutter
    Info << nl << "Creating face " << nIntFaces_ << endl;
    Info << "Original face " << f << endl;
    forAll(f, pI)
        Info << "Point " << f[pI] << " " << mesh_.points()[f[pI]] << endl;
    # endif

    const edgeList& edges = mesh_.addressingData().edges();
    const VRWGraph& faceEdges = mesh_.addressingData().faceEdges();
    labelList fEdges(faceEdges.sizeOfRow(faceI));
    forAllRow(faceEdges, faceI, eI)
        fEdges[eI] = faceEdges(faceI, eI);
    
    const labelList& newPointLabels = *newPointLabelPtr_;
    const edgeList& edgePoint = *edgePointPtr_;

    # ifdef DEBUG_meshSurfaceCutter
    const DynList<point>& newPoints = *newPointsPtr_;
    forAll(f, pI)
        Info << "Original edge " << pI << "  " << edges[fEdges[pI]] << endl;
    # endif
    
    face newF(f.size());
    boolList inward(f.size());
    boolList outward(f.size());
    short i(0);

    //- add labels of points which are in the domain
    forAll(f, fpI)
    {
        //- add the points at which the original face cuts the boundary
        label startL(-1), endL(-1);
        const edge& fe = edges[fEdges[fpI]];
        
        const edge& e = edgePoint[fEdges[fpI]];
        
        if( fe.start() == f[fpI] )
        {
            startL = e.start();
            endL = e.end();
        }
        else if( fe.end() == f[fpI] )
        {
            startL = e.end();
            endL = e.start();
        }
        else
        {
            FatalErrorIn
            (
                "face faceCutter()"
            ) << "Wrong edge has been addressed" << abort(FatalError);
        }
        
        if( newPointLabels[f[fpI]] != -1 )
        {
            # ifdef DEBUG_meshSurfaceCutter
            Info << "1. Adding point " << newPointLabels[f[fpI]] << endl;
            Info << "1. Old point " << f[fpI] << " has coordinates "
                << newPoints[newPointLabels[f[fpI]]] << endl;
            # endif
            outward.newElmt(i) = false;
            inward.newElmt(i) = false;
            newF.newElmt(i++) = newPointLabels[f[fpI]];

            if( (startL != -1) && (endL != -1) )
            {
                outward.newElmt(i) = true;
                inward.newElmt(i) = false;
                newF.newElmt(i++) = startL;

                inward.newElmt(i) = true;
                outward.newElmt(i) = false;
                newF.newElmt(i++) = endL;

                boundaryCut = true;
            }
            else if( startL != -1 )
            {
                outward.newElmt(i) = true;
                inward.newElmt(i) = false;
                newF.newElmt(i++) = startL;

                boundaryCut = true;
            }
            else if( endL != -1 )
            {
                 FatalErrorIn
                 (
                     "face faceCutter()"
                 ) << "Wrong edge has been addressed" << abort(FatalError);
            }
        }
        else
        {
            if( (startL != -1) && (endL != -1) )
            {
                inward.newElmt(i) = true;
                outward.newElmt(i) = false;
                newF.newElmt(i++) = startL;

                outward.newElmt(i) = true;
                inward.newElmt(i) = false;
                newF.newElmt(i++) = endL;

                boundaryCut = true;
            }
            else if( endL != -1 )
            {
                inward.newElmt(i) = true;
                outward.newElmt(i) = false;
                newF.newElmt(i++) = endL;

                boundaryCut = true;
            }
            else if( startL != -1 )
            {
                FatalErrorIn
                (
                    "face faceCutter()"
                ) << "Wrong edge has been addressed" << abort(FatalError);
            }
        }
    }

    # ifdef DEBUG_meshSurfaceCutter
    Info << "newF " << newF << endl;
    Info << "outward " << outward << endl;
    Info << "inward " << inward << endl;
    # endif

    //- check for points on the edges of the boundary regions
    if( i > 1 )
    {
        newF.setSize(i);

        # ifdef DEBUG_meshSurfaceCutter
        Info << "Finding edge points for face " << faceI << endl;
        Info << "i " << i << endl;
        Info << "newF " << newF << endl;
        forAll(newF, npI)
            Info << "Point " << newF[npI] << " is "
                << newPoints[newF[npI]] << endl;

        # endif

        findBoundaryEdgePoints
        (
            faceI,
            newF,
            outward,
            inward
        );

        return boundaryCut;
    }
    else if( i == 0 )
    {
        return boundaryCut;
    }
    else
    {
        FatalErrorIn
        (
            "void meshSurfaceCutter::faceCutter()"
        ) << "It will not be possible to create boundary faces!"
            << " Please reduce the cell size and try again"
            << exit(FatalError);
    }

    return boundaryCut;
}        
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
