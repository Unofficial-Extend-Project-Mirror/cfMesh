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

\*---------------------------------------------------------------------------*/

#include "boundaryOutwardOrientation.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"

//#define DEBUGMorph

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryOutwardOrientation::checkBoundaryOrientation()
{
    polyMeshGenModifier meshModifier_(mesh_);
    faceListPMG& faces_ = meshModifier_.facesAccess();
    const cellListPMG& cells_ = mesh_.cells();
    const labelList& owner_ = mesh_.owner();
    const labelList& neighbour_ = mesh_.neighbour();
    
    if( mesh_.boundaries().size() == 0 )
    {
        WarningIn
        (
            "void boundaryOutwardOrientation::checkBoundaryOrientation()"
        ) << "Boundary data is not yet calculated!! Cannot proceed!"
            << endl;
    
        return;
    }
    
    const label nIntFaces_ = mesh_.nInternalFaces();
    
    for(register label faceI=nIntFaces_;faceI<faces_.size();faceI++)
    {
        const cell& c = cells_[owner_[faceI]];

        face& newBf = faces_[faceI];
        # ifdef DEBUGMorph
        Info << "Face " << faceI << " is " << newBf <<  endl;
        Info << "Face owner " << owner_[faceI] << endl;
        Info << "Centre " << newBf.centre(points_) << endl;
        Info << "Normal " << newBf.normal(points_) << endl;
        forAll(newBf, pI)
            Info << "Vertex " << pI << " has coordinates "
                << points_[newBf[pI]] << endl;
        # endif

        forAll(c, fI)
            if( neighbour_[c[fI]] != -1 )
            {
                const face& fInt = faces_[c[fI]];
                # ifdef DEBUGMorph
                Info << "Internal face " << fInt << endl;
                Info << "Internal face owner " << owner_[c[fI]] << endl;
                Info << "Internal face neighbour "
                    << neighbour_[c[fI]] << endl;
                Info << "Normal of internal face "
                    << fInt.normal(points_) << endl;
                # endif
                edge e = help::sharedEdge(fInt, newBf);

                if( e.start() != -1 )
                {
                    edge eOther(-1, -1);
                    const edgeList newBfEdges = newBf.edges();
                    forAll(newBfEdges, eJ)
                        if( newBfEdges[eJ] == e )
                            eOther = newBfEdges[eJ];
                
                    if(
                        (
                            (owner_[c[fI]] == owner_[faceI] ) &&
                            (e.start() == eOther.start())
                        ) ||
                        (
                            (neighbour_[c[fI]] == owner_[faceI] ) &&
                            (e.start() == eOther.end())
                        )
                    )
                    {
                        # ifdef DEBUGMorph
                        Info << "Flipping face" << endl;
                        # endif
                        newBf = newBf.reverseFace();
                    }

                    break;
                }
            }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

boundaryOutwardOrientation::boundaryOutwardOrientation
(
    polyMeshGen& mesh
)
    :
    mesh_(mesh)
{
    checkBoundaryOrientation();
}

boundaryOutwardOrientation::~boundaryOutwardOrientation()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
