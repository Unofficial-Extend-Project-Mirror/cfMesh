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

#include "tetCreatorOctree.H"
#include "demandDrivenData.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetCreatorOctree::createTetsFromSplitFaces()
{
    Info << "Creating tets from split faces" << endl;
    
    const labelList& cubeLabel = *cubeLabelPtr_;
    const VRWGraph& subNodeLabels = *subNodeLabelsPtr_;
    const FRWGraph<label, 8>& pLeaves = octreeCheck_.nodeLeaves();
    
    forAll(pLeaves, pointI)
    {
        for(label i=0;i<6;++i)
        {
            const label* fNodes =
                meshOctreeCubeCoordinates::faceNodes_[i];
            
            const label cLabel = pLeaves(pointI, fNodes[0]);
            
            if( cLabel < 0 )
                continue;
            if( cubeLabel[cLabel] < 0 )
                continue;
            
            if(
                (cLabel == pLeaves(pointI, fNodes[1])) &&
                (cLabel == pLeaves(pointI, fNodes[2])) &&
                (cLabel == pLeaves(pointI, fNodes[3]))
            )
            {
                //- create 4 tets
                for(label j=0;j<4;++j)
                {
                    checkAndAppendTet
                    (
                        partTet
                        (
                            pointI,
                            subNodeLabels(cLabel, 7-fNodes[j]),
                            subNodeLabels(cLabel, 7-fNodes[(j+1)%4]),
                            cubeLabel[cLabel]
                        )
                    );
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
