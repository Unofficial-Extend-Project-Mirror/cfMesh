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

#include "createFacesFromChain.H"
#include "helperFunctions.H"
#include "VRWGraph.H"
#include "Map.H"

//#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void createFacesFromChain::findPointsBelongingToTheFace
(
    const label currPos,
    boolList& includePoints,
    boolList& endPoints
) const
{
    includePoints.setSize(chainPoints_.size());
    endPoints.setSize(chainPoints_.size());
    endPoints = false;
    includePoints = false;
    
    const label currRegion = regionsForPointAtPosition_[currPos][0];
    if( regionsForPointAtPosition_[currPos].size() != 1 )
    {
        FatalErrorIn
        (
            "void createFacesFromChain::createFacesFromChain::"
            "findVerticesBelongingToTheFace"
        ) << "Trying to create a face from an invalid point!"
            << abort(FatalError);
    }
    
    # ifdef DEBUGCutter
    Info << "Finding vertices belonging to the face from chain "
        << chainPoints_ << endl;
    # endif
    
    forAll(chainPoints_, pI)
    {
        const label pos = (pI+currPos) % chainPoints_.size();
        
        if( regionsForPointAtPosition_[pos].contains(currRegion) )
            includePoints[pos] = true;
        
        if( regionsForPointAtPosition_[pos].size() > 1 )
        {
            endPoints[pos] = true;
            break;
        }
        else if( !includePoints[pos] )
        {
            FatalErrorIn
            (
                "void createFacesFromChain::findVerticesBelongingToTheFace()"
            ) << "Cannot create boundary face"
                << abort(FatalError);
        }
    }
    
    forAllReverse(chainPoints_, pI)
    {
        const label pos = (pI+currPos) % chainPoints_.size();
        if( includePoints[pos] ) break;
        
        if( regionsForPointAtPosition_[pos].contains(currRegion) )
            includePoints[pos] = true;
        
        if( regionsForPointAtPosition_[pos].size() > 1 )
        {
            endPoints[pos] = true;
            break;
        }
        else if( !includePoints[pos] )
        {
            FatalErrorIn
            (
                "void createFacesFromChain::createFacesFromChain::"
                "findVerticesBelongingToTheFace"
            ) << "Cannot create boundary face"
                << abort(FatalError);
        }
    }
    
    # ifdef DEBUGCutter
    Info << "Include points " << includePoints << endl;
    Info << "End points " << endPoints << endl;
    # endif
}

void createFacesFromChain::shrinkTheChain
(
    const label currPos,
    const boolList& includePoints,
    const boolList& endPoints
)
{
    const label currRegion = regionsForPointAtPosition_[currPos][0];
    if( regionsForPointAtPosition_[currPos].size() != 1 )
    {
        FatalErrorIn
        (
            "void createFacesFromChain::createFacesFromChain::"
            "findVerticesBelongingToTheFace"
        ) << "Trying to create a face from an invalid point!"
            << abort(FatalError);
    }
    
    # ifdef DEBUGCutter
    Info << "Shrinking chain " << chainPoints_ << endl;
    # endif
    
    labelList shrinkedChain(chainPoints_.size());
    List<DynList<label> > shrinkedRegions(chainPoints_.size());
    
    direction pI(0);
    
    forAll(chainPoints_, vI)
        if( !includePoints[vI] )
        {
            shrinkedChain[pI] = chainPoints_[vI];
            shrinkedRegions[pI] = regionsForPointAtPosition_[vI];
            ++pI;
        }
        else if( endPoints[vI] )
        {
            shrinkedChain[pI] = chainPoints_[vI];
            forAll(regionsForPointAtPosition_[vI], regI)
                if( regionsForPointAtPosition_[vI][regI] != currRegion )
                    shrinkedRegions[pI].append
                    (
                        regionsForPointAtPosition_[vI][regI]
                    );
            ++pI;
        }
        
    //- set sizes
    shrinkedChain.setSize(pI);
    shrinkedRegions.setSize(pI);
        
    //- store the shrinked lists
    chainPoints_ = shrinkedChain;
    regionsForPointAtPosition_ = shrinkedRegions;
        
    # ifdef DEBUGCutter
    Info << "Shrinked chain " << chainPoints_ << endl;
    Info << "Shrinked regions " << regionsForPointAtPosition_ << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

createFacesFromChain::createFacesFromChain
(
    const labelList& chVertices,
    const VRWGraph& pointRegions
)
:
    chainPoints_(chVertices),
    regionsForPointAtPosition_(chVertices.size()),
    createdFaces_(),
    faceRegions_()
{
    forAll(chVertices, vI)
    {
        const constRow row = pointRegions[chVertices[vI]];
        regionsForPointAtPosition_[vI].clear();
        forAll(row, rI)
            regionsForPointAtPosition_[vI].append(row[rI]);
    }
    
    # ifdef DEBUGCutter
    Info << "Making boundary faces for chain " << chainPoints_ << endl;
    Info << "Regions for chain points " << regionsForPointAtPosition_ << endl;
    # endif
}

createFacesFromChain::createFacesFromChain
(
    const labelList& chVertices,
    const List<DynList<label> >& pointRegions
)
:
    chainPoints_(chVertices),
    regionsForPointAtPosition_(chVertices.size()),
    createdFaces_(),
    faceRegions_()
{
    forAll(chVertices, vI)
    {
        regionsForPointAtPosition_[vI] = pointRegions[chVertices[vI]];
    }
    
    # ifdef DEBUGCutter
    Info << "Making boundary faces for chain " << chainPoints_ << endl;
    Info << "Regions for chain points " << regionsForPointAtPosition_ << endl;
    # endif
}
            
createFacesFromChain::~createFacesFromChain()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Member functions

void createFacesFromChain::createFacesWithoutACorner()
{
    bool found;
    do
    {
        found = false;
        
        boolList includePoints, endPoints;
        
        //- this function removes nodes from the chainPoints_ and creates
        //- new faces. In order to create a valid face, it is important that
        //- that a chain of vertices belonging to a given patch is singly
        //- connected. Non-singly connected chains are not treated until they
        //- become singly connected after elimination of other chains
        
        /*
        boolList usedPoint(chainPoints_.size(), false);
        Map<label> numOfChains;
        forAll(chainPoints_, pI)
            if( !usedPoint[pI] && (regionsForPointAtPosition_[pI].size() == 1) )
            {
                findPointsBelongingToTheFace(pI, includePoints, endPoints);
                
                forAll(includePoints, ipI)
                    if( includePoints[ipI] )
                        usedPoint[ipI] = true;
                    
                if( !numOfChains.found(regionsForPointAtPosition_[pI][0]) )
                {
                    numOfChains.insert
                    (
                        regionsForPointAtPosition_[pI][0],
                        1
                    );
                }
                else
                {
                    numOfChains[regionsForPointAtPosition_[pI][0]]++;
                }
            }
        */
        
        //- start creating faces and eliminating nodes from the chain for
        //- singly connected topologies
        forAll(chainPoints_, pI)
            if(
                (regionsForPointAtPosition_[pI].size() == 1)
            //&& (numOfChains[regionsForPointAtPosition_[pI][0]] == 1)
            )
            {    
                findPointsBelongingToTheFace(pI, includePoints, endPoints);
                
                //- create a new face
                face f(chainPoints_.size());
                direction vrtI(0);
                
                forAll(includePoints, incI)
                    if( includePoints[incI] )
                        f[vrtI++] = chainPoints_[incI];
                    
                DynList<label> facePatches(3);
                forAll(endPoints, epI)
                    if( endPoints[epI] )
                    {
                        const DynList<label>& pr =
                            regionsForPointAtPosition_[epI];
                        
                        forAll(pr, i)
                            facePatches.appendIfNotIn(pr[i]);
                    }
                    
                //- face must contain an additional corner point if the number
                //- of associated patches is greater than 2. Skip creating the
                //- face is this is the case
                if( facePatches.size() > 2 )
                    continue;
                
                found = true;
                    
                if( vrtI > 2 )
                {
                    f.setSize(vrtI);
                    createdFaces_.append(f);
                    faceRegions_.append(regionsForPointAtPosition_[pI][0]);
                
                    # ifdef DEBUGCutter
                    Info << "Created face " << f << endl;
                    # endif
                }

                //- shrink the chainPoints_
                shrinkTheChain(pI, includePoints, endPoints);
                
                break;
            }
    } while( found );
}
            
void createFacesFromChain::createFacesWithACorner
(
    const label cornerLabel
)
{
    label start(-1), facePatch(-1);
    forAll(chainPoints_, cpI)
        if( regionsForPointAtPosition_[cpI].size() == 2 )
        {
            start = cpI;
            
            DynList<label> commonPatches(2);
            const DynList<label>& np =
                regionsForPointAtPosition_[chainPoints_.fcIndex(cpI)];
            
            forAll(np, npI)
                if( regionsForPointAtPosition_[cpI].contains(np[npI]) )
                    commonPatches.append(np[npI]);
                
            if( commonPatches.size() == 1 )
            {
                facePatch = commonPatches[0];
            }
            else
            {
                FatalErrorIn
                (
                    "void createFacesFromChain::createFacesFromChain::"
                    "createFacesWithACorner(const label cornerLabel)"
                ) << "Cannot determine face patch" << abort(FatalError);
            }
            
            break;
        }
        
    //- start creating faces with a corner
    face f(5);
    direction vI(0);
    f[vI++] = chainPoints_[start];
    
    for(label cpI=1;cpI<chainPoints_.size();++cpI)
    {
        const label pos = (cpI + start) % chainPoints_.size();
        
        if( regionsForPointAtPosition_[pos].size() == 1 )
        {
            f.newElmt(vI++) = chainPoints_[pos];
        }
        else if( regionsForPointAtPosition_[pos].size() == 2 )
        {
            //- store old face
            f.newElmt(vI++) = chainPoints_[pos];
            f.newElmt(vI++) = cornerLabel;
            f.setSize(vI);
            vI = 0;
            createdFaces_.append(f);
            faceRegions_.append(facePatch);
            
            //- start creting new face
            f[vI++] = chainPoints_[pos];
            const label ppos =
                regionsForPointAtPosition_[pos].containsAtPosition(facePatch);
            if( ppos == 0 )
            {
                facePatch = regionsForPointAtPosition_[pos][1];
            }
            else
            {
                facePatch = regionsForPointAtPosition_[pos][0];
            }
        }
        else
        {
            FatalErrorIn
            (
                "void createFacesFromChain::createFacesFromChain::"
                "createFacesWithACorner(const label cornerLabel)"
            ) << "Found chain vertex in more than 2 patches!!"
                << abort(FatalError);
        }
    }
    
    //- add the start position into the last face
    f.newElmt(vI++) = chainPoints_[start];
    f.newElmt(vI++) = cornerLabel;
    f.setSize(vI);
    createdFaces_.append(f);
    faceRegions_.append(facePatch);
}

const DynList<face>& createFacesFromChain::createdFaces() const
{
    return createdFaces_;
}

const DynList<label>& createFacesFromChain::faceRegion()const
{
    return faceRegions_;
}

const labelList& createFacesFromChain::unresolvedPoints() const
{
    return chainPoints_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
