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

#include "tetPointsCreator.H"
#include "tetTessellation.H"
#include "meshOctree.H"
#include "triSurface.H"

//#define DEBUGPoints

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetPointsCreator::createPoints()
{
    //- create initial bunch of points from centroids of the existing tets
    //- repeat this until the first internal point is found
    const LongList<point>& tetPoints = tessellation_.points();
    const LongList<tessellationElement>& elmts = tessellation_.elmts();
    LongList<bool> internalPoint(tetPoints.size(), false);

    bool found;
    label nIter(0);
    do
    {
        Info << "Tessellation has " << tetPoints.size() << " points" << endl;
        Info << "Number of tets " << elmts.size() << endl;
        found = false;
        const label nPoints = tetPoints.size();
        const label nElements = elmts.size();
        for(label elmtI=0;elmtI<nElements;++elmtI)
        //forAll(elmts, elmtI)
        {
            if( isElementChanged(elmtI, nPoints) )
                continue;
            
            const tessellationElement& elmt = elmts[elmtI];
            const point p = elmt.centroid(tetPoints);
            
            if( isPointInsideSurface(p) )
            {
                tessellation_.addCentroid(elmtI);
                found = true;
            }
        }
    } while( found && (++nIter < 7) );
    
    tessellation_.checkTessellation();
}

bool tetPointsCreator::isPointInsideSurface(const point& p) const
{
    const label cLabel = octree_.findLeafContainingVertex(p);
    
    if( cLabel == -1 )
        return false;
    
    if(
        (octree_.returnLeaf(cLabel).cubeType() & meshOctreeCubeBasic::INSIDE) ||
        (octree_.returnLeaf(cLabel).cubeType() & meshOctreeCubeBasic::DATA)
    )
        return true;
    
    return false;
}

bool tetPointsCreator::isElementChanged
(
    const label elmtI,
    const label nElements
) const
{
    if( elmtI >= nElements )
        return true;
    
    const tessellationElement& elmt = tessellation_.elmts()[elmtI];
    
    for(direction i=0;i<DIM;++i)
        if( elmt[i] >= nElements )
            return true;
    
    return false;
}

scalar tetPointsCreator::cellSizeAtLocation(const point& p) const
{
    const label cLabel = octree_.findLeafContainingVertex(p);
    
    if( cLabel == -1 )
        return VGREAT;
    
    return octree_.returnLeaf(cLabel).size(octree_.rootBox());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
