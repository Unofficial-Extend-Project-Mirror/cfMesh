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

#include "delaunayTessellation.H"
#include "Random.H"
#include "error.H"

//#define DEBUGTessalation

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool delaunayTessellation::addPoint(const label pI)
{
    # ifdef DEBUGTessalation
    Info << "Adding point " << pI
        << " with coordinates " << points_[pI] << endl;
    # endif
    
    const point& p = points_[pI];
    
    pointOk_ = true;

    nElmts_ = 0;
    
    const label el = findInitialElement(p);

    treeSearch
    (
        el,
        p
    );

    //- skip the point if it produces invalid tessellation
    if( !pointOk_ )
    {
        resetInfluences(el);
        return false;
    }

    makeNewElements(el, pI);
    
    resetInfluences(el);

    return true;
}

void delaunayTessellation::addCentroid(const label elmtI)
{
    const label pI = points_.size();
    points_.append(elmts_[elmtI].centroid(points_));
    const point& p = points_[pI];
    
    do
    {
        pointOk_ = true;
    
        nElmts_ = 0;
    
        treeSearch(elmtI, p);
    
        //- skip the point if it produces invalid tessellation
        if( !pointOk_ )
        {
            resetInfluences(elmtI);
            continue;
        }
    
        makeNewElements(elmtI, pI);
        
        resetInfluences(elmtI);
    }
    while( !pointOk_ );
}

void delaunayTessellation::addCircumCentre(const label elmtI)
{
    const label pI = points_.size();
    points_.append(elmts_[elmtI].crcmCentre(points_));
    const point& p = points_[pI];
    
    do
    {
        pointOk_ = true;
    
        nElmts_ = 0;
    
        treeSearch(elmtI, p);
    
        //- skip the point if it produces invalid tessellation
        if( !pointOk_ )
        {
            resetInfluences(elmtI);
            continue;
        }
    
        makeNewElements(elmtI, pI);
        
        resetInfluences(elmtI);
    }
    while( !pointOk_ );
}

void delaunayTessellation::addEdgeCentre(const label elmtI, const direction eI)
{
    const label pI = points_.size();
    const edge e = elmts_[elmtI].edges()[eI];
    points_.append(0.5*(points_[e[0]]+points_[e[1]]));
    const point& p = points_[pI];
    
    do
    {
        Info << "Adding point " << pI << endl;
        pointOk_ = true;
    
        nElmts_ = 0;
    
        treeSearch(elmtI, p);
    
        //- skip the point if it produces invalid tessellation
        if( !pointOk_ )
        {
            Warning << "Adding point " << pI << " failed" << endl;
            resetInfluences(elmtI);
            continue;
        }
    
        makeNewElements(elmtI, pI);
        
        resetInfluences(elmtI);
    }
    while( !pointOk_ );
}

void delaunayTessellation::checkTessellation() const
{
    forAll(elmts_, elI)
    {
        const tessellationElement& elmt = elmts_[elI];
        
        for(direction i=0;i<DIM1;++i)
            if( elmt.neighbour(i) > elI )
            {
                const tessellationElement& nei = elmts_[elmt.neighbour(i)];
                if( nei.influencedBy(points_, points_[elmt[i]]) > VSM )
                {
                    Info << "Influence "
                        << nei.influencedBy(points_, points_[elmt[i]]) << endl;
                    FatalErrorIn
                    (
                        "delaunayTessellation::checkTessalation()"
                    ) << "This is not a Delaunay hierarchy!!"
                        << abort(FatalError);
                }
            }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
