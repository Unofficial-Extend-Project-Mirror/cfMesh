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

#include "delaunayTessellation.H"
#include "Random.H"
#include "error.H"

//#define DEBUGTessalation

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label tetTessellation::findInitialElement(const point& p)
{
    //- find an tessellationElement with the closest centre to the point
    label el = elmts_.size() - 1;

    //- find an tessellationElement which is influenced by the point
    //- this brutal force search is not to be used for big tessellations
    const label nEl = elmts_.size();
    for(label tI=0;tI<nEl;++tI)
        if( elmts_[tI].influencedBy(points_, p) > VSM )
            return tI;

    return el;
}

void tetTessellation::treeSearch
(
    const label elmtI,
    const point& p
)
{
    tessellationElement& elmt = elmts_[elmtI];
    
    //- stop searching if there exist an element with unknown influence
    if( !pointOk_ ) return;

    //- stop searching it the tessellationElement is marked already
    if( elmt.influence_ & tessellationElement::GOOD ) return;
        
    const scalar infl = elmt.influencedBy(points_, p);

    if( infl > VSM )
    {
        //- set good bit and store element into the list for deletion
        elmt.influence_ |= tessellationElement::GOOD;
        delElmts_[nElmts_++] = elmtI;

        if( nElmts_ >= 65534 )
        {
            pointOk_ = false;
            return;
            FatalErrorIn
            (
                "void subdivisionTessellation::findBouFaces(tessellationElement*)"
            ) << "Exceeded maximum allowed number of deleted elements!"
                << exit(FatalError);
        }
    }
    else if( infl > -VSM )
    {
        point& r = const_cast<point&>(p);
        Random rnd(0);
        
        const point crcm = elmt.crcmCentre(points_);
        point add =
            p -
            1e-7 *
            (
                crcm - p +
                rnd.position(p, crcm)
            );

        r = add;
        pointOk_ = false;
        return;
    }
    else
    {
        elmt.influence_ |= tessellationElement::BOUND;
        return;
    }

    for(direction i=0;i<DIM1;++i)
    {
        label el = elmt.neighbour(i);
        
        if( el != -1 )
        {
            //- search other tessellationElements
            treeSearch
            (
                el,
                p
            );
        }
    }
}

void tetTessellation::resetInfluences(const label elmtI)
{
    tessellationElement& elmt = elmts_[elmtI];
    
    if( !elmt.influence_ ) return;

    elmt.influence_ = tessellationElement::NONE;

    for(direction ineigh=0;ineigh<DIM1;++ineigh)
    {
        label nei = elmt.neighbour(ineigh);

        if( nei != -1 )
            resetInfluences(nei);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
