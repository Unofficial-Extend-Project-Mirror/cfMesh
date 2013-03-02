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

#include "subdivisionTessellation.H"
#include "Random.H"
#include "error.H"
#include "DynList.H"
#include "Map.H"
#include "demandDrivenData.H"

//#define DEBUGTessalation

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool subdivisionTessellation::addPoint(const label pI)
{
    FatalErrorIn
    (
        "bool subdivisionTessellation::addPoint(const label)"
    ) << "Not implemented!" << exit(FatalError);
    
    return false;
}
void subdivisionTessellation::addCentroid(const label elmtI)
{
    const label pI = points_.size();
    points_.append(elmts_[elmtI].centroid(points_));
    
    # ifdef DEBUGTessalation
    Info << "subdivisionTessellation:: Adding centroid " << pI << endl;
    # endif
    
    nElmts_ = 0;

    treeSearch(elmtI);

    makeNewElements(elmtI, pI);
    
    resetInfluences(elmtI);
}

void subdivisionTessellation::addCircumCentre(const label elmtI)
{
    FatalErrorIn
    (
        "bool subdivisionTessellation::addPoint(const label)"
    ) << "Not implemented!" << exit(FatalError);
}

void subdivisionTessellation::addEdgeCentre
(
    const label elmtI,
    const direction eI
)
{
    const label pI = points_.size();
    const edge e = elmts_[elmtI].edges()[eI];
    points_.append(0.5*(points_[e[0]]+points_[e[1]]));
    
    # ifdef DEBUGTessalation
    Info << "subdivisionTessellation:: Adding edge centre "
        << points_[pI] << endl;
    # endif
    
    nElmts_ = 0;

    treeSearch(elmtI, e);

    makeNewElementsBisect(elmtI, e, pI);
    
    resetInfluences(elmtI);
}
    
void subdivisionTessellation::treeSearch(const label elmtI)
{
    tessellationElement& elmt = elmts_[elmtI];
    
    //- stop searching it the tessellationElement is marked already
    if( elmt.influence_ & tessellationElement::GOOD ) return;
        
    # ifdef DEBUGTessalation
    Info << "Setting GOOD at element " << elmtI << endl;
    # endif
    
    elmt.influence_ |= tessellationElement::GOOD;
    delElmts_[nElmts_++] = elmtI;
    
    for(direction ineigh=0;ineigh<DIM1;++ineigh)
    {
        const label el = elmt.neighbour(ineigh);
        
        if( el != -1 )
            elmts_[el].influence_ |= tessellationElement::BOUND;
    }
}

void subdivisionTessellation::treeSearch
(
    const label elmtI,
    const edge& e
)
{
    tessellationElement& elmt = elmts_[elmtI];
    
    //- stop searching if there exist an element with unknown influence
    if( !pointOk_ ) return;

    //- stop searching it the tessellationElement is marked already
    if( elmt.influence_ & tessellationElement::GOOD ) return;
        
    if( (elmt.whichPosition(e[0]) != -1) && (elmt.whichPosition(e[1]) != -1) )
    {
        elmt.influence_ |= tessellationElement::GOOD;
        delElmts_[nElmts_++] = elmtI;
        # ifdef DEBUGTessalation
        Info << "Selecting element " << elmtI << "with nodes " << elmt << endl;
        # endif
    }
    else
    {
        elmt.influence_ |= tessellationElement::BOUND;
        return;
    }

    for(direction i=0;i<DIM1;++i)
    {
        const label el = elmt.neighbour(i);
        
        if( el != -1 )
        {
            //- search other tessellationElements
            treeSearch
            (
                el,
                e
            );
        }
    }
}

void subdivisionTessellation::makeNewElementsBisect
(
    const label elmtI,
    const edge& e,
    const label pI
)
{
    # ifdef DEBUGTessalation
    for(label i=0;i<nElmts_;++i)
        Info << "Element to delete is " << elmts_[delElmts_[i]] << endl;
    # endif
    
    //- create new elements
    DynList<tessellationElement>* newElementsPtr =
        new DynList<tessellationElement>(4*nElmts_);
    for(label i=0;i<nElmts_;++i)
    {
        const tessellationElement& elmt = elmts_[delElmts_[i]];
        
        for(direction j=0;j<DIM1;++j)
        {
            const label nei = elmt.neighbour(j);
            if(
                (nei == -1) ||
                (elmts_[nei].influence_ & tessellationElement::BOUND)
            )
            {
                triFace f = elmt.face(j);
                direction nFound(0);
                forAll(f, pJ)
                    if( e.otherVertex(f[pJ]) != -1 )
                        ++nFound;
                if( nFound != 2 )
                {
                    tessellationElement nelmt(f[0], f[1], f[2], pI);
                    nelmt.setNeighbour(3, nei);
                    nelmt.influence_ = elmt.influence_;
                    newElementsPtr->append(nelmt);
                }
            }
        }
    }
    
    //- create labels of new elements and store them
    labelList newLabels(newElementsPtr->size());
    forAll(newLabels, lI)
        if( lI < nElmts_ )
        {
            elmts_[delElmts_[lI]] = (*newElementsPtr)[lI];
            newLabels[lI] = delElmts_[lI];
        }
        else
        {
            newLabels[lI] = elmts_.size();
            elmts_.append((*newElementsPtr)[lI]);
        }
        
    deleteDemandDrivenData(newElementsPtr);
    # ifdef DEBUGTessalation
    Info << "Labels of new elements " << newLabels << endl;
    # endif
    //- update neighbouring information for BOUND elements
    Map<label> newPointLabel;
    forAll(newLabels, lI)
    {
        const tessellationElement& elmt = elmts_[newLabels[lI]];
        for(direction i=0;i<DIM;++i)
            if( !newPointLabel.found(elmt[i]) )
            {
                const label n = newPointLabel.size();
                newPointLabel.insert(elmt[i], n);
            }
        
        const label nei = elmt.neighbour(3);
        if( nei != -1 )
        {
            tessellationElement& nelmt = elmts_[nei];
            for(direction i=0;i<DIM1;++i)
                if( elmt.whichPosition(nelmt[i]) == -1 )
                {
                    nelmt.setNeighbour(i, newLabels[lI]);
                    break;
                }
        }
    }
    
    //- create neighbours of newly created elements
    List< DynList<label> > nodeElements
    (
        newPointLabel.size(),
        DynList<label>(6)
    );
    
    forAll(newLabels, lI)
    {
        const tessellationElement& elmt = elmts_[newLabels[lI]];
        for(direction i=0;i<DIM;++i)
            nodeElements[newPointLabel[elmt[i]]].append(newLabels[lI]);
    }
        
    # ifdef DEBUGTessalation
    Info << "Node elements " << nodeElements << endl;
    Info <<"New point label " << newPointLabel << endl;
    # endif
        
    forAll(newLabels, lI)
    {
        tessellationElement& elmt = elmts_[newLabels[lI]];
        for(direction i=0;i<DIM;++i)
        {
            const label s = newPointLabel[elmt[(i+1)%3]];
            const label e = elmt[(i+2)%3];
            const DynList<label>& nel = nodeElements[s];
            forAll(nel, elI)
                if( (elmts_[nel[elI]].whichPosition(e) != -1) &&
                    (nel[elI] != newLabels[lI])
                )
                {
                    elmt.setNeighbour(i, nel[elI]);
                    break;
                }
        }
    }
    
    # ifdef DEBUGTessalation
    forAll(newLabels, lI)
    {
        const tessellationElement& elmt = elmts_[newLabels[lI]];
        Info << "New element " << newLabels[lI] << " is " << elmt << endl;
        for(direction i=0;i<DIM1;++i)
        {
            Info << "Neighbour over face " << elmt.face(i) << " is "
                << elmt.neighbour(i) << endl;
                        
            triFace f = elmt.face(i);
            if( elmt.neighbour(i) != -1 )
            {
                const label nei = elmt.neighbour(i);
                bool found(false);
                for(direction j=0;j<DIM1;++j)
                {
                    triFace fnei = elmts_[nei].face(j);
                    if( fnei == f )
                    {
                        found = true;
                        if( elmts_[nei].neighbour(j) != newLabels[lI] )
                            FatalError << "Wrong neighbour addressing"
                                << abort(FatalError);
                    }
                }
                
                if( !found )
                    FatalError << "Cannot find neighbour!" << abort(FatalError);
            }
        }
    }
    # endif
}

void subdivisionTessellation::checkTessellation() const
{
    forAll(elmts_, elI)
    {
        const tessellationElement& elmt = elmts_[elI];
        
        if( elmt.mag(points_) < 0.0 )
        {
            FatalErrorIn
            (
                "subdivisionTessellation::checkTessalation()"
            ) << "Element " << elI << " is inverted!" << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
