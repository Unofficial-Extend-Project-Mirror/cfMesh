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

#include "error.H"
#include "tetTessellation.H"
#include "Map.H"
#include "demandDrivenData.H"
#include "DynList.H"

//#define DEBUGTessalation

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetTessellation::makeNewElements(const label elmtI, const label pI)
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
                tessellationElement nelmt(f[0], f[1], f[2], pI);
                nelmt.setNeighbour(3, nei);
                nelmt.influence_ = elmt.influence_;
                newElementsPtr->append(nelmt);
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
        DynList<label>()
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
                    elmt.setNeighbour(i, nel[elI]);
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
