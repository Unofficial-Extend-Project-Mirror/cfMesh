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

Class
    sortEdgesIntoChains

Description
    Creates boundary faces from a chain of open edges

SourceFiles
    sortEdgesIntoChains.C

\*---------------------------------------------------------------------------*/

#ifndef sortEdgesIntoChains_H
#define sortEdgesIntoChains_H

#include "DynList.H"
#include "labelList.H"
#include "edge.H"
#include "Map.H"
#include "boolList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sortEdgesIntoChains Declaration
\*---------------------------------------------------------------------------*/

class sortEdgesIntoChains
{
    // Members
        const DynList<edge>& bEdges_;
    
        bool openEdges_;
    
        Map<label> newNodeLabel_;
    
        List<DynList<label> > edgesAtPoint_;
    
        DynList<labelList> createdChains_;
    
    // Private member functions
        void createNodeLabels();
    
        bool findPointsBelongingToTheChain
        (
            const label currPos,
            boolList& chainEdges
        ) const;
    
        void shrinkEdges(const boolList& chainEdges);
        
        void createChainFromEdges(const boolList& chainEdges);
        
        void sortEdges();
    
    public:
        
        sortEdgesIntoChains
        (
            const DynList<edge>& bEdges
        );
    
        ~sortEdgesIntoChains();
    
    // Member functions
        
        //- a list of points which have not yet been resolved
        const DynList<labelList>& sortedChains() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
