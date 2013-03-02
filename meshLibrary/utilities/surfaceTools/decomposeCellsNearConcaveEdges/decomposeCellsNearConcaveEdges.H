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
    decomposeCellsNearConcaveEdges

Description
    Decomposes cells near concave edges into pyramids

SourceFiles
    decomposeCellsNearConcaveEdges.C

\*---------------------------------------------------------------------------*/

#ifndef decomposeCellsNearConcaveEdges_H
#define decomposeCellsNearConcaveEdges_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class decomposeCellsNearConcaveEdges Declaration
\*---------------------------------------------------------------------------*/

class decomposeCellsNearConcaveEdges
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;
    
        //- holds data which cells have to be decomposed
        boolList decomposeCell_;

    // Private member functions
 
         //- find problematic cells which need to be decomposed
        void findConcaveFacesAndCells();
        
        //- perform patch correction
        void decomposeConcaveCells();

        //- Disallow default bitwise copy construct
        decomposeCellsNearConcaveEdges(const decomposeCellsNearConcaveEdges&);

        //- Disallow default bitwise assignment
        void operator=(const decomposeCellsNearConcaveEdges&);

public:

    // Constructors

        //- Construct from mesh data
        decomposeCellsNearConcaveEdges(polyMeshGen& mesh);

    // Destructor

        ~decomposeCellsNearConcaveEdges();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
