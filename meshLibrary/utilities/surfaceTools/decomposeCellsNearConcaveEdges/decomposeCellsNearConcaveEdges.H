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
