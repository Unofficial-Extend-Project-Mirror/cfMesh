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
    checkCellConnectionsOverFaces

Description
    Checks cells if all cells can be visited over face neighbours. Deletes

SourceFiles
    checkCellConnectionsOverFaces.C

\*---------------------------------------------------------------------------*/

#ifndef checkCellConnectionsOverFaces_H
#define checkCellConnectionsOverFaces_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class checkCellConnectionsOverFaces Declaration
\*---------------------------------------------------------------------------*/

class checkCellConnectionsOverFaces
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;
    
        //- each cell is assigned a group
        //- cells which can be visited over face neighbours are assigned
        //- to the same group
        labelList cellGroup_;
    
        //- global group label (needed for parallel execution)
        labelList globalGroupLabel_;
    
        //- number of groups
        label nGroups_;

    // Private member functions
        //- decompose marked cells
        void findCellGroups();

        //- Disallow default bitwise copy construct
        checkCellConnectionsOverFaces(const checkCellConnectionsOverFaces&);

        //- Disallow default bitwise assignment
        void operator=(const checkCellConnectionsOverFaces&);

public:

    // Constructors

        //- Construct from polyMeshGen
        checkCellConnectionsOverFaces(polyMeshGen& mesh);

    // Destructor

        ~checkCellConnectionsOverFaces();

    // Member Functions
        //- check if all cells are connected as a single domain
        //- return true if the mesh has been changed
        bool checkCellGroups();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
