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
    surfaceMorpherCells

Description
    Changes the boundary of the mesh template, such that each boundary cell
    has only one boundary face.

SourceFiles
    surfaceMorpherCells.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMorpherCells_H
#define surfaceMorpherCells_H

#include "polyMeshGenModifier.H"
#include "labelLongList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceMorpherCells Declaration
\*---------------------------------------------------------------------------*/

class surfaceMorpherCells
{
    // Private data
        //- mesh modifier
        polyMeshGen& mesh_;

        //- number of internal faces
        label nIntFaces_;
    
        //- true for vertices which are at the boundary
        boolList boundaryVertex_;
    
        //- true for cells which are at the boundary
        List<direction> cellFlags_;
    
        //- new boundary faces and owner cells
        VRWGraph newBoundaryFaces_;
        labelLongList newBoundaryOwners_;
        labelLongList newBoundaryPatches_;
    
    // Private enumerators
        enum cellTypes
        {
            NONE = 0,
            BOUNDARY = 1,
            CHANGED = 2
        };

    // Private member functions
        //- find surface vertices
        void findBoundaryVertices();
        
        //- find boundary cells
        void findBoundaryCells();
        
        //- remove cells with all vertices at the boundary
        bool removeCellsWithAllVerticesAtTheBoundary();
        
        //- morph boundary faces
        bool morphBoundaryFaces();
        
        //- morph internal faces
        bool morphInternalFaces();
        
        //- replace the boundary of the mesh with the newly created one
        void replaceMeshBoundary();

        //- Disallow default bitwise copy construct
        surfaceMorpherCells(const surfaceMorpherCells&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceMorpherCells&);

public:

    // Constructors

        //- Construct from points and cell shapes
        surfaceMorpherCells(polyMeshGen& mesh);

    // Destructor

        ~surfaceMorpherCells();


    // Member Functions

        void morphMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
