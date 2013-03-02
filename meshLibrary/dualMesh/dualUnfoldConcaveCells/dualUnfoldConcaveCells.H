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
    dualUnfoldConcaveCells

Description
    Checks for concave cells near surface edges. Removes vertices and merges
    faces in order to untangle the mesh
    
SourceFiles
    dualUnfoldConcaveCells.C

\*---------------------------------------------------------------------------*/

#ifndef dualUnfoldConcaveCells_H
#define dualUnfoldConcaveCells_H

#include "polyMeshGenModifier.H"
#include "labelListPMG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
class meshSurfaceEngine;
class meshOctree;

/*---------------------------------------------------------------------------*\
                           Class dualUnfoldConcaveCells Declaration
\*---------------------------------------------------------------------------*/

class dualUnfoldConcaveCells
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;
    
        //- reference to the octree
        const meshOctree& octree_;
    
        //- list containing classification of cell type
        List<direction> typeOfCell_;
    
        //- list containing classification of vertex type
        List<direction> typeOfVertex_;
    
        //- new boundary faces and owner cells
        VRWGraph newBoundaryFaces_;
        labelListPMG newBoundaryOwners_;
        labelListPMG newBoundaryPatches_;

    // Private member functions
        //- find edges to remove
        bool findConcaveEdges(const meshSurfaceEngine& mse);
    
        //- mark vertices according to their type
        void markVertexTypes(const meshSurfaceEngine& mse);
        
        //- merge boundary faces of the given cell
        label mergeBoundaryFacesOfCell(const label cellI);
    
        //- merge boundary faces of concave cells
        void storeAndMergeBoundaryFaces(const meshSurfaceEngine& mse);
        
        //- create boundary faces for cells next to concave ones
        void createNeighbouringBoundaryFaces(const meshSurfaceEngine& mse);
        
        //- store remaining boundary faces which have not been altered
        void storeRemainingBoundaryFaces(const meshSurfaceEngine& mse);
    
        //- remove concave vertices from internal faces
        void removeConcaveVerticesFromIntFaces(const meshSurfaceEngine& mse);
    
        //- replace mesh boundary with the newly created boundary faces
        void replaceBoundary();
        
        //- perform a check to determine if all vertices have been removed
        //- correctly, and that all treated cells are topologically closed
        void checkAndRepairBoundary();
        
    // Enumerations
    
        enum vertexType_
        {
            NONE = 0,
            CONCAVE = 1,
            EDGE = 2,
            CORNER = 4,
            MOVED = 8,
            REMOVE = 16
        };
        
        enum cellType_
        {
            INTERNALCELL = 0,
            BOUNDARYCELL = 1,
            CONCAVECELL = 2,
            TREATEDCELL = 4,
            MERGEDCELL = 8
        };

        //- Disallow default bitwise copy construct
        dualUnfoldConcaveCells(const dualUnfoldConcaveCells&);

        //- Disallow default bitwise assignment
        void operator=(const dualUnfoldConcaveCells&);

public:

    // Constructors
        //- Construct from polyMeshGen
        dualUnfoldConcaveCells
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

    // Destructor

        ~dualUnfoldConcaveCells();

    // Member Functions
        //- clean topologically problematic cells
        void unfoldInvalidCells();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
