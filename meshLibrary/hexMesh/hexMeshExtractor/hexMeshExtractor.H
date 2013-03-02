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
    hexMeshExtractor

Description
    Extracts hex mesh from the octree
    
SourceFiles
    hexMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef hexMeshExtractor_H
#define hexMeshExtractor_H

#include "polyMeshGenModifier.H"
#include "meshOctreeAddressing.H"
#include "polyMeshGenModifierAddCellByCell.H"
#include "meshOctree.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class hexMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class hexMeshExtractor
{
    // Private data
        //- octree and addressing
        meshOctreeAddressing octreeAddressing_;
    
        //- reference to the mesh
        polyMeshGen& mesh_;
    
        //- centre node labels
        labelListPMG* centreNodeLabelPtr_;
    
        //- subvertices generated inside octree leaves
        VRWGraph* subVerticesPtr_;
    
        //- vertex type and its orientation
        List<direction> octreeVertexType_;

    // Private member functions
        //- delete freestore data
        void clearOut();
    
        //- classify octree points (CORNER, EDGECENTRE, FACECENTRE)
        void classifyOctreePoints();
    
        //- create mesh points and pointLeaves addressing
        void createPoints();

        //- create mesh data
        void createHexMesh();
        
        //- create hexes from FACECENTRE points
        //void createFaceCentreHexes();
        
        //- create hexes which are dual of octree points
        //void createHexesAtOctreePoints();
        
        //- create hexes from octree edges which contains octree cubes
        //- at the same level, and a point at a greater level than the cubes
        //void createEdgeHexes();
        
    // Nested classes
        class createHexCells
        {
            // Private data
                //- mesh modifier
                polyMeshGenModifierAddCellByCell meshModifier_;
            
                //- const reference to the octree
                const meshOctree& octree_;
    
                //- centre node labels
                const labelListPMG& centreNodeLabel_;
    
                //- subvertices generated inside octree leaves
                const VRWGraph& subVertices_;
    
                //- vertex type and its orientation
                const List<direction>& octreeVertexType_;
            
                //- node labels of octree cubes
                const VRWGraph& nodeLabels_;
            
                //- leaves at octree nodes
                const FRWGraph<label, 8>& nodeLeaves_;
            
                //- mapping between Z-order and cylic order of points in face
                FixedList<label, 4> mapping_;
                
            // Private member functions
                //- create hexes from FACECENTRE points
                inline void createFaceCentreHexes(const label pointI);
        
                //- create hexes which are dual of octree points
                inline void createHexesAtOctreePoints(const label pointI);
                
                //- create hexes from octree edges which is contained in
                //- octree cubes at the same level, and a point
                //- at a greater level than the cubes
                inline void createEdgeHexes(const label pointI);
                
                //- create cell faces and add it to the mesh
                inline void addHexCell(const FixedList<label, 8>& cellPoints);
            
            // Private copy constructor
            
                //- Disallow default bitwise copy construct
                createHexCells(const createHexCells&);
        
                //- Disallow default bitwise assignment
                void operator=(const createHexCells&);
        
        public:
        
            // Constructors
                createHexCells
                (
                    polyMeshGen& mesh,
                    const meshOctree& octree,
                    const labelListPMG& centreNodeLabel,
                    const VRWGraph& subVertices,
                    const List<direction>& octreeVertexType,
                    const VRWGraph& nodeLabels,
                    const FRWGraph<label, 8>& nodeLeaves
                );
            
            // Destructor
        
                ~createHexCells();
        
            // Public member functions
                //- create hex cells which are the dual of the octree points
                void generateCells();
        };
    
    // Enumerators
        enum typeOfVertex_
        {
            NONE = 0,
            CORNER = 1,
            NOSUBVERTICES = 2,
            FACECENTRE = 4,
            BOUNDARY = 8,
            MIXED = 16
        };

    // Private copy constructor
        //- Disallow default bitwise copy construct
        hexMeshExtractor(const hexMeshExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const hexMeshExtractor&);

public:

    // Constructors

        //- Construct from octree and mesh data
        hexMeshExtractor
        (
            const meshOctree& octree,
            const IOdictionary& dict,
            polyMeshGen& mesh
        );

    // Destructor

        ~hexMeshExtractor();


    // Member Functions

    void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hexMeshExtractorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
