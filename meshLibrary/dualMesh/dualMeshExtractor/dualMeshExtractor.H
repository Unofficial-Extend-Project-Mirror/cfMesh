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
    dualMeshExtractor

Description
    Octree for mesh generation

SourceFiles
    dualMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef dualMeshExtractor_H
#define dualMeshExtractor_H

#include "polyMeshGenModifier.H"
#include "meshOctreeAddressing.H"
#include "Map.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;

/*---------------------------------------------------------------------------*\
                           Class dualMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class dualMeshExtractor
{
    // Private data
        //- helper for cell creation
        static const direction faceFlip_[6][4];

        //- reference to the octree
        meshOctreeAddressing octreeAddressing_;
    
        //- reference to the mesh
        polyMeshGen& mesh_;
    
        //- centre node labels
        labelLongList* centreNodeLabelPtr_;

    // Private member functions
        //- delete freestore data
        void clearOut();
    
        //- create vertices and pointLeaves addressing
        void createPoints();

        //- create mesh data
        void createPolyMesh();
        
    // Nested classes
        class decomposeCreatedPoly
        {
            // Private data
                const faceList& cellFaces_;
            
                const Map<direction>& nodeLevel_;
            
                boolListList splitFaceNode_;
            
                //- list of edges over which a face has been split
                DynList<edge> splitEdge_;
            
                //- list containing information if a face
                //- is a product of splitting
                DynList<bool> splitFace_;
            
                //- a list of decomposed faces
                DynList<face> decomposedFaces_;
            
                //- helper data
                labelListList dfe_;
                labelListList def_;
                DynList<edge> de_;
                boolList se_;
            
            // Private member functions
                //- decompose non-planar faces into triangles
                void decomposeFaces();
            
                //- calculate helper data
                void calculateHelperAddressing();
                
                //- find additional splitting edges
                void selectAdditionalSplitEdges();
            
                //- distribute new faces into cells
                void selectFacesForCell(List<faceList>& cFaces);
            
                //- create missing faces to topologically close new cells
                void createMissingFaces(List<faceList>& cFaces);
            
            public:
                
                decomposeCreatedPoly
                (
                    const faceList& cf,
                    const Map<direction>& nl
                );
            
                ~decomposeCreatedPoly();
            
                void decomposeCell(List<faceList>& decCells);
        };

    // Private copy constructor
        //- Disallow default bitwise copy construct
        dualMeshExtractor(const dualMeshExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const dualMeshExtractor&);

public:

    // Constructors

        //- Construct from octree and mesh data
        dualMeshExtractor
        (
            const meshOctree& octree,
            const IOdictionary& dict,
            polyMeshGen& mesh
        );

    // Destructor

        ~dualMeshExtractor();


    // Member Functions

    void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
