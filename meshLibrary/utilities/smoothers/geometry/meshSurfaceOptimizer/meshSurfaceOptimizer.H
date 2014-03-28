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
    meshSurfaceOptimizer

Description
    Surface smoothing without any topological changes

SourceFiles
    meshSurfaceOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceOptimizer_H
#define meshSurfaceOptimizer_H

#include "polyMeshGen.H"
#include "boolList.H"
#include "labelLongList.H"
#include "DynList.H"
#include "triFace.H"
#include "meshSurfaceEngine.H"
#include "parTriFace.H"
#include "meshSurfacePartitioner.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class plane;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceOptimizer Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceOptimizer
{
    // Private data
        //- reference to the mesh surface
        meshSurfaceEngine& surfaceEngine_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- type of surface vertex
        List<direction> vertexType_;

        //- form trihedra
        mutable LongList<triFace>* trianglesPtr_;
        mutable VRWGraph* pointTrianglesPtr_;

        //- surface partitioner
        meshSurfacePartitioner partitioner_;

    // Private member functions
        //- classify surface vertices as PARTITION, EDGE, CORNER
        void classifySurfaceVertices();

        //- calculate surface triangulation
        void calculateTrianglesAndAddressing() const;
        inline const LongList<triFace>& triangles() const;
        inline const VRWGraph& pointTriangles() const;

        //- find vertices which are inverted
        //- mark additional layers of vertices near the inverted ones
        label findInvertedVertices
        (
            boolList& smoothVertex,
            const label nAdditionalLayers = 2
        ) const;

        //- transform into a 2D space in plane
        inline bool transformIntoPlane
        (
            const label bpI,
            const plane& pl,
            vector& vecX,
            vector& vecY,
            DynList<point>& pts,
            DynList<triFace>& trias
        ) const;

        //- transform into a 2D space in plane for parallel boundaries
        inline bool transformIntoPlaneParallel
        (
            const label bpI,
            const plane& pl,
            const std::map<label, DynList<parTriFace> >& m,
            vector& vecX,
            vector& vecY,
            DynList<point>& pts,
            DynList<triFace>& trias
        ) const;

        //- new position of a node after laplacian smoothing
        //- the position is the average of neighbouring vertex positions
        inline point newPositionLaplacian
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- new position of a node after laplacian smoothing
        //- the position is the average of neighbouring face centres
        inline point newPositionLaplacianFC
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- new position of a node after laplacian smoothing
        //- the position is the weighted average of neighbouring face centres
        //- the weights are the magnitudes of neighbouring face area vectors
        inline point newPositionLaplacianWFC
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- new position of a node after using surfaceOptimizer
        inline point newPositionSurfaceOptimizer
        (
            const label bpI,
            const scalar tol = 0.001,
            const bool transformIntoPlane = true
        ) const;

        //- new position of edge vertex
        //- the position is the average of neighbouring edge vertices
        inline point newEdgePositionLaplacian
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- smooth the node using the laplacian smoother
        //- new position is the average of the neighbouring vertices
        void nodeDisplacementLaplacian
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- smooth the node using the laplacian smoother
        //- new position is the average of the centres of faces attached
        //- to the vertex
        void nodeDisplacementLaplacianFC
        (
            const label bpI,
            const bool transformIntoPlane = true
        ) const;

        //- smooth the node using surface optimizer
        void nodeDisplacementSurfaceOptimizer
        (
            const label bpI,
            const scalar tol = 1.0
        );

        //- smooth edge node
        void edgeNodeDisplacement(const label bpI) const;

    // Functions needed for parallel runs

        //- transfer data between processors
        void exchangeData
        (
            const labelLongList& nodesToSmooth,
            std::map<label, DynList<parTriFace> >& m
        ) const;

        //- laplacian smoothing of points at processor boundaries
        //- new position is the average of the neighbouring vertices
        void nodeDisplacementLaplacianParallel
        (
            const labelLongList& nodesToSmooth,
            const bool transformIntoPlane = true
        );

        //- laplacian smoothing of points at processor boundaries
        void nodeDisplacementLaplacianFCParallel
        (
            const labelLongList& nodesToSmooth,
            const bool transformIntoPlane = true
        );

        //- smooth nodes at processor boundaries using surface optimizer
        void nodeDisplacementSurfaceOptimizerParallel
        (
            const labelLongList& nodesToSmooth
        );

        //- smooth edge nodes at processor boundaries
        void edgeNodeDisplacementParallel
        (
            const labelLongList& nodesToSmooth
        );

        //- find invalid faces in the mesh
        label findBadFaces(labelHashSet&, boolList&) const;

        //- Disallow default bitwise copy construct
        meshSurfaceOptimizer(const meshSurfaceOptimizer&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceOptimizer&);

    // Enumerators
        enum vertexTypes
        {
            PARTITION = 1,
            EDGE = 2,
            CORNER = 4,
            PROCBND = 8
        };

public:

    // Constructors

        //- Construct from mesh surface and octree
        meshSurfaceOptimizer
        (
            meshSurfaceEngine& surface,
            const meshOctree& octree
        );


    // Destructor

        ~meshSurfaceOptimizer();

    // Member Functions
        //- runs a surface smoother on the selected boundary points
        bool untangleSurface
        (
            const labelLongList& activeBoundaryPoints,
            const label nAdditionalLayers = 0
        );

        //- checks for invertex surface elements and tries to untangle them
        //- it tries to keep the points on the surface for a couple
        //- of iteration and gives up the final iterations
        //- by default, it smooths two additional layers of elements
        //- around the inverted ones
        bool untangleSurface(const label nAdditionalLayers = 2);

        //- optimize boundary nodes after boundary regions are created
        void optimizeSurface(const label nIterations = 5);

        //- optimize the surface of a 2D mesh
        void optimizeSurface2D(const label nIterations = 5);

        //- untangle the surface of a 2D mesh
        void untangleSurface2D();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSurfaceOptimizerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
