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
    bondaryLayerOptimisation

Description
    A class intended for improving quality of boundary layer cells. It is intended
    to be used on meshes with one boundary layer which can later be refined.

SourceFiles
    boundaryLayerOptimisation.C
    boundaryLayerOptimisationFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryLayerOptimisation_H
#define boundaryLayerOptimisation_H

#include "DynList.H"
#include "polyMeshGenModifier.H"
#include "labelLongList.H"
#include "boolList.H"
#include "edgeLongList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class dictionary;

/*---------------------------------------------------------------------------*\
                Class boundaryLayerOptimisation Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayerOptimisation
{
    // Private data
        //- reference to polyMeshGen
        polyMeshGen& mesh_;

        //- const pointer to meshSurfaceEngine
        mutable const meshSurfaceEngine* meshSurfacePtr_;
        const bool deleteMeshSurface_;

        //- mesh surface partitioner
        mutable const meshSurfacePartitioner* partitionerPtr_;

        //- boundary layer hairs
        edgeLongList hairEdges_;

        //- hair edges attached to a boundary point
        VRWGraph hairEdgesAtBndPoint_;

        //- hair edge to other hair edges
        VRWGraph hairEdgesNearHairEdge_;

        //- is boundary face a base for a prism in the bnd layer
        boolList isBndLayerBase_;

        //- is boundary face part of a layer where a layer exits
        boolList isExitFace_;

        //- classification of hair edges
        List<direction> hairEdgeType_;

        //- stores information where boundary hairs are made thinner
        boolList thinnedHairEdge_;

    // Private member functions
        //- access to mesh surface
        const meshSurfaceEngine& meshSurface() const;

        //- access to meshSurfacePartitioner
        const meshSurfacePartitioner& surfacePartitioner() const;

        //- calculate normal vector for all patches at a point
        //- points can be filtered by a type of hair edges
        typedef std::map<label, std::pair<point, scalar> > patchNormalType;
        typedef std::map<label, patchNormalType> pointNormalsType;
        void calculateNormalVectors
        (
            const direction eType,
            pointNormalsType&
        ) const;

        //- calculate hairEdges
        void calculateHairEdges();

        //- calculate hair vectors at the boundary
        void calculateHairVectorsAtTheBoundary(vectorField&);

        //- optimise hair normals at the boundary
        //- these are the regions where boundary layers exit the domain
        void optimiseHairNormalsAtTheBoundary(const label nIterations);

        //- optimise hair normals inside the mesh
        void optimiseHairNormalsInside(const label nIterations);

        //- calculate the thickness of the layer for edge heI
        //- respective to the height of edge heJ
        scalar calculateThickness
        (
            const label heI,
            const label heJ,
            const scalar tangentTol,
            const scalar featureSizeFactor
        ) const;

        //- optimise thickness variation
        void optimiseThicknessVariation
        (
            const label nIterations,
            const scalar tangentTol,
            const scalar featureSizeFactor,
            const direction edgeType = (INSIDE|BOUNDARY)
        );

        //- optimise layers at exitting faces due to requests from the inside
        bool optimiseLayersAtExittingFaces();

        //- write vectors into a VTK file. Helper for debugging
        static void writeVTK
        (
            const fileName& fName,
            const pointField& origin,
            const vectorField& vecs
        );

        //- write vector correcposing to hair edges. Helper for debugging
        void writeHairEdges
        (
            const fileName& fName,
            const direction eType,
            const vectorField& vecs
        ) const;

        //- write selected hair edges into a file
        void writeHairEdges(const fileName& fName, const direction eType) const;

        //- Disallow default bitwise copy construct
        boundaryLayerOptimisation(const boundaryLayerOptimisation&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryLayerOptimisation&);

    // Private enumerators

        enum hairEdgeTypes_
        {
            NONE = 0,
            ATEDGE = 1,
            ATCORNER = 2,
            BOUNDARY = 4,
            INSIDE = 8,
            FEATUREEDGE = 16
        };

public:

    // Constructors
        //- construct from polyMeshGen
        boundaryLayerOptimisation(polyMeshGen& mesh);

        //- Construct from polyMeshGen and meshSurfaceEngine
        boundaryLayerOptimisation
        (
            polyMeshGen& mesh,
            const meshSurfaceEngine&
        );


    // Destructor

        ~boundaryLayerOptimisation();

    // Member Functions
        //- return hair edges
        const edgeLongList& hairEdges() const;

        //- hair edges attached to a boundary point
        const VRWGraph& hairEdgesAtBndPoint() const;

        //- boundary faces which serve as base faces of a boundary layer
        //- are set to true
        const boolList& isBaseFace() const;

        //- boundary faces where the layers exit at the boundary
        const boolList& isExitFace() const;

        //- optimise boundary layer
        void optimiseLayer
        (
            const label nIterations = 5,
            const scalar thicknessTol = 0.15,
            const scalar featureSizeFactor = 0.3
        );

    // Static member functions
        //- read the settings from dictionary and run the smoother
        static void optimiseLayer
        (
            const dictionary&,
            boundaryLayerOptimisation&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
