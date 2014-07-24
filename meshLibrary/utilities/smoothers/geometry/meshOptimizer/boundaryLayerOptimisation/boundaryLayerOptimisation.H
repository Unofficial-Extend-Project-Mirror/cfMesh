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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class boundaryLayerOptimisation Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayerOptimisation
{
    // Private data
        //- reference to polyMeshGen
        polyMeshGen& mesh_;

        //- const reference to meshSurfaceEngine
        const meshSurfaceEngine& meshSurface_;

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

    // Private member functions
        //- calculate hairEdges
        void calculateHairEdges();

        //- Disallow default bitwise copy construct
        boundaryLayerOptimisation(const boundaryLayerOptimisation&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryLayerOptimisation&);

    // Private enumerators

        enum hairEdgeType_
        {
            INSIDE = 1,
            ATEDGE = 2,
            ATCORNER = 4,
            BOUNDARY = 8,
            FEATUREEDGE = 16
        };

public:

    // Constructors

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
        inline const edgeLongList& hairEdges() const
        {
            return hairEdges_;
        }

        //- hair edges attached to a boundary point
        inline const VRWGraph& hairEdgesAtBndPoint() const
        {
            return hairEdgesAtBndPoint_;
        }

        //- boundary faces which serve as base faces of a boundary layer
        //- are set to true
        inline const boolList& isBaseFace() const
        {
            return isBndLayerBase_;
        }

        //- boundary faces where the layers exit at the boundary
        inline const boolList& isExitFace() const
        {
            return isExitFace_;
        }

        //- optimise hair normals to improve quality of layer cells
        void optimiseHairNormals(const label nIterations=10);

        //- optimise thickness variation
        void optimiseThicknessVariation(const scalar tangentTol=0.15);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
