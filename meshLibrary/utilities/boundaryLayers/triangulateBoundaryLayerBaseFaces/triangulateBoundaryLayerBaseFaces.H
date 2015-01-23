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
    triangulateBoundaryLayerBaseFaces

Description
    Splits selected boundary layer cells into triangular prisms

SourceFiles
    triangulateBoundaryLayerBaseFaces.C
    triangulateBoundaryLayerBaseFacesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triangulateBoundaryLayerBaseFaces_H
#define triangulateBoundaryLayerBaseFaces_H

#include "objectRegistry.H"
#include "Time.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "DynList.H"
#include "labelLongList.H"
#include "labelPair.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
            Class triangulateBoundaryLayerBaseFaces Declaration
\*---------------------------------------------------------------------------*/

class triangulateBoundaryLayerBaseFaces
{
    // Private data
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- layers of boundary layer cells belonging to a single column
        const VRWGraph& layerCellsInColumn_;

        //- boolean list containing information about inverted cells
        boolList invertedCell_;

        //- layer column for owner face
        labelList ownerInColumn_;

        //- layer column for neighbour face
        labelList neiInColumn_;

        //- classification of faces in the mesh
        List<direction> faceType_;

        //- label of new point generated at the face centre
        labelList faceCentreLabel_;


    // Private member functions
        //- classify faces in the mesh
        void classifyMeshFaces();

        //- create new points in the mesh
        bool createNewPoints();

        //- create new cells
        void createNewFacesAndCells();

        //- Disallow bitwise copy construct
        triangulateBoundaryLayerBaseFaces
        (
            const triangulateBoundaryLayerBaseFaces&
        );

        //- Disallow bitwise assignment
        void operator=(const triangulateBoundaryLayerBaseFaces&);

    // Enumerators

        enum faceTypes_
        {
            NONE = 0,
            FACEINCOLUMN = 1,
            FACEBETWEENCOLUMNS = 2,
            BASEFACEATBND = 4,
            BASEFACEINSIDEBND = 8,
            BASEFACE = 16
        };

public:

    // Constructors

        //- Construct from mesh
        triangulateBoundaryLayerBaseFaces
        (
            polyMeshGen& mesh,
            const VRWGraph& layerCellsInColumn
        );

    // Destructor
        ~triangulateBoundaryLayerBaseFaces();

    // Public member functions

        // Settings
        //- set invalid cells
        void setBadCells(const labelHashSet&);

        //- set invalid cells
        template<class labelListType>
        void setBadCells(const labelListType& l)
        {
            forAll(l, i)
                invertedCell_[l[i]] = true;
        }

    // Worker
        //- splits inverted boundary layer prisms into triangular prisms
        void triangulateLayers();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
