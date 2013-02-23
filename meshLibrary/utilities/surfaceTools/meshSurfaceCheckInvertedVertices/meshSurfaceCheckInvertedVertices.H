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
    meshSurfaceCheckInvertedVertices

Description
    Checks if any of the surface is tangled at any of its vertices. This
    is performed by checking the dot product between the vertex normal and
    the normals of the faces attached to the vertex.

SourceFiles
    meshSurfaceCheckInvertedVertices.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCheckInvertedVertices_H
#define meshSurfaceCheckInvertedVertices_H

#include "polyMeshGenModifier.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
            Class meshSurfaceCheckInvertedVertices Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceCheckInvertedVertices
{
    // Private data
        //- shall the mesh surface be deleted or not
        bool deleteSurface_;
    
        //- mesh surface
        const meshSurfaceEngine* surfaceEnginePtr_;
    
        //- active surface points
        const boolList* activePointsPtr_;

        //- set of inverted vertices
        labelHashSet invertedVertices_;
    
    // Private member functions
        //- check vertices by calculating dot products
        void checkVertices();
        
        //- Disallow default bitwise copy construct
        meshSurfaceCheckInvertedVertices
        (
            const meshSurfaceCheckInvertedVertices&
        );

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceCheckInvertedVertices&);

public:

    // Constructors

        //- Construct from meshSurfaceEngine
        meshSurfaceCheckInvertedVertices
        (
            const meshSurfaceEngine&,
            const boolList* activePointsPtr = NULL
        );

        //- Construct from the mesh
        meshSurfaceCheckInvertedVertices(polyMeshGen&);

    // Destructor

        ~meshSurfaceCheckInvertedVertices();

    // Member Functions
        //- return the labels of inverted vertices
        inline const labelHashSet& invertedVertices() const
        {
            return invertedVertices_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
