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
    tetPointsCreator

Description
    A class which extracts tet mesh out of a Delaunay structure

SourceFiles
    tetPointsCreator.C

\*---------------------------------------------------------------------------*/

#ifndef tetPointsCreator_H
#define tetPointsCreator_H

#include "bool.H"
#include "scalar.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class tetTessellation;

/*---------------------------------------------------------------------------*\
                           Class tetPointsCreator Declaration
\*---------------------------------------------------------------------------*/

class tetPointsCreator
{
    // Private data
        //- reference to the Delaunay tessellation
        tetTessellation& tessellation_;
    
        //- reference to the octree
        const meshOctree& octree_;

    // Private member functions
        //- create vertices and pointLeaves addressing
        void createPoints();
    
    // Helper functions
        //- checks if the point is inside the surface
        bool isPointInsideSurface(const point&) const;
    
        //- find out if this element has been changed already
        bool isElementChanged(const label elmtI, const label nElements) const;
    
        //- return cell size at the given location
        scalar cellSizeAtLocation(const point&) const;

    // Private copy constructor
        //- Disallow default bitwise copy construct
        tetPointsCreator(const tetPointsCreator&);

        //- Disallow default bitwise assignment
        void operator=(const tetPointsCreator&);

public:

    // Constructors

        //- Construct from octree and mesh data
        tetPointsCreator
        (
            tetTessellation& tessellation,
            const meshOctree& octree
        );

    // Destructor

        ~tetPointsCreator();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
