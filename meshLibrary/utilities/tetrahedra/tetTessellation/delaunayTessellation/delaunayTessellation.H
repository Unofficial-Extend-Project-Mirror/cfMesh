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
    delaunayTesselation

Description
    A class for creating tetrahedrisations

SourceFiles
    

\*---------------------------------------------------------------------------*/

#ifndef delaunayTessellation_H
#define delaunayTessellation_H

#include "tetTessellation.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class delaunayTessellation Declaration
\*---------------------------------------------------------------------------*/

class delaunayTessellation
: public tetTessellation
{
public:

    // Constructors

        //- Construct from list of points
        delaunayTessellation(const LongList<point>&);

        //- Construct from triangulated surface
        delaunayTessellation(const triSurface& surf);

        //- Construct from boundBox
        delaunayTessellation(const boundBox& bb);

    // Destructor

        ~delaunayTessellation();

    // Member functions

        //- add point which is already in the list of points
        bool addPoint(const label pI);

        //- create a new point from the centroid of the given element
        //- and add it into the tessellation
        void addCentroid(const label elmtI);

        //- create a new point from the circumcentre of the given element
        //- and add it into the tessellation
        void addCircumCentre(const label elmtI);
        
        //- create a new point at the centre of the given edge of the element
        //- and add it into the tessellation
        void addEdgeCentre(const label elmtI, const direction eI); 

        //- check if the tesselation obeys Delaunay properties
        void checkTessellation() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
