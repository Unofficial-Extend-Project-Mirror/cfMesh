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
    delaunayTesselation

Description
    A class for creating tetrahedrisations

SourceFiles
    

\*---------------------------------------------------------------------------*/

#ifndef subdivisionTessellation_H
#define subdivisionTessellation_H

#include "tetTessellation.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class subdivisionTessellation Declaration
\*---------------------------------------------------------------------------*/

class subdivisionTessellation
: public tetTessellation
{
    // Private member functions
        //- mark all tetrahedra sharing edge e (used in addEdgeCentre)
        void treeSearch(const label elmtI, const edge& e);
    
        //- mark tetrahedron (used in addCentroid);
        void treeSearch(const label elmtI);
    
        //- (used in addEdgeCentre);
        void makeNewElementsBisect
        (
            const label elmtI,
            const edge& e,
            const label pI
        );
        
public:

    // Constructors

        //- Construct from list of points
        subdivisionTessellation(const LongList<point>&);

        //- Construct from triangulated surface
        subdivisionTessellation(const triSurface& surf);

        //- Construct from boundBox
        subdivisionTessellation(const boundBox& bb);

    // Destructor

        ~subdivisionTessellation();

    // Member functions
        //- this function is not implemented
        bool addPoint(const label);

        //- create a new point from the centroid of the given element
        //- and add it into the tessellation
        void addCentroid(const label elmtI);

        //- this is not implemented
        void addCircumCentre(const label);
        
        //- create a new point at the centre of the given edge of the element
        //- and add it into the tessellation
        void addEdgeCentre(const label elmtI, const direction eI); 

        //- check if all elements have positive volume
        void checkTessellation() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
