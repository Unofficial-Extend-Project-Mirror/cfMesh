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
    tetTessellation

Description
    A class for creating tetrahedrisations

SourceFiles
    

\*---------------------------------------------------------------------------*/

#ifndef tetTessellation_H
#define tetTessellation_H

#include "tessellationElement.H"
#include "tessellationDimSpace.H"
#include "point.H"
#include "LongList.H"
#include "boundBox.H"

namespace Foam
{
    
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class tetTessellation Declaration
\*---------------------------------------------------------------------------*/

class tetTessellation
{
protected:
    
    // Protected data
        //- The points
        LongList<point> points_;
    
        //- labels of boundary tet points
        //- boundary tet encompasses the whole structure and all elements of
        //- interest are contained inside it. The labels of these points are
        //- stored such that they can be removed from the final mesh
        label boundaryPoints_[8];

        //- The tessellationElements
        LongList<tessellationElement> elmts_;

        //- helper
        bool pointOk_;

        label delElmts_[1000];
        unsigned short nElmts_;

        point max_;
        point min_;

    // Protected member functions
        label findInitialElement(const point&);

        void treeSearch(const label elmtI, const point&);

        void resetInfluences(const label elmtI);

        void makeNewElements(const label elmtI, const label pI);
        
        void createInitialTets();
        
public:

    // Constructors

        //- Construct from list of points
        tetTessellation(const LongList<point>&);

        //- Construct from triangulated surface
        tetTessellation(const triSurface& surf);

        //- Construct from boundBox
        tetTessellation(const boundBox& bb);

    // Destructor

        virtual ~tetTessellation();

    // Member functions

        //- add point which is already in the list of points
        virtual bool addPoint(const label pI) = 0;

        //- create a new point from the centroid of the given element
        //- and add it into the tessellation
        virtual void addCentroid(const label elmtI) = 0;

        //- create a new point from the circumcentre of the given element
        //- and add it into the tessellation
        virtual void addCircumCentre(const label elmtI) = 0;
        
        //- create a new point at the centre of the given edge of the element
        //- and add it into the tessellation
        virtual void addEdgeCentre(const label elmtI, const direction eI) = 0; 

        const LongList<point>& points() const
        {
            return points_;
        };
        
        const LongList<tessellationElement>& elmts() const
        {
            return elmts_;
        }

        const point& max() const
        {
            return max_;
        }

        const point& min() const
        {
            return min_;
        }

        virtual void checkTessellation() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
