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
    meshSurfaceCutter

Description
    Cut mesh template with the triangulate

SourceFiles
    meshSurfaceCutter.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCutter_H
#define meshSurfaceCutter_H

#include "polyMeshGenModifier.H"
#include "triSurface.H"
#include "boolList.H"
#include "DynList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meshSurfaceCutter Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceCutter
{
    // private data
        //- triangulated surface
        const triSurface& surface_;

        //- mesh data
        polyMeshGen& mesh_;
    
        //- helpers
        label nIntFaces_;
        label nPoints_;
        List<direction> nFacesInCell_;
        boolList problematicTopology_;
        boolList boundaryCell_;
    
        edgeList* edgePointPtr_;
        DynList<point>* newPointsPtr_;
        labelList* newPointLabelPtr_;
        faceListList* facesFromFacePtr_;
        List< DynList<label> > pointTriIndex_;

    // private member functions
        //- return patches for the given point
        labelList patchesForPoint(const label) const;
        
        //- find internal points and intersections with the surface
        void findInternalPointsAndIntersections();

        //- find vertices on boundary edges
        void findBoundaryEdgePoints
        (
            const label,
            const face&,
            const boolList&,
            const boolList&
        );

        //- cuts the face
        bool faceCutter
        (
            const label
        );
        
        void createInternalMeshPointsAndFaces();
        
        void createBoundaryFaces();
        //void checkCellTopology();
        //void checkFaceDirections();

        //- disallows bitwise construct
        void operator=(const meshSurfaceCutter&);

        //- copy constructor
        meshSurfaceCutter(const meshSurfaceCutter&);

public:

    // Constructors

        //- construct components
        meshSurfaceCutter
        (
            polyMeshGen& mesh,
            const triSurface& surface
        );

        //- Destructor
        ~meshSurfaceCutter();

    // Member functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
