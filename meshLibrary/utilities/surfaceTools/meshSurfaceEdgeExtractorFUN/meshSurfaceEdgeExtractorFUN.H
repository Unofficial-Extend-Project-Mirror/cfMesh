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
    meshSurfaceEdgeExtractorFUN

Description
    Stores boundary faces into patches and captures edges and corners
    by inserting fundamental mesh sheets

SourceFiles
    meshSurfaceEdgeExtractorFUN.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractorFUN_H
#define meshSurfaceEdgeExtractorFUN_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractorFUN Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractorFUN
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        //- octree
        const meshOctree& meshOctree_;

        //- mesh surface pointer
        meshSurfaceEngine* surfaceEnginePtr_;

    // Private member functions

        //- return reference to surface engine
        meshSurfaceEngine& surfaceEngine();

        //- clear mesh suirface engine
        void clearOut();

        //- distribute boundary faces into patches
        void distributeBoundaryFaces();

        //- check whether all corners in the surface mesh are present
        //- in the volume mesh
        void reviseCorners();

        //- check whether all edges in the surface mesh are present
        //- in the volume mesh
        void reviseEdges();

        //- find corner vertices and correct patches
        void findCornersAndCorrectPatches();

        //- create fundamental sheets
        void createBasicFundamentalSheets();

        //- smooth the surface
        void smoothMeshSurface();

        //- modify fundamental sheets
        void improveQualityOfFundamentalSheets();

        //- re-map points after edges have been extracted
        void remapBoundaryPoints();

        //- Disallow default construct
        meshSurfaceEdgeExtractorFUN();

        //- Disallow default bitwise copy construct
        meshSurfaceEdgeExtractorFUN(const meshSurfaceEdgeExtractorFUN&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractorFUN&);

public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractorFUN
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

    // Destructor

        ~meshSurfaceEdgeExtractorFUN();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
