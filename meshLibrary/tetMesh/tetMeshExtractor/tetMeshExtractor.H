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
    tetMeshExtractor

Description
    A class which extracts tet mesh out of a Delaunay structure

SourceFiles
    tetMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshExtractor_H
#define tetMeshExtractor_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class tetTessellation;

/*---------------------------------------------------------------------------*\
                           Class tetMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class tetMeshExtractor
{
    // Private data
    //- reference to the Delaunay tessellation
    const tetTessellation& tessellation_;

    //- reference to the octree
    const meshOctree& octree_;

    //- reference to the mesh
    polyMeshGen& mesh_;

    //- include DATA boxes
    bool useBoundaryTets_;

    //- list of patches for which data boxes should be used as mesh cells
    labelList patchesForBoundaryTet_;

    //- a list which contains information about which Delaunay elements
    //- will be used
    boolList useElement_;

    // Private member functions
    //- find elements which will be used as mesh cells
    void selectElements();

    //- create vertices and pointLeaves addressing
    void createPoints();

    //- create mesh data
    void createPolyMesh();

    // Private copy constructor
    //- Disallow default bitwise copy construct
    tetMeshExtractor ( const tetMeshExtractor& );

    //- Disallow default bitwise assignment
    void operator= ( const tetMeshExtractor& );

public:

    // Constructors

    //- Construct from octree and mesh data
    tetMeshExtractor
    (
        const tetTessellation& tessellation,
        const meshOctree& octree,
        polyMeshGen& mesh
    );

    // Destructor

    ~tetMeshExtractor();


    // Member Functions
    void useTetsIntersectingBoundary();
    void useTetsIntersectingBoundaryPatches ( const wordList& patchNames );
    void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
