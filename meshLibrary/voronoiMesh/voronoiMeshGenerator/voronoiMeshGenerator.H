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
    voronoiMeshGenerator

Description
    Creates voronoi mesh template from the octree

SourceFiles
    voronoiMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef voronoiMeshGenerator_H
#define voronoiMeshGenerator_H

#include "polyMeshGen.H"
#include "IOdictionary.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;
class meshOctree;
class Time;

/*---------------------------------------------------------------------------*\
                     Class voronoiMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class voronoiMeshGenerator
{
    // Private data
        //- reference to Time
        const Time& runTime_;

        //- pointer to the surface
        const triSurf* surfacePtr_;

        //- pointer to the octree
        meshOctree* octreePtr_;

        //- pointer to the list patches for boundary vertices
        labelList* pointRegionsPtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- mesh
        polyMeshGen mesh_;

    // Private member functions
        //- create voronoi mesh
        void createVoronoiMesh();

        //- prepare mesh surface
        void surfacePreparation();

        //- map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- capture edges and corners
        void mapEdgesAndCorners();

        //- optimise surface mesh
        void optimiseMeshSurface();

        //- decompose concave cells
        void decomposeConcaveCells();

        //- add boundary layers
        void generateBoudaryLayers();

        //- mesh optimisation
        void optimiseFinalMesh();

        //- replace boundaries
        void replaceBoundaries();

        //- renumber the mesh
        void renumberMesh();

        //- generate mesh
        void generateMesh();

        //- Disallow default bitwise copy construct
        voronoiMeshGenerator(const voronoiMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const voronoiMeshGenerator&);

public:

    // Constructors

        //- Construct from time
        voronoiMeshGenerator(const Time&);

        //- Construct from time and desired cell size
        //voronoiMeshGenerator(const Time&, const volScalarField&);

    // Destructor

        ~voronoiMeshGenerator();


    // Member Functions

        //- write the mesh
        void writeMesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
