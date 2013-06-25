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
    hexMeshGenerator

Description
    Generator of hexahedral meshes

SourceFiles
    hexMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef hexMeshGenerator_H
#define hexMeshGenerator_H

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
                      Class hexMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class hexMeshGenerator
{
    // Private data
        //- reference to Time
        const Time& runTime_;

        //- pointer to the surface
        const triSurf* surfacePtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- pointer to the octree
        meshOctree* octreePtr_;

        //- mesh
        polyMeshGen mesh_;

    // Private member functions
        //- generate octree based on user's refinement settings
        void generateOctree();

        //- generate dual mesh
        void generateDualMesh();

        //- prepare mesh surface
        void surfacePreparation();

        //- map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- optimise surface of the mesh
        void optimiseMeshSurface();

        //- add boundary layers
        void generateBoundaryLayers();

        //- mesh optimisation
        void optimiseFinalMesh();

        //- replace boundaries
        void replaceBoundaries();

        //- renumber the mesh
        void renumberMesh();

        //- generate mesh
        void generateMesh();

        //- Disallow default bitwise copy construct
        hexMeshGenerator(const hexMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const hexMeshGenerator&);

public:

    // Constructors

        //- Construct from time
        hexMeshGenerator(const Time&);

        //- Construct from time and desired cell size
        //hexMeshGenerator(const Time&, const volScalarField&);

    // Destructor

        ~hexMeshGenerator();


    // Member Functions

        //- write the mesh
        void writeMesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
