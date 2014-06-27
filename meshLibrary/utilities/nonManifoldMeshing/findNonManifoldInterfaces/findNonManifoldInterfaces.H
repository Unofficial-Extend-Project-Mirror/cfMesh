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
    findNonManifoldInterfaces

Description
    Aclass for detecting non-manifold interfaces in the mesh

SourceFiles
    findNonManifoldInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef findNonManifoldInterfaces_H
#define findNonManifoldInterfaces_H

#include "polyMeshGenModifier.H"
#include "VRWGraph.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class findNonManifoldInterfaces Declaration
\*---------------------------------------------------------------------------*/

class findNonManifoldInterfaces
{
        //- reference to mesh
        polyMeshGen& mesh_;

        //- const reference to meshOctree
        const meshOctree& octree_;

        //- graph of mesh cells intesected by surface facets
        VRWGraph patchesIntersectingCells_;

        //- cell groups assigned to neighbouring cells
        //- of the intersected cells
        VRWGraph neiGroups_;

        //- store information about face candidates
        //- and the patches they may be mapped to
        VRWGraph internalFacePatches_;

    // Private member functions
        //- find intersected cells
        void findIntersectedCells();

        //- find neighbouring groups
        void findNeighbouringGroups();

        //- find faces which are candidates to be constrained
        //- to the internal interfaces
        bool findFaceCandidates();

        //- extract non-manifold interfaces
        bool extractInterfaces();

        //- Disallow default bitwise copy construct
        findNonManifoldInterfaces(const findNonManifoldInterfaces&);

        //- Disallow default bitwise assignment
        void operator=(const findNonManifoldInterfaces&);

public:

    // Constructors

        //- Construct from mesh and octree
        findNonManifoldInterfaces(polyMeshGen& mesh, const meshOctree& octree);

    // Destructor

        ~findNonManifoldInterfaces();

    // Public member functions
        //- create interfaces
        void createNonManifoldInterfaces();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
