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
    triangulateNonPlanarBoundaryFaces

Description
    Octree for mesh generation

SourceFiles
    triangulateNonPlanarBoundaryFaces.C

\*---------------------------------------------------------------------------*/

#ifndef triangulateNonPlanarBoundaryFaces_H
#define triangulateNonPlanarBoundaryFaces_H

#include "polyMeshGenModifier.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class triangulateNonPlanarBoundaryFaces Declaration
\*---------------------------------------------------------------------------*/

class triangulateNonPlanarBoundaryFaces
{
    // Private data
        //- mesh modifier
        polyMeshGen& mesh_;

        boolList decomposeCell_;
    
        VRWGraph pointPatches_;

    // Private member functions
        //- find surface vertices
        void findBoundaryCellsToDecompose();

        //- Disallow default bitwise copy construct
        triangulateNonPlanarBoundaryFaces(const triangulateNonPlanarBoundaryFaces&);

        //- Disallow default bitwise assignment
        void operator=(const triangulateNonPlanarBoundaryFaces&);

public:

    // Constructors

        //- Construct from points and cell shapes
        triangulateNonPlanarBoundaryFaces(polyMeshGen& mesh);

    // Destructor

        ~triangulateNonPlanarBoundaryFaces();


    // Member Functions

        void decomposeMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
