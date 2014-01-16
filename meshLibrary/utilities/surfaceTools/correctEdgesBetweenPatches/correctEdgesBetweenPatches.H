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
    correctEdgesBetweenPatches

Description
    Changes mesh surface until every boundary face shares only one edge with
    the faces in other boundary patch

SourceFiles
    correctEdgesBetweenPatches.C

\*---------------------------------------------------------------------------*/

#ifndef correctEdgesBetweenPatches_H
#define correctEdgesBetweenPatches_H

#include "polyMeshGenModifier.H"
#include "labelLongList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class correctEdgesBetweenPatches Declaration
\*---------------------------------------------------------------------------*/

class correctEdgesBetweenPatches
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        //- pointer to surface engine
        mutable const meshSurfaceEngine* msePtr_;

        //- boundary data
        wordList patchNames_;
        VRWGraph newBoundaryFaces_;
        labelLongList newBoundaryOwners_;
        labelLongList newBoundaryPatches_;

        //- holds data which cells have to be decomposed
        boolList decomposeCell_;
        bool decompose_;

    // Private member functions

        //- construct and return mesh surface
        const meshSurfaceEngine& meshSurface() const;

        //- delete mesh surface
        void clearMeshSurface();

         //- replace boundary
        void replaceBoundary();

        //- decompose corrected cells
        void decomposeCorrectedCells();

        //- decompose bnd faces for which the feature edges are not connected
        //- over vertices and internal faces which have more than one feature
        //- edge at the boundary
        void decomposeProblematicFaces();

        //- perform decompose remaining faces having more than one feature edge
        void patchCorrection();

        //- Disallow default bitwise copy construct
        correctEdgesBetweenPatches(const correctEdgesBetweenPatches&);

        //- Disallow default bitwise assignment
        void operator=(const correctEdgesBetweenPatches&);

public:

    // Constructors

        //- Construct from mesh data
        correctEdgesBetweenPatches(polyMeshGen& mesh);

    // Destructor

        ~correctEdgesBetweenPatches();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
