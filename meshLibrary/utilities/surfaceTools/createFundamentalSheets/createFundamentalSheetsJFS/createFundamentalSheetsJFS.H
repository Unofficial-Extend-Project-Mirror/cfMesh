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
    createFundamentalSheetsJFS

Description
    Inserts sheets at the boundary of the mesh to capture all feature edges

SourceFiles
    createFundamentalSheetsJFS.C

\*---------------------------------------------------------------------------*/

#ifndef createFundamentalSheetsJFS_H
#define createFundamentalSheetsJFS_H

#include "polyMeshGenModifier.H"
#include "createFundamentalSheets.H"
#include "labelListPMG.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class createFundamentalSheetsJFS Declaration
\*---------------------------------------------------------------------------*/

class createFundamentalSheetsJFS
:
    public createFundamentalSheets
{
    // Private member functions
        //- create inital sheet from all boundary faces of the surface mesh
        void createInitialSheet();

        //- create fundamental sheets for all feature edges
        void createSheetsAtFeatureEdges();

        //- Disallow default construct
        createFundamentalSheetsJFS();

        //- Disallow default bitwise copy construct
        createFundamentalSheetsJFS(const createFundamentalSheetsJFS&);

        //- Disallow default bitwise assignment
        void operator=(const createFundamentalSheetsJFS&);

public:

    // Runtime type information
        TypeName("Shepherd");

    // Constructors

        //- Construct from mesh data
        createFundamentalSheetsJFS(polyMeshGen& mesh);

    // Destructor

        ~createFundamentalSheetsJFS();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
