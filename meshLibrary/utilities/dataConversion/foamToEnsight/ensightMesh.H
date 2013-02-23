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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ensightMesh

Description

SourceFiles
    ensightMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "polyMeshGen.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightMesh(const ensightMesh&);

        //- Disallow default bitwise assignment
        void operator=(const ensightMesh&);

        void writePoints
        (
            const scalarField& pointsComponent,
            OFstream& ensightGeometryFile
        ) const;

        void writePolys
        (
            const cellListPMG& cellFaces,
            const faceListPMG& faces,
            OFstream& ensightGeometryFile
        ) const;
		
		void writeFaces
		(
			const faceList& patchFaces,
			OFstream& ensightGeometryFile
		) const;
		
		void writePatches
		(
			OFstream& ensightGeometryFile
		) const;

public:

    // Public data

        const polyMeshGen& mesh_;

    // Constructors

        //- Construct from polyMeshGen
        ensightMesh(const polyMeshGen& mesh);


    // Destructor

        ~ensightMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            void write
            (
                OFstream& ensightGeometryFile
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
