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
    fpmaMesh

Description

SourceFiles
    flmaMesh.C

\*---------------------------------------------------------------------------*/

#ifndef flmaMesh_H
#define flmaMesh_H

#include "polyMeshGen.H"
#include "OFstream.H"
#include "LongList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class flmaMesh Declaration
\*---------------------------------------------------------------------------*/

class flmaMesh
{
    // Private data
    
        //- Reference to the mesh
        const polyMeshGen& mesh_;
    
        //- list of cell types
        labelList cellType_;
    
        //- additional points
        LongList<point> additionalPoints_;
    
        //- centre labels of decomposed faces
        Map<label> faceCentreLabel_;
        
        //- centre labels of decomposed cells
        Map<label> cellCentreLabel_;

    // Private Member Functions
    
        //- Create points for decomposed cells
        void createPointsForDecomposition();

        //- Disallow default bitwise copy construct
        flmaMesh(const flmaMesh&);

        //- Disallow default bitwise assignment
        void operator=(const flmaMesh&);

        void writePoints
        (
            OFstream& flmaGeometryFile
        ) const;
	
		void writeCellTypes
		(
			OFstream& flmaGeometryFile
		) const;

        void writeCells
        (
            OFstream& flmaGeometryFile
        ) const;
		
		void writeSelections
		(
			OFstream& flmaGeometryFile
		) const;

public:

    // Constructors

        //- Construct from polyMeshGen
        flmaMesh(const polyMeshGen& mesh);


    // Destructor

        ~flmaMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            void write(OFstream& flmaGeometryFile) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
