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
    octreeMesh

Description

SourceFiles
    octreeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef octreeMeshFPMA_H
#define octreeMeshFPMA_H

#include "meshOctree.H"
#include "meshOctreeAddressing.H"
#include "polyMeshGen.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class octreeMeshFPMA Declaration
\*---------------------------------------------------------------------------*/

class octreeMeshFPMA
{
    // Private data
        //- octree addressing
        meshOctreeAddressing addressing_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        octreeMeshFPMA(const octreeMeshFPMA&);

        //- Disallow default bitwise assignment
        void operator=(const octreeMeshFPMA&);

        void writePoints
        (
            OFstream& fpmaGeometryFile
        ) const;
    
        void writeFaces
        (
            OFstream& fpmaGeometryFile
        ) const;

        void writeCells
        (
            OFstream& fpmaGeometryFile
        ) const;
        
        void writeSelections
        (
            OFstream& fpmaGeometryFile
        ) const;
        
public:

    // Constructors

        //- Construct from meshOctree
        octreeMeshFPMA(const meshOctree& octree, const dictionary& dict);

    // Destructor

        ~octreeMeshFPMA();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            void write
            (
                OFstream& fpmaGeometryFile
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
