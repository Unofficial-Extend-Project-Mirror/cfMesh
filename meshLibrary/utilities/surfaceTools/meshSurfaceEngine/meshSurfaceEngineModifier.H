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
    meshSurfaceEngineModifier

Description
    Modifier for the meshSurfaceEngine

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngineModifier_H
#define meshSurfaceEngineModifier_H

#include "meshSurfaceEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
					Class meshSurfaceEngineModifier Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngineModifier
{
    // Private data
        //- reference to the meshSurfaceEngine
		meshSurfaceEngine& surfaceEngine_;

        //- Disallow default bitwise copy construct
        meshSurfaceEngineModifier(const meshSurfaceEngineModifier&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngineModifier&);

public:
	
    // Constructors

        //- Construct from reference to meshSurfaceEngine
        meshSurfaceEngineModifier(meshSurfaceEngine& surfaceEngine);

		//- Construct from const reference to meshSurfaceEngine
		meshSurfaceEngineModifier(const meshSurfaceEngine& surfaceEngine);

    // Destructor

        ~meshSurfaceEngineModifier();

    // Member Functions
		
		//- relocate the selected boundary vertex
		void moveBoundaryVertexNoUpdate(const label bpI, const point& newP);
		
        //- relocate the selected boundary vertex and update geometry data
		void moveBoundaryVertex(const label bpI, const point& newP);
        
        //- update normals of boundary vertices at processor boundaries
        inline void updateVertexNormals()
        {
            if( surfaceEngine_.pointNormalsPtr_ )
                surfaceEngine_.updatePointNormalsAtProcBoundaries();
        }
        
        //- updates faceCentres, face normal, and point normals
        //- of the boundary faces attached to the selected boundary points
        void updateGeometry(const labelListPMG&);
        void updateGeometry();
        
        //- makes sure that all surface vertices at parallel boundaries
        //- have the same coordinates
        void syncVerticesAtParallelBoundaries();
        void syncVerticesAtParallelBoundaries(const labelListPMG&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
