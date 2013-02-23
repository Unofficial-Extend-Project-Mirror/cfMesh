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
    meshSurfaceEdgeExtractorNonTopo

Description
    Stores boundary faces into patches and captures edges and corners

SourceFiles
    meshSurfaceEdgeExtractorNonTopo.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractorNonTopo_H
#define meshSurfaceEdgeExtractorNonTopo_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractorNonTopo Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractorNonTopo
{
    // Private data
		//- mesh
		polyMeshGen& mesh_;
	
		//- octree
		const meshOctree& meshOctree_;

    // Private member functions
	
		//- distribute boundary faces into patches
		void distributeBoundaryFaces();
	
		//- re-map points after edges have been extracted
		void remapBoundaryPoints();

        //- Disallow default bitwise copy construct
        meshSurfaceEdgeExtractorNonTopo(const meshSurfaceEdgeExtractorNonTopo&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractorNonTopo&);

public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractorNonTopo
        (
			polyMeshGen& mesh,
            const meshOctree& octree
        );

    // Destructor

        ~meshSurfaceEdgeExtractorNonTopo();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
