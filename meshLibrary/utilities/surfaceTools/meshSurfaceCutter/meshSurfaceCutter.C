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

Description

\*---------------------------------------------------------------------------*/

#include "meshSurfaceCutter.H"

//#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
meshSurfaceCutter::meshSurfaceCutter
(
	polyMeshGen& mesh,
	const triSurface& surface
)
:
	surface_(surface),
	mesh_(mesh),
	nIntFaces_(0),
	nPoints_(mesh.points().size()),
	nFacesInCell_(mesh.cells().size(), direction(0)),
	problematicTopology_(mesh.cells().size(), false),
	boundaryCell_(mesh.cells().size(), false),
	edgePointPtr_(NULL),
	newPointsPtr_(NULL),
	newPointLabelPtr_(NULL),
	facesFromFacePtr_(NULL),
	pointTriIndex_()

{
    # ifdef DEBUGCutter
    Info << "Starting creating mesh points and faces " << endl;
    #endif
    createInternalMeshPointsAndFaces();

    # ifdef DEBUGCutter
    Info << "Starting creating boundary faces " << endl;
    #endif
    createBoundaryFaces();

    # ifdef DEBUGCutter
    Info << "Checking for cell consisting of two or more closed parts" << endl;
    # endif
    //checkCellTopology();
	
    # ifdef DEBUGCutter
    Info << "Checking directions of face normals" << endl;
    # endif
    //checkFaceDirections();
}

// Destructor
meshSurfaceCutter::~meshSurfaceCutter()
{
    deleteDemandDrivenData(facesFromFacePtr_);
}
        

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
