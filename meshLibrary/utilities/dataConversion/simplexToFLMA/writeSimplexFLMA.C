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

Description
    Translates FOAM mesh to EnSight format

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "Time.H"
#include "polyMeshGen.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "fileName.H"

#include "simplexMesh.H"
#include "writeSimplexFLMA.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeSimplexFLMA(const partTetMeshSimplex& simplex, const word& fName)
{
    const word postProcDir = "FLMA";

    fileName postProcPath = postProcDir;

	if( !Foam::isDir(postProcPath) )
	{
		mkDir(postProcPath);
	}

	// Open the flma file
	const fileName flmaFileName = fName + ".flma";
	
	Info << "Writting mesh into " << flmaFileName << endl;

	OFstream flmaGeometryFile(postProcPath/flmaFileName);
	
    // Construct the flma mesh
    simplexMesh Mesh(simplex);
	Mesh.write(flmaGeometryFile);
}

}

// ************************************************************************* //
