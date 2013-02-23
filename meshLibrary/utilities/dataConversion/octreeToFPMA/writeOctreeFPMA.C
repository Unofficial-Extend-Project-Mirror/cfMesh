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
    Translates octree to EnSight format

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "meshOctree.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "fileName.H"

#include "octreeMeshFPMA.H"
#include "writeOctreeFPMA.H"

#include <sstream>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOctreeFPMA(const meshOctree& octree, const word& fName)
{
    const word postProcDir = "FPMA";
    
    const word prepend = fName + '.';
	Serr << "Writting mesh into " << postProcDir << endl;

    fileName postProcPath;
	if( Pstream::parRun() )
	{
		std::ostringstream ss;
		ss << Pstream::myProcNo();
		postProcPath = "processor"+ss.str()/postProcDir;
	}
	else
	{
		postProcPath = postProcDir;
	}

	if( !Foam::isDir(postProcPath) )
	{
		mkDir(postProcPath);
	}

	// Open the Case file
	const fileName fpmaFileName = fName + ".fpma";
	
	Info << "Writting octree into " << fpmaFileName << endl;

	OFstream fpmaGeometryFile(postProcPath/fpmaFileName);
	
    // Construct the octree mesh
    dictionary dict;
    octreeMeshFPMA Mesh(octree, dict);
	Mesh.write(fpmaGeometryFile);
}

}

// ************************************************************************* //
