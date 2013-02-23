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

#include "ensightMesh.H"
#include "writeMeshEnsight.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeMeshEnsight(const polyMeshGen& mesh, const word& fName)
{
	const Time& time = mesh.returnTime();
	
    const word postProcDir = "EnSight"+fName;
	const word prepend = fName + '.';
	Info << "Writting mesh into " << postProcDir << endl;

    fileName postProcPath = time.path()/postProcDir;

	if( Foam::isDir(postProcPath) )
	{
		Foam::rmDir(postProcPath);
	}

	Foam::mkDir(postProcPath);

	// Open the Case file
	fileName ensightCaseFileName = prepend + "case";

/*	OFstream ensightCaseFile
	(
		postProcPath/ensightCaseFileName,
		IOstream::ASCII,
		IOstream::currentVersion,
		IOstream::UNCOMPRESSED
	);
*/
	OFstream ensightCaseFile(postProcPath/ensightCaseFileName);

	Info<< nl << "Case file is " << ensightCaseFileName << endl;

/*	OFstream ensightGeometryFile
	(
		postProcPath/prepend+"000.mesh",
		IOstream::ASCII,
		IOstream::currentVersion,
		IOstream::UNCOMPRESSED
	);
*/
	OFstream ensightGeometryFile(postProcPath/prepend+"000.mesh");
	
    // Construct the EnSight mesh
    ensightMesh eMesh(mesh);
	eMesh.write(ensightGeometryFile);

	//- write the case file
    ensightCaseFile << "FORMAT" << nl;
    ensightCaseFile << "type: ensight gold" << nl << nl;

    word geomCaseFileName = prepend + "000";
    if( Pstream::master() )
    {
        // test pre check variable if there is a moving mesh
        ensightCaseFile
            << "GEOMETRY" << nl
            << "model:        1     "
            << (geomCaseFileName + ".mesh").c_str() << nl;
    }
	
    ensightCaseFile << nl << "TIME" << nl
        << "time set:                      " << 1 << nl
        << "number of steps:               " << 1 << nl
        << "filename start number:         " << 0 << nl
        << "filename increment:            " << 1 << nl;

    ensightCaseFile << "time values:" << nl;

    ensightCaseFile.setf(ios_base::scientific, ios_base::floatfield);
    ensightCaseFile.precision(5);
}

}

// ************************************************************************* //
