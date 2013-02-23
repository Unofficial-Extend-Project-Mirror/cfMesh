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

\*---------------------------------------------------------------------------*/

#include "ensightMesh.H"
#include "IOmanip.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMeshGen
Foam::ensightMesh::ensightMesh(const polyMeshGen& mesh)
:
    mesh_(mesh)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightMesh::~ensightMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::writePoints
(
    const Foam::scalarField& pointsComponent,
    Foam::OFstream& ensightGeometryFile
) const
{
    forAll(pointsComponent, pointi)
    {
        ensightGeometryFile << setw(12) << pointsComponent[pointi] << nl;
    }
}

void ensightMesh::writePolys
(
    const cellListPMG& cellFaces,
    const faceListPMG& faces,
    OFstream& ensightGeometryFile
) const
{
	ensightGeometryFile << "nfaced" << nl << setw(10) << cellFaces.size() << nl;

	forAll(cellFaces, cI)
	{
		ensightGeometryFile << setw(10) << cellFaces[cI].size() << nl;
	}

	forAll(cellFaces, cI)
	{
		const labelList& cf = cellFaces[cI];

		forAll(cf, faceI)
		{
			ensightGeometryFile  << setw(10) << faces[cf[faceI]].size() << nl;
		}
	}

	forAll(cellFaces, cI)
	{
		const labelList& cf = cellFaces[cI];
        
		forAll(cf, faceI)
		{
			const face& f = faces[cf[faceI]];

			forAll(f, pointI)
			{
				ensightGeometryFile << setw(10) << f[pointI] + 1;
			}
			ensightGeometryFile << nl;
		}
	}
}

void Foam::ensightMesh::writeFaces
(
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    if (patchFaces.size())
    {
		ensightGeometryFile
			<< "nsided" << nl << setw(10) << patchFaces.size() << nl;

		forAll(patchFaces, i)
		{
			ensightGeometryFile
				<< setw(10) << patchFaces[i].size() << nl;
		}

        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];
                
            forAll(patchFace, pointi)
            {
                ensightGeometryFile << setw(10) << patchFace[pointi] + 1;
            }
            ensightGeometryFile << nl;
        }
    }
}

void Foam::ensightMesh::writePatches
(
	Foam::OFstream& ensightGeometryFile
) const
{
	const faceListPMG& faces = mesh_.faces();
	const PtrList<writePatch>& boundaries = mesh_.boundaries();
	
    forAll(boundaries, patchI)
    {
		ensightGeometryFile
			<< "part" << nl
			<< setw(10) << patchI + 2 << nl
			<< boundaries[patchI].patchName() << nl
			<< "coordinates" << nl
			<< setw(10) << mesh_.points().size()
			<< endl;
		
			for (direction d=0; d<vector::nComponents; d++)
			{
				writePoints
				(
					mesh_.points().component(d),
					ensightGeometryFile
				);
			}

			faceList patchFaces(boundaries[patchI].patchSize());
			const label start = boundaries[patchI].patchStart();
			forAll(patchFaces, faceI)
				patchFaces[faceI] = faces[start+faceI];

			writeFaces
			(
				patchFaces,
				ensightGeometryFile
			);
    }
}


void ensightMesh::write
(
    OFstream& ensightGeometryFile
) const
{
    const pointFieldPMG& points = mesh_.points();
    const cellListPMG& cellFaces = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();

	// Set Format
	ensightGeometryFile.setf
	(
		ios_base::scientific,
		ios_base::floatfield
	);
	ensightGeometryFile.precision(5);

	ensightGeometryFile
		<< "FOAM Geometry File " << nl << nl
		<< "node id assign" << nl
		<< "element id assign" << nl;
	
	ensightGeometryFile
		<< "part" << nl
		<< setw(10) << 1 << nl
		<< "FOAM cells" << nl
		<< "coordinates" << nl
		<< setw(10) << points.size()
		<< endl;

	for (direction d=0; d<vector::nComponents; d++)
	{
        scalarField comp = points.component(d)();
		writePoints(comp, ensightGeometryFile);
	}

	writePolys
	(
		cellFaces,
		faces,
		ensightGeometryFile
	);

	writePatches(ensightGeometryFile);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
