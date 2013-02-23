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

Application
    Test for creation of boundary faces from chains

Description
    

\*---------------------------------------------------------------------------*/

#include "boundaryFacesGenerator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
	labelList chainPoints(10);
	forAll(chainPoints, pI)
		chainPoints[pI] = pI;
	
	List<DynList<label> > pointPatches(10);
	pointPatches[0].append(1);
	pointPatches[0].append(2);
	
	pointPatches[1].append(2);
	
	pointPatches[2].append(2);
	pointPatches[2].append(3);
	
	pointPatches[3].append(3);
	
	pointPatches[4].append(2);
	pointPatches[4].append(3);
	
	pointPatches[5].append(2);
	
	pointPatches[6].append(0);
	pointPatches[6].append(2);
	
	pointPatches[7].append(0);

	pointPatches[8].append(0);
	pointPatches[8].append(1);
	
	pointPatches[9].append(1);
	
	boundaryFacesGenerator::createFacesFromChain cffc
	(
		chainPoints,
		pointPatches
	);
	
	cffc.createFacesWithoutACorner();
	
	Info << "1. Created faces are " << cffc.createdFaces() << endl;
	Info << "1. Face patches " << cffc.faceRegion() << endl;
	Info << "1. Unresolved points " << cffc.unresolvedPoints() << endl;
	
	cffc.createFacesWithACorner(10);
	
	Info << "2. Created faces are " << cffc.createdFaces() << endl;
	Info << "2. Face patches " << cffc.faceRegion() << endl;
	Info << "2. Unresolved points " << cffc.unresolvedPoints() << endl;	
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
