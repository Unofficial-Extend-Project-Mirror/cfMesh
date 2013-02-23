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
    Prints topology of cells in a given cell set

Description
    

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "cellSet.H"
#include "polyMesh.H"
#include "IOobjectList.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("cellSet");
	
#   include "setRootCase.H"
#   include "createTime.H"
	
	Info<< "Reading mesh for time = " << runTime.value() << endl;
	polyMesh mesh(runTime);

    word setName(args.args()[3]);

    Info<< "Reading cell set from " << setName << endl << endl;
	cellSet cs(mesh, setName);

	//- printCells contained in the cellSet
	const cellList& cells = mesh.cells();
	const faceList& faces = mesh.faces();
	
	forAll(cells, cellI)
		if( cs.found(cellI) )
		{
			const cell& c = cells[cellI];
			Info << "Cell " << cellI << " consists of faces " << c << endl;
			
			forAll(c, fI)
				Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;
			
			Info << nl << endl;
		}

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
