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

#include "sortEdgesIntoChains.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
	DynList<edge> bEdges;
	bEdges.append(edge(33, 34));
	bEdges.append(edge(33, 12));
	bEdges.append(edge(34, 12));
	bEdges.append(edge(12, 35));
	bEdges.append(edge(40, 41));
	bEdges.append(edge(40, 12));
	bEdges.append(edge(41, 36));
	bEdges.append(edge(35, 36));
	
	bEdges.append(edge(41, 50));
	bEdges.append(edge(48, 57));
	bEdges.append(edge(50, 57));
	bEdges.append(edge(41, 48));
	
	sortEdgesIntoChains seic(bEdges);
	
	Info << "Sorted chains are " << seic.sortedChains() << endl;
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
