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
    Test for Delaunay code

Description
    

\*---------------------------------------------------------------------------*/

#include "delaunayTessellation.H"
#include "Random.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
	Random ranGen(label(0));
	
	DynListHP<point> points;
	
	scalar scale = 5;
    vector refPoint(10, -20, 0);
	for(label i=0;i<16;++i)
		points.append(scale*(ranGen.vector01() + refPoint));
	
	cpuTime executionTime;
	delaunayTessellation dt(points);
	
	for(label i=0;i<16;++i)
		while( !dt.addPoint(i) )
		{}
	
	Info<< "ExecutionTime = "
        << executionTime.elapsedCpuTime()
		<< " s\n" << endl << endl;
	
	dt.checkTessellation();
	
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
