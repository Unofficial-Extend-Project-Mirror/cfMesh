/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

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
    //argList::validArgs.append("cellSet");
    
#   include "setRootCase.H"
#   include "createTime.H"
    
    Info<< "Reading mesh for time = " << runTime.value() << endl;
    polyMesh mesh(runTime);

    cellSet cs
    (
        IOobject
        (
            "c0",
            runTime.timeName(),
            "polyMesh/sets",
            runTime,
            IOobject::NO_READ
        )
    );
    
    vector n, origin;
    Info << "Normal x coord" ;
    cin >> n.x();
    Info << endl << "Normal y coord ";
    cin >> n.y();
    Info << endl << "Normal z coord";
    cin >> n.z();
    
    Info << "Origin x coord" ;
    cin >> origin.x();
    Info << endl << "Origin y coord ";
    cin >> origin.y();
    Info << endl << "Origin z coord";
    cin >> origin.z();

    //Info<< "Reading cell set from " << setName << endl << endl;
    //cellSet cs(mesh, setName);

    //- printCells contained in the cellSet
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();
    const labelListList& edgeCells = mesh.edgeCells();
    
    forAll(edges, eI)
    {
        const edge& e = edges[eI];
        bool visibleS, visibleE;
        if( ((points[e.start()] - origin) & n) >= 0.0 )
        {
            visibleS = true;
        }
        else
        {
            visibleS = false;
        }
        
        if( ((points[e.end()] - origin) & n) >= 0.0 )
        {
            visibleE = true;
        }
        else
        {
            visibleE = false;
        }
        
        if( visibleS ^ visibleE )
        {
            forAll(edgeCells[eI], ecI)
                cs.insert(edgeCells[eI][ecI]);
        }
    }
    
    cs.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
