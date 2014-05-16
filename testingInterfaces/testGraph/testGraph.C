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
    Test for VCWGraph

Description
    - creates an octree and calculates its dual

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "VRWGraph.H"
#include "VRWGraphModifier.H"

#include <omp.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    
//#   include "setRootCase.H"
//#   include "createTime.H"

    VRWGraph origGraph(40000);
    forAll(origGraph, i)
    {
        origGraph.setRowSize(i, 20);

        forAllRow(origGraph, i, j)
            origGraph(i, j) = i+j;
    }
    
    for(label i=0;i<10;++i)
    {
        omp_set_num_threads(1);
        
        Info << "Starting transpose" << endl;
        
        scalar before = omp_get_wtime();
        {
            VRWGraph reverse;
            reverse.reverseAddressing(origGraph);
        
            Info << "Single thread run " << omp_get_wtime() - before << endl;
        }
        
        omp_set_num_threads(2);
        {
            VRWGraph reverse;
            before = omp_get_wtime();
            VRWGraphModifier(reverse).reverseAddressing(origGraph);
        
            Info << "Parallel run " << omp_get_wtime() - before << endl;
            
            VRWGraph rev;
            rev.reverseAddressing(origGraph);
            
            forAll(reverse, rowI)
            {
                if( reverse.sizeOfRow(rowI) != rev.sizeOfRow(rowI) )
                    FatalError << "Drek!" << abort(FatalError);

                forAllRow(reverse, rowI, i)
                    if( reverse(rowI, i) != rev(rowI, i) )
                        FatalError << "Drek 2.!" << abort(FatalError);
            }
        }
    }
    
    //Info << "Orig graph " << origGraph << endl;
    //Info << "Reversed graph " << reverse << endl;
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
