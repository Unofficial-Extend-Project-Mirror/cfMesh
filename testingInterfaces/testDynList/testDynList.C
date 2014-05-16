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
    Test of the octree dual

Description
    - creates an octree and calculates its dual

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "DynList.H"
#include "point.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    
//#   include "setRootCase.H"
//#   include "createTime.H"
    
    //objectRegistry registry(runTime);
    
    for(label i=0;i<2;++i)
    {
        DynList<label> dlist;
        
        dlist.setSize(6);
        
        dlist = 11;
        
        Info << dlist << endl;
        
        for(label j=0;j<150;j++)
            dlist.append(j);
        
        Info << dlist.size() << " " << dlist.containsAtPosition(88) << endl;
        
        dlist.shrink();
        
        dlist.setSize(10);
        
        Info << "dlist" << dlist << endl;
    }
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
