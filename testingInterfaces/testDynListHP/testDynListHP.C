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
    Test of the octree dual

Description
    - creates an octree and calculates its dual

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "labelListPMG.H"
#include "pointFieldPMG.H"
#include "DynList.H"
#include "FCWGraph.H"
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
    

    /*
    for(label i=0;i<2;i++)
    {
        DynListHP<point*> l1;

        labelListPMG l2(1000000, 15);
        
        l2.append(13);
        l2.append(18);
        l2.shrink();
        
        for(label i=0;i<10000000;++i)
            l2.append(i);
        
        //Info << "l2 " << l2 << endl;
        
        DynListHP<point> points(100000);
        forAll(points, pI)
            points[pI] = point(pI+1, pI+2, pI+3);
        
        points.append(point(1, 0, 0));
        points.append(point(15, 5, 3));
        points.append(point(2, 3, 1));
        
        FCWGraph<label, 8> lg(10000000);
        forAll(lg, i)
        {
            FixedList<label, 8> lst(10);
            lg.setRow(i, lst);
        }
        
        forAll(points, pI)
            l1.append(&points[pI]);
    }
    
    if( Pstream::parRun() )
    {
        Pout << "Here" << endl;
        labelListPMG data;
        if( Pstream::myProcNo() == 0 )
        {
            data.setSize(100);
            data = 10;
            OPstream toOtherProc(Pstream::nonBlocking, 1, data.byteSize());
            toOtherProc << data;
        }
        else if( Pstream::myProcNo() == 1 )
        {
            labelList d;
            IPstream fromOtherProc(Pstream::nonBlocking, 0);
            fromOtherProc >> d;
            
            Pout << "Received data has size " << d.size() << endl;
            Pout << "Values " << d << endl;
        }
    }
    */
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
