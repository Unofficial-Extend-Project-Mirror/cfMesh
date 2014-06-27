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
#include "matrix3D.H"

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
        matrix3D mat;
        mat[0][0] = 1.0;
        mat[0][1] = 0.5;
        mat[0][2] = 0.5;
        mat[1][0] = 0.5;
        mat[1][1] = 1.0;
        mat[1][2] = 0.5;
        mat[2][0] = 0.5;
        mat[2][1] = 0.5;
        mat[2][2] = 1.0;

        Info << "Determinant " << mat.determinant() << endl;

        Info << "Inverse matrix " << mat.inverse() << endl;
        FixedList<scalar, 3> source;
	source[0] = 0;
        source[1] = 1;
        source[2] = 2;

        FixedList<scalar, 3> sol = mat.solve(source);

        Info << "Solution " << sol << endl;
    }
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
