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
    Test of surface morpher

Description
    - morphs mesh surface to make it smooth for mapping

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMeshGenModifier.H"
#include "polyMeshGenAddressing.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
    
    objectRegistry registry(runTime);
    
    polyMeshGen pmg(registry);
    pmg.read();
    
    const labelList& owner = pmg.owner();
    const labelList& neighbour = pmg.neighbour();
    
    label band(0);
    forAll(owner, fI)
    {
        const label diff = neighbour[fI] - owner[fI];
        if( diff > band )
            band = diff;
    }
    
    Info << "Band before renumbering is " << band << endl;
    
    polyMeshGenModifier(pmg).renumberMesh();
    
    const labelList& newOwner = pmg.owner();
    const labelList& newNeighbour = pmg.neighbour();
    band = 0;
    forAll(newOwner, fI)
    {
        const label diff = newNeighbour[fI] - newOwner[fI];
        if( diff > band )
            band = diff;
    }
    
    Info << "Band after renumbering is " << band << endl;
    
    pmg.addressingData().checkMesh(true);
    
    //pmg.write();
    
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
