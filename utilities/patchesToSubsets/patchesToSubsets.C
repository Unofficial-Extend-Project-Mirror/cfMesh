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

Description
	Reads the AVL's surface mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "triSurf.H"
#include "triFaceList.H"
#include "labelListPMG.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
	argList::validArgs.append("output subset file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
	fileName subsetFileName(args.args()[2]);
	
	triSurf origSurf(inFileName);
	
    wordList patchNames(origSurf.patches().size());
    forAll(origSurf.patches(), patchI)
        patchNames[patchI] = origSurf.patches()[patchI].name();
    
    forAll(patchNames, patchI)
    {
        labelListPMG subsetFacets;
        forAll(origSurf, triI)
        {
            if( origSurf[triI].region() == patchI )
                subsetFacets.append(triI);
        }
        
        origSurf.addFacetsToSubset(patchNames[patchI], subsetFacets);
    }
	
	origSurf.writeFaceSubsets(subsetFileName);
	
    Info << "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //
