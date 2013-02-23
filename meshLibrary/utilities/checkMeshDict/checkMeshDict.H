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

Class
    checkMeshDict

Description
    Check whether the meshDict file is set correctly

SourceFiles
    checkMeshDict.C

\*---------------------------------------------------------------------------*/

#ifndef checkMeshDict_H
#define checkMeshDict_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class checkMeshDict Declaration
\*---------------------------------------------------------------------------*/

class checkMeshDict
{
		//- Reference to the mesh
        const IOdictionary& meshDict_;
	
	// Private member functions
        //- check patchCellSize entry
        void checkPatchCellSize();
    
        //- check keepCellsIntersectingPatches entry
        void checkKeepCellsIntersectingPatches();
    
        //- check removeCellsIntersectingPatches entry
        void checkRemoveCellsIntersectingPatches();
    
        //- check objectRefinements entry
        void checkObjectRefinements();
    
        //- check entry for boundary layers
        void checkBoundaryLayers();
    
        //- check renameBoundary entry
        void checkRenameBoundary();
    
        //- perform all checks
        void checkEntries();
	
	public:
		
	// Constructors
	
		//- Construct from IOdictionary
		checkMeshDict(const IOdictionary& meshDict);
	
	// Destructor
		~checkMeshDict();
	
	// Public member functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
