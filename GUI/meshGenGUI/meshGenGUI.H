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
    meshGenGUI

Description
    Implementation of the user inteface

SourceFiles
    meshGenGUI.C

\*---------------------------------------------------------------------------*/

#ifndef meshGenGUI_H
#define meshGenGUI_H

#include "IOdictionary.H"
#include "wordList.H"
#include "patchRefinement.H"
#include "objectRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class meshGenGUI Declaration
\*---------------------------------------------------------------------------*/

class meshGenGUI
{
    // Private data

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;
	
		//- pointer to the triangulated surface
		mutable triSurface* surfacePtr_;
	
	// Private member functions
		//- create surface from fileName
		void createSurface(const fileName&) const;

        //- Disallow default bitwise copy construct
        meshGenGUI(const meshGenGUI&);

        //- Disallow default bitwise assignment
        void operator=(const meshGenGUI&);

public:

    // Constructors

        //- Construct from registry
        meshGenGUI(const objectRegistry&);

    // Destructor

        ~meshGenGUI();

    // Member Functions

		//- return surface
		const triSurface& surface() const;
	
	// General meshDict settings
	
		//- set and get surface file name
		void setSurfaceFileName(const fileName&);
		fileName surfaceFileName() const;
		
		//- set and get max cell size
		void setMaxCellSize(const scalar s);
		scalar maxCellSize() const;
		
		//- set and get boundary cell size
		bool boundaryCellSizeEntryExist() const;
		void removeBoundaryCellSize();
		void setBoundaryCellSize(const scalar s);
		scalar boundaryCellSize() const;
		
		//- set and get cell size for automatic refinement
		bool minCellSizeEntryExist() const;
		void removeMinCellSize();
		void setMinCellSize(const scalar s);
		scalar minCellSize() const;
		
		//- set and get if boundary cells should be used
		bool keepCellsIntersectingBoundaryEntryExist() const;
		void setKeepCellsIntersectingBoundary();
		void removeKeepCellsIntersectingBoundary();
		bool keepCellsIntersectingBoundary() const;
		
		//- set and get the check for glued mesh should be used
		bool checkForGluedMeshEntryExist() const;
		void setCheckForGluedMesh();
		void removeCheckForGluedMesh();
		bool checkForGluedMesh() const;
		
	// Settings for patch-wise refinement
		//- set and get patch refinement
		bool patchCellSizeEntryExist() const;
		void addPatchCellSize(const patchRefinement&);
		void removePatchCellSize(const word&);
		List<patchRefinement> patchCellSize() const;
		
	// Settings for keep cells intersecting patches
		//- set end get keep cells intersecting patches
		bool keepCellsIntersectingPatchesEntryExist() const;
		void addKeepCellsIntersectingPatches(const word&);
		void removeKeepCellsIntersectingPatches(const word&);
		wordList keepCellsIntersectingPatches() const;
		
	// Settings for patch-wise boundary layers
		//- set and get patches for which the layers should be created
		bool bndLayersEntryExist() const;
		void addBndLayersPatches(const word&);
		void removeBndLayersPatches(const word&);
		wordList keepBndLayersPatches() const;
	
	// Settings for subset-wise refinement
		//- set and get functions for face-subset refinement
		void setSubsetFileName(const fileName&);
		fileName subsetFileName() const;
		bool subsetCellSizeEntryExist() const;
		void addSubsetCellSize(const patchRefinement&);
		void removeSubsetCellSize(const word&);
		List<patchRefinement> subsetCellSize() const;
		
	// object refinement
		//- general object refinement stuff
		bool objectRefinementEntryExist() const;
		PtrList<entry> objectRefinements() const;
		void addObjectRefinement(const objectRefinement&);
		void removeObjectRefinement(const word& name);
		
	// rename boundaries
		bool renameBoundaryEntryExist() const;
		bool defaultPatchNameEntryExist() const;
		word defaultPatchName() const;
		void setDefaultPatchName(const word&);
		void removeDefaultPatchName();
		bool defaultPatchTypeEntryExist() const;
		word defaultPatchType() const;
		void setDefaultPatchType(const word&);
		void removeDefaultPatchType();
		bool newPatchNamesEntryExist() const;
		PtrList<entry> newPatchNames() const;
		void addNewPatchName(const word&, const dictionary&);
		void removePatchName(const word&);
		
	// write the dictionary
        void writeDict() const
		{
			meshDict_.regIOobject::write();
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
