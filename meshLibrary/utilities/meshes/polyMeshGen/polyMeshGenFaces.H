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
    polyMeshGenFaces

Description
    Mesh with geometry and faces. It can be used for volume and surface meshes

SourceFiles
    polyMeshGenFaces.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenFaces_H
#define polyMeshGenFaces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPoints.H"
#include "faceListPMG.H"
#include "labelIOList.H"
#include "writeProcessorPatch.H"
#include "writePatch.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenFaces
    : public polyMeshGenPoints
{
protected:
    
    // Protected data
		//- list of faces
		faceListPMG faces_;
    
        //- inter-processor boundaries. These faces are located after the other
        //- boundary faces in the list of faces. The processor
        //- boundaries are internal faces in the end.
        PtrList<writeProcessorPatch> procBoundaries_;
	
		//- boundary data
		PtrList<writePatch> boundaries_;

        //- face subsets
        std::map<label, meshSubset> faceSubsets_;
	
	// Addressing data
		//- number of internal faces, owner and neighbour
		mutable label nIntFaces_;
		mutable labelIOList* ownerPtr_;
		mutable labelIOList* neighbourPtr_;
		
	// Private member functions
    
        //- calculate owner and neighbour addressing
        virtual void calculateOwnersAndNeighbours() const = 0;
		
		//- clear all pointer data
		void clearOut() const;
	
	// Disallow bitwise assignment
		void operator=(const polyMeshGenFaces&);
	
		polyMeshGenFaces(const polyMeshGenFaces&);

public:
	
	friend class polyMeshGenModifier;

	// Constructors
		//- Null constructor
		polyMeshGenFaces(const Time&);

		//- Construct from components without the boundary
		polyMeshGenFaces
		(
			const Time&,
			const pointField& points,
			const faceList& faces
		);

		//- Construct from components with the boundary
		polyMeshGenFaces
		(
			const Time&,
			const pointField& points,
			const faceList& faces,
			const wordList& patchNames,
			const labelList& patchStart,
			const labelList& nFacesInPatch
		);
		
	// Destructor
		virtual ~polyMeshGenFaces();
		
	// Member functions

		//- access to faces
		inline const faceListPMG& faces() const;
		
		//- return number of internal faces
		inline label nInternalFaces() const;
		
		//- owner and neighbour cells for faces
		inline const labelList& owner() const;
		inline const labelList& neighbour() const;
		
    // Boundary data
    
        //- inter-processor boundaries
        inline const PtrList<writeProcessorPatch>& procBoundaries() const;
        
        //- return processor patch label for the given face label
		label faceIsInProcPatch(const label faceLabel) const;
        
		//- ordinary boundaries
		inline const PtrList<writePatch>& boundaries() const;
		
		//- return patch label for the given face label
		label faceIsInPatch(const label faceLabel) const;
    
    // Subsets
    
        label addFaceSubset(const word&);
        void removeFaceSubset(const label);
        word faceSubsetName(const label) const;
        label faceSubsetIndex(const word&) const;
        inline void addFaceToSubset(const label, const label);
        inline void removeFaceFromSubset(const label, const label);
        inline void faceInSubsets(const label, DynList<label>&) const;
        inline void faceSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void facesInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateFaceSubsets(const ListType&);
        inline void updateFaceSubsets(const VRWGraph&);
		
	// Read mesh
		void read();
		
	// Write mesh
		void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenFacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
