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
    createFacesFromChain

Description
    Creates boundary faces from a chain of open edges

SourceFiles
    createFacesFromChain.C

\*---------------------------------------------------------------------------*/

#ifndef createFacesFromChain_H
#define createFacesFromChain_H

#include "DynList.H"
#include "labelList.H"
#include "face.H"
#include "boolList.H"

namespace Foam
{
	
class VRWGraph;

/*---------------------------------------------------------------------------*\
                           Class createFacesFromChain Declaration
\*---------------------------------------------------------------------------*/

class createFacesFromChain
{
	// Members
		labelList chainPoints_;
	
		List<DynList<label> > regionsForPointAtPosition_;
	
		DynList<face> createdFaces_;
		DynList<label> faceRegions_;
	
	// Private member functions
		void findPointsBelongingToTheFace
		(
			const label currPos,
			boolList& includePoints,
			boolList& endPoints
		) const;
	
		void shrinkTheChain
		(
			const label currPos,
			const boolList& includePoints,
			const boolList& endPoints
		);
	
	public:
		
		createFacesFromChain
		(
			const labelList& chPoints,
			const VRWGraph& pointRegions
		);
	
		createFacesFromChain
		(
			const labelList& chPoints,
			const List<DynList<label> >& pointRegions
		);
	
		~createFacesFromChain();
	
	// Member functions
		//- create faces which do not require an additional corner
		void createFacesWithoutACorner();
	
		//- create faces with a corner
		void createFacesWithACorner(const label);
	
		//- return created faces and the corresponding patches
		const DynList<face>& createdFaces() const;
		const DynList<label>& faceRegion() const;
	
		//- a list of points which have not yet been resolved
		const labelList& unresolvedPoints() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
