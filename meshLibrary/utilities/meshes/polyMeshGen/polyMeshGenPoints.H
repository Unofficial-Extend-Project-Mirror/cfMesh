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
    polyMeshGenPoints

Description
    Mesh with geometry information

SourceFiles
    polyMeshGenPoints.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenPoints_H
#define polyMeshGenPoints_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Time.H"
#include "meshSubset.H"
#include "pointFieldPMG.H"
#include "DynList.H"

#include <map>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenPoints
{
protected:
    
    // Protected data
		//- reference to the Time registry
		const Time& runTime_;

		//- list of vertices
		pointFieldPMG points_;

        //- map of point subsets
        std::map<label, meshSubset> pointSubsets_;
	
	// Disallow bitwise assignment
		void operator=(const polyMeshGenPoints&);
	
		polyMeshGenPoints(const polyMeshGenPoints&);

public:
	
	friend class polyMeshGenModifier;

	// Constructors
		//- Null constructor
		polyMeshGenPoints(const Time&);

		//- Construct from components without the boundary
		polyMeshGenPoints
		(
			const Time&,
			const pointField& points
		);
		
	// Destructor
		~polyMeshGenPoints();
		
	// Member functions
		//- access to Time
		inline const Time& returnTime() const;
		
		//- access to points
		inline const pointFieldPMG& points() const;
        
        //- non-const access to points
        inline pointFieldPMG& points();
        
        //- append a vertex to the end of the list
        inline void appendVertex(const point& p);
        
        //- point subsets
        label addPointSubset(const word&);
        void removePointSubset(const label);
        word pointSubsetName(const label) const;
        label pointSubsetIndex(const word&) const;
        inline void addPointToSubset(const label, const label);
        inline void removePointFromSubset(const label, const label);
        inline void pointInSubsets(const label, DynList<label>&) const;
        inline void pointSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void pointsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updatePointSubsets(const ListType&);
		
	// Read mesh
		void read();
		
	// Write mesh
		void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
