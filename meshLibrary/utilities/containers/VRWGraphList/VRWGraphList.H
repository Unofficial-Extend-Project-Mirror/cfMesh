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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    VRWGraphList

Description
    This class is an implementation of a list of graphs
	with variable column width. The implementation is memory efficient.

SourceFiles
    VRWGraphListI.H
    VRWGraphList.C

\*---------------------------------------------------------------------------*/

#ifndef VRWGraphList_H
#define VRWGraphList_H

#include "VRWGraph.H"
#include "subGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VRWGraphList Declaration
\*---------------------------------------------------------------------------*/

class VRWGraphList
{
		
    // Private data
		//- graph containing the data
		VRWGraph data_;
	
		//- number of rows
		LongList<rowElement> rows_;

	// Private member functions
		//- check index
		inline void checkIndex
		(
			const label i,
			const label j,
			const label k
		) const;
	
public:

    // Constructors

        //- Construct null
        inline VRWGraphList();

        //- Copy contructor
		inline VRWGraphList(const VRWGraphList&);

	// Destructor

		inline ~VRWGraphList();

    // Member Functions

        // Access

            //- Returns the number of graphs
            inline label size() const;
			
			//- Returns the number of rows in the graph at that position
			inline label sizeOfGraph(const label posI) const;
			
			//- Return the number of element in the row at the given position
			inline label sizeOfRow(const label posI, const label rowI) const;

            //- Clear the graph
            inline void clear();

    // Member Operators

        //- Append a graph at the end of the graphList
		template<class GraphType>
        inline void appendGraph(const GraphType& l);
		
		//- get and set operators
		inline label operator()
		(
			const label i,
			const label j,
			const label k
		) const;
		
		inline label& operator()(const label i, const label j, const label k);
		
		inline const subGraph<const VRWGraph> operator[](const label i) const;
		
		//- Assignment operator
		inline void operator=(const VRWGraphList&);


    // IOstream operators

        // Write VRWGraphList to Ostream.
        friend Ostream& operator<<(Ostream&, const VRWGraphList&);

        //- Read from Istream, discarding contents of existing VRWGraphList.
/*        friend Istream& operator>> <T, width>
        (
            Istream&,
            VRWGraphList<T, width>&
        );
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VRWGraphListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
