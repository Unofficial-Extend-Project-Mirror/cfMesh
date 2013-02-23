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
    meshSubset

Description
    A class containing labels of elements

SourceFiles
    meshSubsetI.H

\*---------------------------------------------------------------------------*/

#ifndef meshSubset_H
#define meshSubset_H

#include "label.H"
#include "HashSet.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;
class VRWGraph;

/*---------------------------------------------------------------------------*\
                        class meshSubset Declaration
\*---------------------------------------------------------------------------*/

class meshSubset
{
	// Private data
        //- name of the given subset
        word name_;
    
        //- type of subset
        label type_;
    
        //- labels of elements
        HashSet<label> data_;

public:
    
    // Enumerators

        enum subsetType_
        {
            CELLSUBSET = 1,
            FACESUBSET = 2,
            POINTSUBSET = 4
        };

    // Constructors

        //- Construct from name. Creates an empty subset
        inline meshSubset(const word& name, const subsetType_&);

        //- Construct from name and element labels
        template<class ListType>
        inline meshSubset
        (
            const word& name,
            const subsetType_&,
            const ListType& elements
        );

	// Destructor
		~meshSubset();


    // Member Functions

        // Access

            //- Return name
            inline const word& name() const;
        
            //- Type of subset
            inline label type() const;

            //- elements contained in the subset
            template<class ListType>
            inline void containedElements(ListType&) const;

        // Edit
        
            //- add element label to subset
            inline void addElement(const label);
            
            //- remove element from subset
            inline void removeElement(const label);
            
            //- update subset after mesh modification
            //- the list contains new labels of mesh elements
            //- removed elements have a negative label
            template<class ListType>
            inline void updateSubset(const ListType&);
            
            //- update subset after modification
            //- this modifier is used in case when elements are decomposed
            inline void updateSubset(const VRWGraph&);
			
		// Searches
			//- find if the element exists in the subset
			inline bool contains(const label) const;

    // IOstream operators

			friend Ostream& operator<<(Ostream&, const meshSubset&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSubsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
