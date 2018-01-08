/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Module::DynList

Description
    A dynamic list is a 1-D vector of objects of type T which resizes
    itself as necessary to accept the new objects.  Internal storage
    is a compact array and the list can be shrunk to compact storage.
    The increase of list size is controlled by three template parameters,
    which allows the list storage to either increase by the given increment
    or the given multiplier and divider(allowing non-integer multiples).

SourceFiles
    DynListI.H
    DynList.C

\*---------------------------------------------------------------------------*/

#ifndef DynList_H
#define DynList_H

#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

// Forward declarations

template<class T, int StaticSize> class DynList;

template<class T, int StaticSize>
Ostream& operator<<
(
    Ostream&,
    const DynList<T, StaticSize>&
);
template<class T, int StaticSize>
Istream& operator>>
(
    Istream&,
    DynList<T, StaticSize>&
);


/*---------------------------------------------------------------------------*\
                           Class DynList Declaration
\*---------------------------------------------------------------------------*/

template<class T, int StaticSize = 16>
class DynList
{
    // Private data

        //- pointer to the data
        T* dataPtr_;

        //- size of the allocated data
        label nAllocated_;

        //- statically allocated data (used for short lists)
        T staticData_[StaticSize];

        //- Number of next free element
        label nextFree_;


    // Private member functions

        //- access to the data pointer
        inline T* data();

        //- const access to the data pointer
        inline const T* data() const;

        //- allocate list size
        inline void allocateSize(const label);

        //- check if index is inside the scope (used for debugging only)
        inline void checkIndex(const label) const;

        //- check if nAllocated_ is greater or equal to nextFree_
        inline void checkAllocation() const;


public:

    // Constructors

        //- Construct null
        inline DynList();

        //- Construct given size
        explicit inline DynList(const label);

        //- Construct given integer size
        #if WM_LABEL_SIZE == 64
        explicit inline DynList(const int32_t nElem)
        :
            DynList(label(nElem))
        {}
        #endif

        //- Construct from given size and defualt value
        explicit inline DynList(const label, const T&);

        //- Construct from UList. nextFree_ set to size().
        explicit inline DynList(const UList<T>&);

        //- Construct from other ListType
        template<class ListType>
        inline DynList(const ListType&);

        //- Copy constructor
        inline DynList(const DynList<T, StaticSize>&);

        //- Construct from Istream. nextFree_ set to size().
        explicit DynList(Istream&);

    //- Destructor
    inline ~DynList();


    // Member Functions

        // Access

            //- Size of the active part of the list.
            //- Direct over-ride of list size member function
            inline label size() const;

            //- Number of bytes used by the active part of the list
            //- Direct over-ride of list byteSize member function
            inline label byteSize() const;


        // Edit

            //- Reset size of List.
            void setSize(const label);

            //- Clear the list, i.e. set next free to zero.
            //  Allocated size does not change
            void clear();

            //- Shrink the List<T> to the number of elements used
            void shrink();


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Append an element at the end of the list if it is not yet
        //- present in the list (takes linear time)
        inline void appendIfNotIn(const T& e);

        //- Find and return location of element in list
        inline label find(const T& val) const;

        //- Check if the element is in the list
        inline bool found(const T& e) const;

        //- return a const reference to the first element
        inline const T& first() const;

        //- return a const reference to the last element
        inline const T& last() const;

        //- Return and remove the last element
        inline T remove();

        //- Return and remove the element
        inline T removeElement(const label i);

        //- return a reference to the element. Resize the list if necessary
        inline T& newElmt(const label);

        //- Return non-const access to an element,
        //-  resizing the list if necessary
        inline T& operator()(const label);

        //- return access to an element
        inline const T& operator[](const label) const;
        inline T& operator[](const label);

        //- return forward and reverse circular indices
        inline label fcIndex(const label index) const;
        inline label rcIndex(const label index) const;

        //- return forward and reverse circular elements
        inline const T& fcValue(const label index) const;
        inline const T& rcValue(const label index) const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);

        //- Copy of another list
        inline void operator=(const DynList<T, StaticSize>&);

        //- Compare the list with the another one
        inline bool operator==(const DynList<T, StaticSize>&) const;
        inline bool operator!=(const DynList<T, StaticSize>&) const;


    // IOstream operators

        // Write DynList to Ostream.
        friend Ostream& operator<< <T, StaticSize>
        (
            Ostream&,
            const DynList<T, StaticSize>&
        );

        //- Read from Istream, discarding contents of existing DynList.
        friend Istream& operator>> <T, StaticSize>
        (
            Istream&,
            DynList<T, StaticSize>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
