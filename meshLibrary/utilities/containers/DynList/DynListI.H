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

\*---------------------------------------------------------------------------*/

template<class T, int StaticSize>
inline T* Foam::Module::DynList<T, StaticSize>::data()
{
    return dataPtr_;
}


template<class T, int StaticSize>
inline const T* Foam::Module::DynList<T, StaticSize>::data() const
{
    return dataPtr_;
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::allocateSize(const label s)
{
    checkAllocation();

    if (s > StaticSize)
    {
        if (s > nAllocated_)
        {
            //- allocates enough space for the elements
            T* newData = new T[s];

            for (label i = 0; i < nextFree_; ++i)
            {
                newData[i] = this->operator[](i);
            }

            if (nAllocated_ > StaticSize)
            {
                delete[] dataPtr_;
            }

            dataPtr_ = newData;
            nAllocated_ = s;
        }
        else if (s < nAllocated_)
        {
            //- shrinks the list
            T* newData = new T[s];

            for (label i = 0; i < s; ++i)
            {
                newData[i] = this->operator[](i);
            }

            delete[] dataPtr_;

            dataPtr_ = newData;
            nAllocated_ = s;
        }
    }
    else
    {
        if (nAllocated_ > StaticSize)
        {
            //- delete dynamically allocated data
            for (label i = 0; i < s; ++i)
            {
                staticData_[i] = dataPtr_[i];
            }

            delete[] dataPtr_;
        }

        dataPtr_ = staticData_;
        nAllocated_ = StaticSize;
    }
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::checkIndex
(
    const label i
) const
{
    if ((i < 0) || (i >= nextFree_))
    {
        FatalErrorInFunction
            << "Index " << i << " is not in range " << 0
            << " and " << nextFree_ << abort(FatalError);
    }
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::checkAllocation() const
{
    if (nextFree_ > nAllocated_)
    {
        FatalErrorInFunction
            << "nextFree_ is out of scope 0 " << " and " << nAllocated_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::DynList()
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(0);

    # ifdef DEBUG
    checkAllocation();
    # endif

}


template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::DynList(const label s)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(s);

    # ifdef DEBUG
    checkAllocation();
    # endif

}


template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::DynList
(
    const label s,
    const T& val
)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(s);

    for (label i = 0; i < s; ++i)
    {
        this->operator[](i) = val;
    }

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::DynList(const UList<T>& ul)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(ul.size());

    forAll(ul, i)
    {
        this->operator[](i) = ul[i];
    }

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
template<class ListType>
inline Foam::Module::DynList<T, StaticSize>::DynList(const ListType& l)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(l.size());
    for (label i = 0; i < nextFree_; ++i)
    {
        this->operator[](i) = l[i];
    }

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::DynList
(
    const DynList<T, StaticSize>& dl
)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    setSize(dl.size());
    for (label i = 0; i < nextFree_; ++i)
    {
        this->operator[](i) = dl[i];
    }

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline Foam::Module::DynList<T, StaticSize>::~DynList()
{
    allocateSize(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int StaticSize>
inline Foam::label Foam::Module::DynList<T, StaticSize>::size() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return nextFree_;
}


template<class T, int StaticSize>
inline Foam::label Foam::Module::DynList<T, StaticSize>::byteSize() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (!contiguous<T>())
    {
        FatalErrorInFunction
            << "Cannot return the binary size of a list of "
            << "non - primitive elements"
            << abort(FatalError);
    }

    return nextFree_*sizeof(T);
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::setSize(const label s)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    allocateSize(s);
    nextFree_ = s;

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::clear()
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    nextFree_ = 0;
}


template<class T, int StaticSize>
void Foam::Module::DynList<T, StaticSize>::shrink()
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    allocateSize(nextFree_);

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::append(const T& e)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (nextFree_ >= nAllocated_)
    {
        const label newSize = 2*nAllocated_ + 2;
        allocateSize(newSize);
    }

    # ifdef DEBUG
    checkAllocation();
    # endif

    this->operator[](nextFree_++) = e;
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::appendIfNotIn(const T& e)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (!found(e))
    {
         append(e);
    }

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int StaticSize>
inline Foam::label Foam::Module::DynList<T, StaticSize>::find
(
    const T& e
) const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    for (label i = 0; i < nextFree_; ++i)
    {
        if (this->operator[](i) == e)
        {
            return i;
        }
    }

    return -1;
}


template<class T, int StaticSize>
inline bool Foam::Module::DynList<T, StaticSize>::found(const T& e) const
{
    return this->find(e) != -1;
}



template<class T, int StaticSize>
inline const T& Foam::Module::DynList<T, StaticSize>::first() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](0);
}


template<class T, int StaticSize>
inline const T& Foam::Module::DynList<T, StaticSize>::last() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](nextFree_ - 1);
}


template<class T, int StaticSize>
inline T Foam::Module::DynList<T, StaticSize>::remove()
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (nextFree_ == 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    T el = this->operator[](nextFree_ - 1);
    --nextFree_;
    return el;
}


template<class T, int StaticSize>
inline T Foam::Module::DynList<T, StaticSize>::removeElement(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (nextFree_ == 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    T el = this->operator[](i);
    this->operator[](i) = this->operator[](nextFree_ - 1);
    --nextFree_;

    # ifdef DEBUG
    checkAllocation();
    # endif

    return el;
}


template<class T, int StaticSize>
inline T& Foam::Module::DynList<T, StaticSize>::newElmt(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator()(i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int StaticSize>
inline T& Foam::Module::DynList<T, StaticSize>::operator()(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    nextFree_ = Foam::max(nextFree_, i + 1);

    if (nextFree_ >= nAllocated_)
    {
        allocateSize(2*nextFree_ + 1);
    }

    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](i);
}


template<class T, int StaticSize>
inline const T& Foam::Module::DynList<T, StaticSize>::operator[]
(
    const label i
) const
{
    # ifdef FULLDEBUG
    checkAllocation();
    checkIndex(i);
    # endif

    return dataPtr_[i];
}


template<class T, int StaticSize>
inline T& Foam::Module::DynList<T, StaticSize>::operator[](const label i)
{
    # ifdef FULLDEBUG
    checkAllocation();
    checkIndex(i);
    # endif

    return dataPtr_[i];
}


template<class T, int StaticSize>
inline Foam::label Foam::Module::DynList<T, StaticSize>::fcIndex
(
    const label index
) const
{
    return (index + 1) % nextFree_;
}


template<class T, int StaticSize>
inline Foam::label Foam::Module::DynList<T, StaticSize>::rcIndex
(
    const label index
) const
{
    return (index + nextFree_ - 1) % nextFree_;
}


template<class T, int StaticSize>
inline const T& Foam::Module::DynList<T, StaticSize>::fcValue
(
    const label index
) const
{
    return operator[](fcIndex(index));
}


template<class T, int StaticSize>
inline const T& Foam::Module::DynList<T, StaticSize>::rcValue
(
    const label index
) const
{
    return operator[](rcIndex(index));
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::operator=(const T& t)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    for (label i = 0; i < nextFree_; ++i)
    {
        operator[](i) = t;
    }
}


template<class T, int StaticSize>
inline void Foam::Module::DynList<T, StaticSize>::operator=
(
    const DynList<T, StaticSize>& dl
)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    allocateSize(dl.size());
    nextFree_ = dl.size();

    # ifdef DEBUG
    checkAllocation();
    # endif

    for (label i = 0; i < nextFree_; ++i)
    {
        this->operator[](i) = dl[i];
    }
}


template<class T, int StaticSize>
inline bool Foam::Module::DynList<T, StaticSize>::operator==
(
    const DynList<T, StaticSize>& DL
) const
{
    if (nextFree_ != DL.nextFree_)
    {
        return false;
    }

    forAll(DL, i)
    {
        if (this->operator[](i) != DL[i])
        {
            return false;
        }
    }

    return true;
}


template<class T, int StaticSize>
inline bool Foam::Module::DynList<T, StaticSize>::operator!=
(
    const DynList<T, StaticSize>& DL
) const
{
    return !operator==(DL);
}


// ************************************************************************* //
