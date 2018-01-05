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

template<class T, Foam::label staticSize>
inline T* Foam::Module::DynList<T, staticSize>::data()
{
    return dataPtr_;
}


template<class T, Foam::label staticSize>
inline const T* Foam::Module::DynList<T, staticSize>::data() const
{
    return dataPtr_;
}


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::allocateSize(const label s)
{
    checkAllocation();

    if (s > staticSize)
    {
        if (s > nAllocated_)
        {
            //- allocates enough space for the elements
            T* newData = new T[s];

            for (label i = 0; i < nextFree_; ++i)
            {
                newData[i] = this->operator[](i);
            }

            if (nAllocated_ > staticSize)
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
        if (nAllocated_ > staticSize)
        {
            //- delete dynamically allocated data
            for (label i = 0; i < s; ++i)
            {
                staticData_[i] = dataPtr_[i];
            }

            delete[] dataPtr_;
        }

        dataPtr_ = staticData_;
        nAllocated_ = staticSize;
    }
}


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::checkIndex
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


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::checkAllocation() const
{
    if (nextFree_ > nAllocated_)
    {
        FatalErrorInFunction
            << "nextFree_ is out of scope 0 " << " and " << nAllocated_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::DynList()
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


template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::DynList(const label s)
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


template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::DynList
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


template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::DynList(const UList<T>& ul)
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


template<class T, Foam::label staticSize>
template<class ListType>
inline Foam::Module::DynList<T, staticSize>::DynList(const ListType& l)
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


template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::DynList
(
    const DynList<T, staticSize>& dl
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


template<class T, Foam::label staticSize>
inline Foam::Module::DynList<T, staticSize>::~DynList()
{
    allocateSize(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, Foam::label staticSize>
inline Foam::label Foam::Module::DynList<T, staticSize>::size() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return nextFree_;
}


template<class T, Foam::label staticSize>
inline Foam::label Foam::Module::DynList<T, staticSize>::byteSize() const
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


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::setSize(const label s)
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


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::clear()
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    nextFree_ = 0;
}


template<class T, Foam::label staticSize>
void Foam::Module::DynList<T, staticSize>::shrink()
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    allocateSize(nextFree_);

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::append(const T& e)
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


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::appendIfNotIn(const T& e)
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


template<class T, Foam::label staticSize>
inline Foam::label Foam::Module::DynList<T, staticSize>::find
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


template<class T, Foam::label staticSize>
inline bool Foam::Module::DynList<T, staticSize>::found(const T& e) const
{
    return this->find(e) != -1;
}



template<class T, Foam::label staticSize>
inline const T& Foam::Module::DynList<T, staticSize>::first() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](0);
}


template<class T, Foam::label staticSize>
inline const T& Foam::Module::DynList<T, staticSize>::last() const
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](nextFree_ - 1);
}


template<class T, Foam::label staticSize>
inline T Foam::Module::DynList<T, staticSize>::remove()
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


template<class T, Foam::label staticSize>
inline T Foam::Module::DynList<T, staticSize>::removeElement(const label i)
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


template<class T, Foam::label staticSize>
inline T& Foam::Module::DynList<T, staticSize>::newElmt(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator()(i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, Foam::label staticSize>
inline T& Foam::Module::DynList<T, staticSize>::operator()(const label i)
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


template<class T, Foam::label staticSize>
inline const T& Foam::Module::DynList<T, staticSize>::operator[]
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


template<class T, Foam::label staticSize>
inline T& Foam::Module::DynList<T, staticSize>::operator[](const label i)
{
    # ifdef FULLDEBUG
    checkAllocation();
    checkIndex(i);
    # endif

    return dataPtr_[i];
}


template<class T, Foam::label staticSize>
inline Foam::label Foam::Module::DynList<T, staticSize>::fcIndex
(
    const label index
) const
{
    return (index + 1) % nextFree_;
}


template<class T, Foam::label staticSize>
inline Foam::label Foam::Module::DynList<T, staticSize>::rcIndex
(
    const label index
) const
{
    return (index + nextFree_ - 1) % nextFree_;
}


template<class T, Foam::label staticSize>
inline const T& Foam::Module::DynList<T, staticSize>::fcValue
(
    const label index
) const
{
    return operator[](fcIndex(index));
}


template<class T, Foam::label staticSize>
inline const T& Foam::Module::DynList<T, staticSize>::rcValue
(
    const label index
) const
{
    return operator[](rcIndex(index));
}


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::operator=(const T& t)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    for (label i = 0; i < nextFree_; ++i)
    {
        operator[](i) = t;
    }
}


template<class T, Foam::label staticSize>
inline void Foam::Module::DynList<T, staticSize>::operator=
(
    const DynList<T, staticSize>& dl
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


template<class T, Foam::label staticSize>
inline bool Foam::Module::DynList<T, staticSize>::operator==
(
    const DynList<T, staticSize>& DL
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


template<class T, Foam::label staticSize>
inline bool Foam::Module::DynList<T, staticSize>::operator!=
(
    const DynList<T, staticSize>& DL
) const
{
    return !operator==(DL);
}


// ************************************************************************* //
