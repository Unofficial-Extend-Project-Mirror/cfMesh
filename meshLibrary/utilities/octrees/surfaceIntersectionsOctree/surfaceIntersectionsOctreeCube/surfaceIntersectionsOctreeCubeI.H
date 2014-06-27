/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description


\*---------------------------------------------------------------------------*/

#include "surfaceIntersectionsOctreeCube.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline bool surfaceIntersectionsOctreeCube::isLeaf() const
{
    if( subCubesPtr_ )
        return false;

    return true;
}

inline FixedList<surfaceIntersectionsOctreeCube*, 8>* surfaceIntersectionsOctreeCube::subCubes() const
{
    return subCubesPtr_;
}

inline const boundBox& surfaceIntersectionsOctreeCube::bb() const
{
    return cubeBox_;
}

inline void surfaceIntersectionsOctreeCube::append(const label tI)
{
    containedElements_.append(tI);
}

inline const SLList<label>& surfaceIntersectionsOctreeCube::containedElements() const
{
    return containedElements_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
