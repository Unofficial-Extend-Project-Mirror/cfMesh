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
