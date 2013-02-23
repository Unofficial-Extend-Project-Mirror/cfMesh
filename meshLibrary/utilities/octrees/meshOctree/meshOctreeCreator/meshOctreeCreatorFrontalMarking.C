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

#include "meshOctreeCreator.H"
#include "meshOctreeInsideOutside.H"

//#define DEBUGSearch
# ifdef DEBUGSearch
#include "writeOctreeEnsight.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCreator::createInsideOutsideInformation()
{
    Info << "Marking inside/outside." << endl;

    meshOctreeInsideOutside insideOutside(octree_);

    # ifdef DEBUGSearch
    writeOctreeEnsight(octree_, "octreeInternal", meshOctreeCubeBasic::INSIDE);
    writeOctreeEnsight(octree_, "octreeOutside", meshOctreeCubeBasic::OUTSIDE);
    writeOctreeEnsight(octree_, "octreeData", meshOctreeCubeBasic::DATA);
    writeOctreeEnsight(octree_, "octreeUnknown", meshOctreeCubeBasic::UNKNOWN);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
