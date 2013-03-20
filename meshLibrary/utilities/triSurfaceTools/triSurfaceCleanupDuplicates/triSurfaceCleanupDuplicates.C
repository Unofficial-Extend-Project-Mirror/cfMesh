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

#include "triSurfaceCleanupDuplicates.H"
#include "meshOctree.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceCleanupDuplicates::triSurfaceCleanupDuplicates
(
    const meshOctree& octree,
    const scalar tol
)
:
    tolerance_(tol),
    surf_(const_cast<triSurf&>(octree.surface())),
    octree_(octree),
    newTriangleLabel_(),
    done_(false)
{}

triSurfaceCleanupDuplicates::~triSurfaceCleanupDuplicates()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceCleanupDuplicates::mergeIdentities()
{
    if( Pstream::parRun() )
        FatalError << "Material detection does not run in parallel"
            << exit(FatalError);

    if( done_ )
    {
        WarningIn("void triSurfaceCleanupDuplicates::mergeIdentities()")
            << "Operation is already performed" << endl;
        return;
    }

    newTriangleLabel_.setSize(surf_.size());
    forAll(newTriangleLabel_, triI)
        newTriangleLabel_[triI] = triI;

    bool finished;
    do
    {
        finished = true;

        if( checkDuplicateTriangles() )
            finished = false;
        if( mergeDuplicatePoints() )
            finished = false;
    } while( !finished );

    done_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
