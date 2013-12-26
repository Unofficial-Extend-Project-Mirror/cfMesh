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

#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGen2DEngine::findActiveFacesAndPoints()
{
    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& points = mesh_.points();

    //- classify points
    const scalar tZ = 0.05 * (bb_.max().z() - bb_.min().z());

    activePoint_.setSize(points.size());
    label counter(0);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50) reduction(+ : counter)
    # endif
    forAll(points, pointI)
    {
        if( Foam::mag(points[pointI].z() - bb_.min().z()) < tZ )
        {
            activePoint_[pointI] = true;
            ++counter;
        }
        else
        {
            activePoint_[pointI] = false;
        }
    }

    if( 2 * counter != points.size() )
    {
        FatalErrorIn
        (
            "void polyMeshGen2DEngine::findActiveFacesAndPoints()"
        ) << "The number of points at smallest x coordinate is"
          << " not half of the total number of points."
          << " This is not a 2D mesh or is not aligned with the z axis"
          << exit(FatalError);
    }

    activePoints_.setSize(counter);
    offsetPoints_.setSize(counter);
    counter = 0;
    forAll(activePoint_, pointI)
        if( activePoint_[pointI] )
            activePoints_[counter++] = pointI;

    const VRWGraph& pointPoints = mesh_.addressingData().pointPoints();
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(activePoints_, pI)
    {
        const label pointI = activePoints_[pI];

        label nInactive(0), offsetPoint(-1);
        forAllRow(pointPoints, pointI, ppI)
        {
            if( !activePoint_[pointPoints(pointI, ppI)] )
            {
                ++nInactive;
                offsetPoint = pointPoints(pointI, ppI);
            }
        }

        if( nInactive == 1 )
        {
            offsetPoints_[pI] = offsetPoint;
        }
        else
        {
            FatalErrorIn
            (
                "void polyMeshGen2DEngine::findActiveFacesAndPoints()"
            ) << "This cannot be a 2D mesh" << exit(FatalError);
        }
    }

    //- find active faces
    activeFace_.setSize(faces.size());

    counter = 0;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50) reduction(+ : counter)
    # endif
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        bool allActive(true);
        forAll(f, pI)
        {
            if( !activePoint_[f[pI]] )
            {
                allActive = false;
                break;
            }
        }

        if( allActive )
        {
            activeFace_[faceI] = true;
            ++counter;
        }
        else
        {
            activeFace_[faceI] = false;
        }
    }

    //- set active faces
    activeFaces_.setSize(counter);
    forAll(activeFace_, faceI)
    {
        if( activeFace_[faceI] )
            activeFaces_[counter++] = faceI;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGen2DEngine::polyMeshGen2DEngine(polyMeshGen& mesh)
:
    mesh_(mesh),
    bb_(mesh.points()),
    activeFace_(),
    activeFaces_(),
    activePoint_(),
    activePoints_(),
    offsetPoints_()
{}

polyMeshGen2DEngine::~polyMeshGen2DEngine()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGen2DEngine::correctPoints()
{
    pointFieldPMG& points = mesh_.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(activePoints_, apI)
    {
        point& p = points[activePoints_[apI]];
        point& op = points[offsetPoints_[apI]];

        op.x() = p.x();
        op.y() = p.y();
        p.z() = bb_.min().z();
        op.z() = bb_.max().z();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
