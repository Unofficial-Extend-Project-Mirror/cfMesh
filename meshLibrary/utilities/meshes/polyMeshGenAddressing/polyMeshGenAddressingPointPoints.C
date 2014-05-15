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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcPointPoints() const
{
    if( ppPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcPointPoints() const")
            << "pointPoints already calculated"
            << abort(FatalError);
    }
    else
    {
        ppPtr_ = new VRWGraph();
        VRWGraph& pp = *ppPtr_;

        const faceListPMG& faces = mesh_.faces();
        const VRWGraph& pointFaces = this->pointFaces();

        labelList nPoints(pointFaces.size());

        # ifdef USE_OMP
        const label nThreads = 3 * omp_get_num_procs();
        # pragma omp parallel num_threads(nThreads) if( nPoints.size() > 10000 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(nPoints, i)
                nPoints[i] = 0;

            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(pointFaces, pointI)
            {
                DynList<label, 32> helper;

                forAllRow(pointFaces, pointI, pfI)
                {
                    const face& f = faces[pointFaces(pointI, pfI)];

                    const label pos = f.which(pointI);
                    helper.appendIfNotIn(f.prevLabel(pos));
                    helper.appendIfNotIn(f.nextLabel(pos));
                }

                nPoints[pointI] = helper.size();
            }

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp master
            # endif
            VRWGraphSMPModifier(pp).setSizeAndRowSize(nPoints);

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp for schedule(static)
            # endif
            forAll(pointFaces, pointI)
            {
                DynList<label, 32> helper;

                forAllRow(pointFaces, pointI, pfI)
                {
                    const face& f = faces[pointFaces(pointI, pfI)];

                    const label pos = f.which(pointI);
                    const label pLabel = f.prevLabel(pos);
                    const label nLabel = f.nextLabel(pos);

                    helper.appendIfNotIn(nLabel);
                    helper.appendIfNotIn(pLabel);
                }

                pp.setRow(pointI, helper);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::pointPoints() const
{
    if( !ppPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::pointPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcPointPoints();
    }

    return *ppPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
