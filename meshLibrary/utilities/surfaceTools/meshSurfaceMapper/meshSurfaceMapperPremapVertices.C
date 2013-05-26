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

#include "demandDrivenData.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshOctree.H"
#include "refLabelledPoint.H"
#include "helperFunctionsPar.H"

#include <omp.h>

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::preMapVertices(const label nIterations)
{
    Info << "Smoothing mesh surface before mapping. Iteration:" << flush;

    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();

    List<labelledPoint> preMapPositions(boundaryPoints.size());

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        //- use the shrinking laplace first
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40)
        # endif
        forAll(pointFaces, bpI)
        {
            labelledPoint lp(0, vector::zero);

            forAllRow(pointFaces, bpI, bfI)
            {
                ++lp.pointLabel();
                lp.coordinates() += faceCentres[pointFaces(bpI, bfI)];
            }

            preMapPositions[bpI] = lp;
        }

        if( Pstream::parRun() )
        {
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();

            //- collect data to be sent to other processors
            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(surfaceEngine_.bpNeiProcs(), i)
                exchangeData.insert
                (
                    std::make_pair
                    (
                        surfaceEngine_.bpNeiProcs()[i],
                        LongList<refLabelledPoint>()
                    )
                );

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append
                    (
                        refLabelledPoint
                        (
                            globalPointLabel[bpI],
                            preMapPositions[bpI]
                        )
                    );
                }
            }

            //- exchange data with other processors
            LongList<refLabelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- combine collected data with the available data
            forAll(receivedData, i)
            {
                const refLabelledPoint& rlp = receivedData[i];
                const labelledPoint& lps = rlp.lPoint();

                const label bpI = globalToLocal[rlp.objectLabel()];

                labelledPoint& lp = preMapPositions[bpI];
                lp.pointLabel() += lps.pointLabel();
                lp.coordinates() += lps.coordinates();
            }
        }

        //- calculate coordinates of points for searching
        # ifdef USE_OMP
        # pragma omp parallel for shared(size)
        # endif
        forAll(preMapPositions, bpI)
        {
            labelledPoint& lp = preMapPositions[bpI];

            if( lp.pointLabel() == 0 )
            {
                Warning << "Surface point " << bpI
                    << " has no supporting faces" << endl;
                continue;
            }

            lp.coordinates() /= lp.pointLabel();
        }

        //- create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

        # ifdef USE_OMP
        # pragma omp parallel for if( size > 1000 ) shared(size) \
        schedule(dynamic, Foam::max(100, size / (3 * omp_get_max_threads())))
        # endif
        forAll(boundaryPoints, bpI)
        {
            const point& p = points[boundaryPoints[bpI]];

            label patch;
            point pMap = p;
            scalar dSq;

            meshOctree_.findNearestSurfacePoint
            (
                pMap,
                dSq,
                patch,
                preMapPositions[bpI].coordinates()
            );

            const point newP = 0.5 * (pMap + p);

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
        }

        surfaceModifier.updateGeometry();
        surfaceModifier.syncVerticesAtParallelBoundaries();

        Info << "." << flush;
    }

    Info << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
