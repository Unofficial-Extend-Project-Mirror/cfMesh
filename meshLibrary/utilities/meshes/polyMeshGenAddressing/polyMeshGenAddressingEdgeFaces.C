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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcEdgeFaces() const
{
    if( efPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcEdges() const")
            << "edgeFaces already calculated"
            << abort(FatalError);
    }
    else
    {
        const faceListPMG& faces = mesh_.faces();
        const VRWGraph& pointFaces = this->pointFaces();
        const edgeList& edges = this->edges();

        efPtr_ = new VRWGraph();
        VRWGraph& edgeFaceAddr = *efPtr_;

        labelList nef(edges.size());

        # ifdef USE_OMP
        const label nThreads = 3 * omp_get_num_procs();
        # endif

        # ifdef USE_OMP
        # pragma omp parallel num_threads(nThreads) if( edges.size() > 10000 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(nef, edgeI)
                nef[edgeI] = 0;

            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(edges, edgeI)
            {
                const edge& ee = edges[edgeI];
                const label s = ee.start();

                forAllRow(pointFaces, s, pfI)
                {
                    const label faceI = pointFaces(s, pfI);

                    const face& f = faces[faceI];

                    forAll(f, eI)
                    {
                        if( f.faceEdge(eI) == ee )
                        {
                            ++nef[edgeI];
                            break;
                        }
                    }
                }
            }

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp master
            # endif
            VRWGraphSMPModifier(edgeFaceAddr).setSizeAndRowSize(nef);

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp for schedule(static)
            # endif
            forAll(edges, edgeI)
            {
                const edge& ee = edges[edgeI];
                const label s = ee.start();

                DynList<label> eFaces;
                forAllRow(pointFaces, s, pfI)
                {
                    const label faceI = pointFaces(s, pfI);

                    const face& f = faces[faceI];

                    forAll(f, eI)
                    {
                        if( f.faceEdge(eI) == ee )
                        {
                            eFaces.append(faceI);
                            break;
                        }
                    }
                }

                edgeFaceAddr.setRow(edgeI, eFaces);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::edgeFaces() const
{
    if( !efPtr_ )
        calcEdgeFaces();

    return *efPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
