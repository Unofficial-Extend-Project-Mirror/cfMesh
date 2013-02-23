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

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceEdges() const
{
    if( fePtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcFaceEdges() const")
            << "faceEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        fePtr_ = new VRWGraph();
        VRWGraph& faceEdgesAddr = *fePtr_;
        
        const edgeList& edges = this->edges();
        
        const VRWGraph& pointFaces = this->pointFaces();
        const faceListPMG& faces = mesh_.faces();
        
        labelList nfe(faces.size());
        
        const label nThreads = 3 * omp_get_num_procs();
        
        # pragma omp parallel num_threads(nThreads) if( faces.size() > 10000 )
        {
            # pragma omp for schedule(static)
            forAll(faces, faceI)
                nfe[faceI] = faces[faceI].size();
            
            # pragma omp barrier
            
            # pragma omp master
            VRWGraphSMPModifier(faceEdgesAddr).setSizeAndRowSize(nfe);
            
            # pragma omp barrier
            
            # pragma omp for schedule(static)
            forAll(edges, edgeI)
            {
                const edge ee = edges[edgeI];
                const label s = ee.start();
                
                forAllRow(pointFaces, s, pfI)
                {
                    const label faceI = pointFaces(s, pfI);
                    
                    const face& f = faces[faceI];
                    forAll(f, eI)
                    {
                        if( f.faceEdge(eI) == ee )
                        {
                            faceEdgesAddr[faceI][eI] = edgeI;
                            break;
                        }
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::faceEdges() const
{
    if( !fePtr_ )
        calcFaceEdges();

    return *fePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
