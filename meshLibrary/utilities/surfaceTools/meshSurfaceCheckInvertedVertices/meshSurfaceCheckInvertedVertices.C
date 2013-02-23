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

#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfaceEngine.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "helperFunctionsPar.H"

#include <map>
#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
void meshSurfaceCheckInvertedVertices::checkVertices()
{
    const meshSurfaceEngine& mse = *surfaceEnginePtr_;
    const pointFieldPMG& points = mse.points();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& pointInFaces = mse.pointInFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const vectorField& fCentres = mse.faceCentres();
    
    invertedVertices_.clear();
    
    # pragma omp parallel for if( pointFaces.size() > 100 ) \
    schedule(dynamic, 20)
    forAll(pointFaces, bpI)
    {
        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI) )
            continue;
        
        forAllRow(pointFaces, bpI, pfI)
        {
            const label pI = pointInFaces(bpI, pfI);
            const label bfI = pointFaces(bpI, pfI);
            
            const face& bf = bFaces[bfI];
            
            //- chech the first triangle (with the next node)
            triangle<point, point> triNext
            (
                points[bf[pI]],
                points[bf.nextLabel(pI)],
                fCentres[bfI]
            );
            
            vector n = triNext.normal();
            scalar m = mag(n);
            if( m < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            else
            {
                n /= m;
            }
            
            if( magSqr(triNext.a() - triNext.b()) < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            if( magSqr(triNext.c() - triNext.a()) < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            
            if( (n & pNormals[bp[bf[pI]]]) < 0.0 )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            
            //- check the second triangle (with previous node)
            triangle<point, point> triPrev
            (
                points[bf[pI]],
                fCentres[bfI],
                points[bf.prevLabel(pI)]
            );
            
            n = triPrev.normal();
            m = mag(n);
            if( m < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            else
            {
                n /= m;
            }
            
            if( magSqr(triPrev.a() - triPrev.b()) < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            if( magSqr(triPrev.c() - triPrev.a()) < VSMALL )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            
            if( (n & pNormals[bp[bf[pI]]]) < 0.0 )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
                
                continue;
            }
            
            //- check whether the normals of both triangles
            //- point in the same direction
            if( (triNext.normal() & triPrev.normal()) < 0.0 )
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);
            }
        }
    }
    
    if( Pstream::parRun() )
    {
        //- exchange global labels of inverted points
        const labelList& bPoints = mse.boundaryPoints();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        
        std::map<label, labelListPMG> shareData;
        forAll(neiProcs, i)
            shareData.insert(std::make_pair(neiProcs[i], labelListPMG()));
        
        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();
            
            if( !invertedVertices_.found(bPoints[bpI]) )
                continue;
            
            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);
                
                if( neiProc == Pstream::myProcNo() )
                    continue;
                
                shareData[neiProc].append(globalPointLabel[bpI]);
            }
        }
        
        //- exchange data with other processors
        labelListPMG receivedData;
        help::exchangeMap(shareData, receivedData);
            
        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i]];
            invertedVertices_.insert(bPoints[bpI]);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfaceEngine& mse,
    const boolList* activePointsPtr
)
:
    deleteSurface_(false),
    surfaceEnginePtr_(&mse),
    activePointsPtr_(activePointsPtr),
    invertedVertices_()
{
	checkVertices();
}

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    polyMeshGen& mesh
)
:
    deleteSurface_(true),
    surfaceEnginePtr_(new meshSurfaceEngine(mesh)),
    activePointsPtr_(NULL),
    invertedVertices_()
{
	checkVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::~meshSurfaceCheckInvertedVertices()
{
    if( deleteSurface_ )
        deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
