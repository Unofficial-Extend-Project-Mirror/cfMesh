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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

#include <omp.h>

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
meshSurfaceEngine::meshSurfaceEngine(polyMeshGen& mesh)
:
    mesh_(mesh),
    boundaryPointsPtr_(NULL),
    boundaryFacesPtr_(NULL),
    boundaryFacePatchPtr_(NULL),
    boundaryFaceOwnersPtr_(NULL),
    pointFacesPtr_(NULL),
    pointInFacePtr_(NULL),
    pointPatchesPtr_(NULL),
    bppPtr_(NULL),
    pointPointsPtr_(NULL),
    edgesPtr_(NULL),
    bpEdgesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    faceFacesPtr_(NULL),
    pointNormalsPtr_(NULL),
    faceNormalsPtr_(NULL),
    faceCentresPtr_(NULL),
    
    globalBoundaryPointLabelPtr_(NULL),
    globalBoundaryPointToLocalPtr_(NULL),
    bpProcsPtr_(NULL),
    bpNeiProcsPtr_(NULL),
    globalBoundaryEdgeLabelPtr_(NULL),
    globalBoundaryEdgeToLocalPtr_(NULL),
    beProcsPtr_(NULL),
    beNeiProcsPtr_(NULL),
    otherEdgeFaceAtProcPtr_(NULL),
    otherEdgeFacePatchPtr_(NULL),
    globalBoundaryFaceLabelPtr_(NULL)
{
    calculateBoundaryFaces();
    calculateBoundaryNodes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEngine::~meshSurfaceEngine()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const polyMeshGen& meshSurfaceEngine::mesh() const
{
    return mesh_;
}

const pointFieldPMG& meshSurfaceEngine::points() const
{
    return mesh_.points();
}

const faceListPMG& meshSurfaceEngine::faces() const
{
    return mesh_.faces();
}

const cellListPMG& meshSurfaceEngine::cells() const
{
    return mesh_.cells();
}

const labelList& meshSurfaceEngine::bp() const
{
    if( !bppPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::bp() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateBoundaryFaces();
        calculateBoundaryNodes();
    }

    return *bppPtr_;
}

const labelList& meshSurfaceEngine::boundaryPoints() const
{
    if( !boundaryPointsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::boundaryPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateBoundaryNodes();
    }

    return *boundaryPointsPtr_;
}

const faceList::subList& meshSurfaceEngine::boundaryFaces() const
{
    if( !boundaryFacesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const faceList::subList&"
                "meshSurfaceEngine::boundaryFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateBoundaryFaces();
    }

    return *boundaryFacesPtr_;
}

const labelList& meshSurfaceEngine::boundaryFacePatches() const
{
    if( !boundaryFacePatchPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::boundaryFacePatches() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateBoundaryFacePatches();
    }
    
    return *boundaryFacePatchPtr_;
}

const labelList& meshSurfaceEngine::faceOwners() const
{
    if( !boundaryFaceOwnersPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::faceOwners() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateBoundaryOwners();
    }

    return *boundaryFaceOwnersPtr_;
}

const VRWGraph& meshSurfaceEngine::pointFaces() const
{
    if( !pointFacesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculatePointFaces();
    }

    return *pointFacesPtr_;
}

const VRWGraph& meshSurfaceEngine::pointInFaces() const
{
    if( !pointInFacePtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointInFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculatePointFaces();
    }

    return *pointInFacePtr_;
}

const VRWGraph& meshSurfaceEngine::pointPatches() const
{
    if( !pointPatchesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointPatches() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculatePointPatches();
    }
    
    return *pointPatchesPtr_;
}

const VRWGraph& meshSurfaceEngine::pointPoints() const
{
    if( !pointPointsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculatePointPoints();
    }

    return *pointPointsPtr_;
}

const vectorField& meshSurfaceEngine::pointNormals() const
{
    if( !pointNormalsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::pointNormals() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}

const vectorField& meshSurfaceEngine::faceNormals() const
{
    if( !faceNormalsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::faceNormals() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateFaceNormals();
    }
    
    return *faceNormalsPtr_;
}

const vectorField& meshSurfaceEngine::faceCentres() const
{
    if( !faceCentresPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::faceCentres() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateFaceCentres();
    }
    
    return *faceCentresPtr_;
}

const edgeList& meshSurfaceEngine::edges() const
{
    if( !edgesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const edgeList& meshSurfaceEngine::edges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateEdgesAndAddressing();
    }

    return *edgesPtr_;
}

const VRWGraph& meshSurfaceEngine::boundaryPointEdges() const
{
    if( !bpEdgesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::boundaryPointEdges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateEdgesAndAddressing();
    }
    
    return *bpEdgesPtr_;
}

const VRWGraph& meshSurfaceEngine::edgeFaces() const
{
    if( !edgeFacesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::edgeFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateEdgeFacesAddressing();
    }

    return *edgeFacesPtr_;
}

const VRWGraph& meshSurfaceEngine::faceEdges() const
{
    if( !faceEdgesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::faceEdges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateFaceEdgesAddressing();
    }

    return *faceEdgesPtr_;
}

const VRWGraph& meshSurfaceEngine::faceFaces() const
{
    if( !faceFacesPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::faceFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calculateFaceFacesAddressing();
    }

    return *faceFacesPtr_;
}

const labelList& meshSurfaceEngine::globalBoundaryPointLabel() const
{
    if( !globalBoundaryPointLabelPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryPointLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryPointLabels();
    }
    
    return *globalBoundaryPointLabelPtr_;
}

const Map<label>& meshSurfaceEngine::globalToLocalBndPointAddressing() const
{
    if( !globalBoundaryPointToLocalPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::globalToLocalBndPointAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryPointLabels();
    }
    
    return *globalBoundaryPointToLocalPtr_;
}

const VRWGraph& meshSurfaceEngine::bpAtProcs() const
{
    if( !globalBoundaryPointLabelPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::bpAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryPointLabels();
    }
    
    return *bpProcsPtr_;
}

const DynList<label>& meshSurfaceEngine::bpNeiProcs() const
{
    if( !bpNeiProcsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& meshSurfaceEngine::bpNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryPointLabels();
    }
    
    return *bpNeiProcsPtr_;
}

const labelList& meshSurfaceEngine::globalBoundaryEdgeLabel() const
{
    if( !globalBoundaryEdgeLabelPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryEdgeLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryEdgeLabels();
    }
    
    return *globalBoundaryEdgeLabelPtr_;
}

const Map<label>& meshSurfaceEngine::globalToLocalBndEdgeAddressing() const
{
    if( !globalBoundaryEdgeToLocalPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::globalToLocalBndEdgeAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryEdgeLabels();
    }
        
    return *globalBoundaryEdgeToLocalPtr_;
}

const VRWGraph& meshSurfaceEngine::beAtProcs() const
{
    if( !globalBoundaryEdgeLabelPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::beAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryEdgeLabels();
    }
    
    return *beProcsPtr_;
}

const DynList<label>& meshSurfaceEngine::beNeiProcs() const
{
    if( !beNeiProcsPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& meshSurfaceEngine::beNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryEdgeLabels();
    }
    
    return *beNeiProcsPtr_;
}

const Map<label>& meshSurfaceEngine::otherEdgeFaceAtProc() const
{
    if( !otherEdgeFaceAtProcPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::otherEdgeFaceAtProc() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcAddressingForProcEdges();
    }
    
    return *otherEdgeFaceAtProcPtr_;
}

const Map<label>& meshSurfaceEngine::otherEdgeFacePatch() const
{
    if( !otherEdgeFacePatchPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::otherEdgeFacePatch() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcAddressingForProcEdges();
    }
    
    return *otherEdgeFacePatchPtr_;
}

const labelList& meshSurfaceEngine::globalBoundaryFaceLabel() const
{
    if( !globalBoundaryFaceLabelPtr_ )
    {
        # ifdef FULLDEBUG
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryFaceLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif
        
        calcGlobalBoundaryFaceLabels();
    }
    
    return *globalBoundaryFaceLabelPtr_;
}

void meshSurfaceEngine::clearOut()
{
    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(boundaryFacesPtr_);
    deleteDemandDrivenData(boundaryFacePatchPtr_);
    deleteDemandDrivenData(boundaryFaceOwnersPtr_);
    deleteDemandDrivenData(pointFacesPtr_);
    deleteDemandDrivenData(pointInFacePtr_);
    deleteDemandDrivenData(pointPatchesPtr_);
    deleteDemandDrivenData(bppPtr_);
    deleteDemandDrivenData(pointPointsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(bpEdgesPtr_);
    deleteDemandDrivenData(edgeFacesPtr_);
    deleteDemandDrivenData(faceEdgesPtr_);
    deleteDemandDrivenData(faceFacesPtr_);
    
    deleteDemandDrivenData(globalBoundaryPointLabelPtr_);
    deleteDemandDrivenData(globalBoundaryPointToLocalPtr_);
    deleteDemandDrivenData(bpProcsPtr_);
    deleteDemandDrivenData(bpNeiProcsPtr_);
    deleteDemandDrivenData(globalBoundaryEdgeLabelPtr_);
    deleteDemandDrivenData(globalBoundaryEdgeToLocalPtr_);
    deleteDemandDrivenData(beProcsPtr_);
    deleteDemandDrivenData(beNeiProcsPtr_);
    deleteDemandDrivenData(otherEdgeFaceAtProcPtr_);
    deleteDemandDrivenData(otherEdgeFacePatchPtr_);
    deleteDemandDrivenData(globalBoundaryFaceLabelPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
