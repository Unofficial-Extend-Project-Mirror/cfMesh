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

Class
    meshSurfaceEngine

Description
    Calculates surface of the mesh

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngine_H
#define meshSurfaceEngine_H

#include "polyMeshGenModifier.H"
#include "SubList.H"
#include "boolList.H"
#include "DynList.H"
#include "VRWGraph.H"
#include "triangle.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meshSurfaceEngine Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngine
{
    // Private data
        //- mesh
		polyMeshGen& mesh_;

        //- boundary points
        mutable labelList* boundaryPointsPtr_;

        //- boundary faces
        mutable faceList::subList* boundaryFacesPtr_;
	
		//- patches boundary faces are in
		mutable labelList* boundaryFacePatchPtr_;

        //- face owners
        mutable labelList* boundaryFaceOwnersPtr_;

        //- point faces addressing
        mutable VRWGraph* pointFacesPtr_;
        mutable VRWGraph* pointInFacePtr_;
    
        //- point-patches addressing
        mutable VRWGraph* pointPatchesPtr_;

        //- pointBoundaryPoint addressing
        mutable labelList* bppPtr_;

        //- point points addressing
        mutable VRWGraph* pointPointsPtr_;

        //- edges 
        mutable edgeList* edgesPtr_;
        
        //- boundary point-edges addressing
        mutable VRWGraph* bpEdgesPtr_;

        //- edge faces addressing
        mutable VRWGraph* edgeFacesPtr_;
		
		//- face edges addressing
		mutable VRWGraph* faceEdgesPtr_;
		
		//- face-faces addressing
		mutable VRWGraph* faceFacesPtr_;

        //- point normals
        mutable vectorField* pointNormalsPtr_;
		
		//- face normals
		mutable vectorField* faceNormalsPtr_;
        
        //- face centres
        mutable vectorField* faceCentresPtr_;
        
    // Private data for parallel execution
    
        //- global boundary point label
        mutable labelList* globalBoundaryPointLabelPtr_;
        
        //- global boundary point to local addressing
        mutable Map<label>* globalBoundaryPointToLocalPtr_;
        
        //- boundary point-processors addressing
        mutable VRWGraph* bpProcsPtr_;
        
        //- neighbour processors for communication when sending point data
        mutable DynList<label>* bpNeiProcsPtr_;
        
        //- global boundary edge label
        mutable labelList* globalBoundaryEdgeLabelPtr_;
        
        //- global boundary edge to local addressing
        mutable Map<label>* globalBoundaryEdgeToLocalPtr_;
        
        //- boundary edge-processors addressing
        mutable VRWGraph* beProcsPtr_;
        
        //- neighbour processors for communication when sending edge data
        mutable DynList<label>* beNeiProcsPtr_;
        
        //- processor containing other face and face-patch addressing
        mutable Map<label>* otherEdgeFaceAtProcPtr_;
        mutable Map<label>* otherEdgeFacePatchPtr_;
        
        //- global label for boundary faces
        mutable labelList* globalBoundaryFaceLabelPtr_;

    // Private member functions

        //- calculate boundary nodes, faces and addressing
        void calculateBoundaryFaces() const;
		void calculateBoundaryOwners() const;
		void calculateBoundaryNodes() const;
		
		void calculateBoundaryFacePatches() const;

        void calculatePointFaces() const;
        
        void calculatePointPatches() const;

        void calculatePointPoints() const;

        void calculatePointNormals() const;
		void calculateFaceNormals() const;
        void calculateFaceCentres() const;
        
        void updatePointNormalsAtProcBoundaries() const;

        void calculateEdgesAndAddressing() const;
        void calculateFaceEdgesAddressing() const;
		void calculateEdgeFacesAddressing() const;
		void calculateFaceFacesAddressing() const;
        
    // Private member functions for parallel runs
    
        void calcGlobalBoundaryPointLabels() const;
        void calcGlobalBoundaryEdgeLabels() const;
        void calcAddressingForProcEdges() const;
        void calcGlobalBoundaryFaceLabels() const;
        
    // Private member functions

        void clearOut();

        //- Disallow default bitwise copy construct
        meshSurfaceEngine(const meshSurfaceEngine&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngine&);

public:

		friend class meshSurfaceEngineModifier;
	
    // Constructors

        //- Construct from mesh
        meshSurfaceEngine(polyMeshGen& mesh);

    // Destructor

        ~meshSurfaceEngine();

    // Member Functions
        const polyMeshGen& mesh() const;
        const pointFieldPMG& points() const;
        const faceListPMG& faces() const;
        const cellListPMG& cells() const;

        const labelList& boundaryPoints() const;

        const faceList::subList& boundaryFaces() const;
		
		//- patch label for each boundary face
		const labelList& boundaryFacePatches() const;

        const labelList& faceOwners() const;

        const VRWGraph& pointFaces() const;
        const VRWGraph& pointInFaces() const;
        
        const VRWGraph& pointPatches() const;

        const VRWGraph& pointPoints() const;

        const vectorField& pointNormals() const;
		
		const vectorField& faceNormals() const;
        
        const vectorField& faceCentres() const;

        const labelList& bp() const;

        const edgeList& edges() const;
        
        const VRWGraph& boundaryPointEdges() const;

        const VRWGraph& edgeFaces() const;
		
		const VRWGraph& faceEdges() const;
		
		const VRWGraph& faceFaces() const;
        
    // Functions for parallel runs
        //- global boundary point label
        const labelList& globalBoundaryPointLabel() const;
        
        //- global point label to local label. Only for processors points
        const Map<label>& globalToLocalBndPointAddressing() const;
        
        //- processors which contain the vertex
        const VRWGraph& bpAtProcs() const;
        
        //- communication matrix for sending point data
        const DynList<label>& bpNeiProcs() const;
        
        //- global boundary edge label
        const labelList& globalBoundaryEdgeLabel() const;
        
        //- global boundary edge label to local label. Only for processor edges
        const Map<label>& globalToLocalBndEdgeAddressing() const;
        
        //- processors which contain the edges
        const VRWGraph& beAtProcs() const;
        
        //- communication matrix for sending edge data
        const DynList<label>& beNeiProcs() const;
        
        //- patch label and processor label at which the other surface face
        //- sharing processor edge is located. Only for processor edges
        const Map<label>& otherEdgeFaceAtProc() const;
        const Map<label>& otherEdgeFacePatch() const;
        
        //- global boundary face label
        const labelList& globalBoundaryFaceLabel() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
