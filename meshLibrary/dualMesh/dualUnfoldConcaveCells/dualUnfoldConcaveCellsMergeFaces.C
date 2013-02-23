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
#include "dualUnfoldConcaveCells.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"
#include "primitiveMesh.H"
#include "tetrahedron.H"

//#define DEBUGEdges

# ifdef DEBUGEdges
#include "cellSet.H"
#include "objectRegistry.H"
#include "Time.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool dualUnfoldConcaveCells::findConcaveEdges(const meshSurfaceEngine& mse)
{
	const pointFieldPMG& points = mesh_.points();
	const faceList::subList& bFaces = mse.boundaryFaces();
	const edgeList& edges = mse.edges();
	const VRWGraph& edgeFaces = mse.edgeFaces();
	const labelList& faceOwner = mse.faceOwners();
	const labelList& facePatch = mse.boundaryFacePatches();
	
	//- mark all boundary cells
	forAll(faceOwner, bfI)
		typeOfCell_[faceOwner[bfI]] |= BOUNDARYCELL;
	
	//- find concave cells
	bool foundConcave(false);
	forAll(edgeFaces, edgeI)
	{
		const label fOwn = edgeFaces(edgeI, 0);
		const label fNei = edgeFaces(edgeI, 1);
		
		if(
			(faceOwner[fOwn] == faceOwner[fNei]) &&
			(facePatch[fOwn] != facePatch[fNei])
		)
		{
			if( !help::isSharedEdgeConvex(points, bFaces[fOwn], bFaces[fNei]) )
			{
				foundConcave = true;
					
				const edge& e = edges[edgeI];
				typeOfVertex_[e.start()] |= CONCAVE;
				typeOfVertex_[e.end()] |= CONCAVE;
				
				typeOfCell_[faceOwner[fOwn]] |= CONCAVECELL;
			}
		}
	}
	
	# ifdef DEBUGEdges
	cellSet badCells
	(
		IOobject
		(
			"concaveCells",
			mesh_.returnRegistry().time().constant(),
			"polyMesh/sets",
			mesh_.returnRegistry()
		)
	);
	
	forAll(typeOfCell_, cI)
		if( typeOfCell_[cI] & CONCAVECELL )
			badCells.insert(cI);
		
	badCells.write();
	//::exit(1);
	# endif
	
	return foundConcave;
}

void dualUnfoldConcaveCells::markVertexTypes(const meshSurfaceEngine& mse)
{	
	const labelList& bPoints = mse.boundaryPoints();
	const VRWGraph& pointFaces = mse.pointFaces();
	const labelList& facePatch = mse.boundaryFacePatches();
	
	forAll(pointFaces, bpI)
	{
		DynList<label> patches(5);
		forAllRow(pointFaces, bpI, pfI)
			patches.appendIfNotIn(facePatch[pointFaces(bpI, pfI)]);
		
		if( patches.size() > 2 )
		{
			typeOfVertex_[bPoints[bpI]] |= CORNER;
		}
		else if( patches.size() == 2 )
		{
			typeOfVertex_[bPoints[bpI]] |= EDGE;
		}
	}
}

label dualUnfoldConcaveCells::mergeBoundaryFacesOfCell(const label cellI)
{
	if( typeOfCell_[cellI] & MERGEDCELL )
	{
		FatalErrorIn
		(
			"label dualUnfoldConcaveCells::"
			"mergeBoundaryFacesOfCell(const label cellI)"
		) << "Cell " << cellI << " has already been merged!"
			<< abort(FatalError);
	}
	else
	{
		typeOfCell_[cellI] |= MERGEDCELL;
	}
	
	const pointFieldPMG& points = mesh_.points();
	const faceListPMG& faces = mesh_.faces();
	const cell& c = mesh_.cells()[cellI];
	
	DynList<face> bFaces(5);
	DynList<label> facePatch(5);
	
	//- find boundary faces and corresponding patches
	forAll(c, fI)
	{
		const label fPatch = mesh_.faceIsInPatch(c[fI]);
		
		if( fPatch != -1 )
		{
			bFaces.append(faces[c[fI]]);
			facePatch.append(fPatch);
		}
	}
	
	# ifdef DEBUGEdges
	Info << "Boundary faces for cell " << cellI << " are " << bFaces << endl;
	Info << "Patches for boundary faces are " << facePatch << endl;
	# endif
	
	//- merge all boundary faces into one face
	//- this assumes that there was originally only one
	//- boundary face per cell. This is the cases after running
	//- surface morpher and topological cleaner
	boolList isFaceMerged(bFaces.size(), false);
	face mf(bFaces[0]);
	isFaceMerged[0] = true;
	bool merged;
	do
	{
		merged = false;
		for(label i=1;i<bFaces.size();++i)
			if(
				!isFaceMerged[i] &&
				help::shareAnEdge(mf, bFaces[i])
			)
			{
				merged = true;
				isFaceMerged[i] = true;
				mf = help::mergeTwoFaces(mf, bFaces[i]);
			}
	} while( merged );
	
	# ifdef DEBUGEdges
	Info << "Merged face for cell is " << mf << endl;
	# endif
	
	DynList<label> shrinkedFace(mf.size());
	forAll(mf, pI)
		if( !(typeOfVertex_[mf[pI]] & EDGE) )
		{
			shrinkedFace.append(mf[pI]);
		}
		else
		{
			# ifdef DEBUGEdges
			Info << "Marking vertex " << mf[pI] << " for removal!" << endl;
			# endif
			
			typeOfVertex_[mf[pI]] |= REMOVE;
		}
		
	shrinkedFace.shrink();
		
	//- face will be stored into a patch corresponding
	//- to the largest face
	label fPatch(facePatch[0]);
	scalar size(0.0);
	forAll(bFaces, bfI)
	{
		# ifdef DEBUGEdges
		Info << "Sizes of bundary face " << bfI << " is "
			<< bFaces[bfI].mag(points) << endl;
		# endif
		
		if( bFaces[bfI].mag(points) > size )
		{
			fPatch = facePatch[bfI];
			size = bFaces[bfI].mag(points);
		}
	}
		
	//- store shrinked face
	# ifdef DEBUGEdges
	Info << "Storing face " << shrinkedFace << " in patch"
		<< fPatch << endl;
	#endif
	
	newBoundaryFaces_.appendList(face(shrinkedFace));
	newBoundaryOwners_.append(cellI);
	newBoundaryPatches_.append(fPatch);
	
	return fPatch;
}

void dualUnfoldConcaveCells::storeAndMergeBoundaryFaces
(
	const meshSurfaceEngine& mse
)
{
	bool addSelected;
	do
	{
		addSelected = false;
		
		//- merge boundary faces for concave cells
		forAll(typeOfCell_, cellI)
			if( (typeOfCell_[cellI] & CONCAVECELL) &&
				!(typeOfCell_[cellI] & TREATEDCELL)
			)
			{
				# ifdef DEBUGEdges
				Info << "Merging faces of concave cell " << cellI << endl;
				# endif
				
				mergeBoundaryFacesOfCell(cellI);
	
				//- set flag that the cell has been treated
				typeOfCell_[cellI] |= TREATEDCELL;
			}
		
		const VRWGraph& edgeFaces = mse.edgeFaces();
		const edgeList& edges = mse.edges();
		const labelList& faceOwner = mse.faceOwners();
		const labelList& facePatch = mse.boundaryFacePatches();		
		forAll(edgeFaces, edgeI)
		{
			const label fOwn = edgeFaces(edgeI, 0);
			const label fNei = edgeFaces(edgeI, 1);
			
			if(
				(faceOwner[fOwn] == faceOwner[fNei]) &&
				(facePatch[fOwn] != facePatch[fNei])
			)
			{
				const edge& e = edges[edgeI];
				if(
					(typeOfVertex_[e[0]] & REMOVE) &&
					(typeOfVertex_[e[1]] & REMOVE) &&
					!(typeOfCell_[faceOwner[fOwn]] & CONCAVECELL)
				)
				{
					typeOfCell_[faceOwner[fOwn]] |= CONCAVECELL;
					addSelected = true;
					
					# ifdef DEBUGEdges
					Info << "Cell " << faceOwner[fOwn] << " is also concave"
						<< endl;
					# endif
				}
			}
		}
	}
	while( addSelected );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
