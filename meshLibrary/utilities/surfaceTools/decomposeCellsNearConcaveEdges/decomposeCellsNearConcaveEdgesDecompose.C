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

#include "decomposeCellsNearConcaveEdges.H"
#include "demandDrivenData.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"

//#define DEBUGDec

# ifdef DEBUGDec
#include "Time.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void decomposeCellsNearConcaveEdges::findConcaveFacesAndCells()
{
	meshSurfaceEngine mse(mesh_);
	const faceList::subList& bFaces = mse.boundaryFaces();
	const labelList& facePatch = mse.boundaryFacePatches();
	const VRWGraph& edgeFaces = mse.edgeFaces();
	const edgeList& edges = mse.edges();
	const pointFieldPMG& points = mesh_.points();
	
	boolList concaveVertex(points.size(), false);
	forAll(edgeFaces, eI)
	{
		const label bfI = edgeFaces[eI][0];
		const label bfJ = edgeFaces[eI][1];
		
		const face& f0 = bFaces[bfI];
		const face& f1 = bFaces[bfJ];
		
		if(
			(facePatch[bfI] != facePatch[bfJ]) &&
			!help::isSharedEdgeConvex(points, f0, f1)
		)
		{
			const edge& e = edges[eI];
			concaveVertex[e[0]] = true;
			concaveVertex[e[1]] = true;
		}
	}
	
	//- decompose concave internal faces
	decomposeFaces(mesh_).decomposeConcaveInternalFaces(concaveVertex);
	
	//- find which cells need to be decomposed
	# ifdef DEBUGDec
	const label id = mesh_.addCellSubset("decomposeCells");
	# endif
    
	const faceListPMG& faces = mesh_.faces();
	const labelList& owner = mesh_.owner();
	const labelList& neighbour = mesh_.neighbour();
	const label nIntFaces = mesh_.nInternalFaces();
	for(label faceI=0;faceI<nIntFaces;++faceI)
	{
		const face& f = faces[faceI];
		
		forAll(f, pI)
        {
            if( f[pI] >= concaveVertex.size() )
                continue;
            
			if( concaveVertex[f[pI]] )
			{
				decomposeCell_[owner[faceI]] = true;
				decomposeCell_[neighbour[faceI]] = true;
				
				# ifdef DEBUGDec
				mesh_.addCellToSubset(id, owner[faceI]);
				mesh_.addCellToSubset(id, neighbour[faceI]);
				# endif
				
				break;
			}
        }
	}
	
	# ifdef DEBUGDec
    Info << "Writting mesh" << endl;
	mesh_.write();
	::exit(1);
	# endif
}	

void decomposeCellsNearConcaveEdges::decomposeConcaveCells()
{
	decomposeCells dc(mesh_);
	dc.decomposeMesh(decomposeCell_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
