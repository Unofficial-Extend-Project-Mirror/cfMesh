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

#include "voronoiMeshExtractor.H"
#include "tessellationElement.H"
#include "demandDrivenData.H"

#define DEBUGVoronoi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createAddressing() const
{
	if( pointEdgesPtr_ )
		return;	
	pointEdgesPtr_ = new VRWGraph(tetCreator_.tetPoints().size());
	VRWGraph& pointEdges = *pointEdgesPtr_;
	
	if( edgeTetsPtr_ )
		return;
	edgeTetsPtr_ = new VRWGraph();
	VRWGraph& edgeTets = *edgeTetsPtr_;
	
	if( boundaryEdgePtr_ )
		return;
	boundaryEdgePtr_ = new boolList();
	boolList& boundaryEdges = *boundaryEdgePtr_;
	
	if( edgesPtr_ )
		return;
	edgesPtr_ = new LongList<edge>();
	LongList<edge>& edges = *edgesPtr_;
	
	//- create edges, pointEdges and edgeTets
	const LongList<partTet>& tets = tetCreator_.tets();
	forAll(tets, tetI)
	{
		const edgeList tetEdges = tets[tetI].edges();
		
		forAll(tetEdges, eI)
		{
			const edge& e = tetEdges[eI];
			const label start = e.start();
			
			bool store(true);
			
			forAllRow(pointEdges, start, peI)
			{
				const label edgeI = pointEdges(start, peI);
				
				if( edges[edgeI] == e )
				{
					store = false;
					
					edgeTets.append(edgeI, tetI);
				}
			}
			
			if( store )
			{
				pointEdges.append(start, edges.size());
				pointEdges.append(e.end(), edges.size());
			
				FixedList<label, 1> helper;
				helper[0] = tetI;
				edgeTets.appendList(helper);
				
				edges.append(e);
			}
		}
	}
	
	pointEdges.optimizeMemoryUsage();
	edgeTets.optimizeMemoryUsage();
	
	//- sort edge-tets in circular order
	boundaryEdges.setSize(edgeTets.size());
	boundaryEdges = false;
	
	forAll(edgeTets, edgeI)
	{
		labelList eTets(edgeTets.sizeOfRow(edgeI));
		forAll(eTets, etI)
			eTets[etI] = edgeTets(edgeI, etI);
		
		labelList nFound(eTets.size(), 0);
		List<FixedList<label, 2> > elNeighbours(eTets.size());
		forAll(elNeighbours, tetI)
			elNeighbours[tetI] = -1;
		
		forAll(eTets, tetI)
		{
			if( nFound[tetI] == 2 )
				continue;
			
			const partTet& pt = tets[eTets[tetI]];
			
			for(label i=tetI+1;i<eTets.size();++i)
			{
				const partTet& ptNei = tets[eTets[i]];
				
				label nShared(0);
				for(label j=0;j<4;++j)
					for(label k=0;k<4;++k)
						if( pt[j] == ptNei[k] )
						{
							++nShared;
							break;
						}
				
				if( nShared == 3 )
				{
					elNeighbours[tetI][nFound[tetI]++] = i;
					elNeighbours[i][nFound[i]++] = tetI;
				}
			}
		}
		
		bool sort(true);
		forAll(nFound, tetI)
			if( nFound[tetI] != 2 )
			{
				boundaryEdges[edgeI] = true;
				sort = false;
				break;
			}

		if( sort )
		{
			labelList sortedTets(eTets.size());
			
			label currI(0), nextI(elNeighbours[0][0]), counter(0);
			sortedTets[counter++] = eTets[currI];

			do
			{
				sortedTets[counter++] = eTets[nextI];
				
				//- find the next element
				if( elNeighbours[nextI][0] == currI )
				{
					currI = nextI;
					nextI = elNeighbours[nextI][1];
				}
				else
				{
					currI = nextI;
					nextI = elNeighbours[nextI][0];
				}
			} while( counter < eTets.size() );
			
			edgeTets.setRow(edgeI, sortedTets);
		}
	}
}

const VRWGraph& voronoiMeshExtractor::pointEdges() const
{
	if( !pointEdgesPtr_ )
		createAddressing();
	
	return *pointEdgesPtr_;
}

const LongList<edge>& voronoiMeshExtractor::edges() const
{
	if( !edgesPtr_ )
		createAddressing();
	
	return *edgesPtr_;
}

const VRWGraph& voronoiMeshExtractor::edgeTets() const
{
	if( !edgeTetsPtr_ )
		createAddressing();
	
	return *edgeTetsPtr_;
}

const boolList& voronoiMeshExtractor::boundaryEdge() const
{
	if( !boundaryEdgePtr_ )
		createAddressing();
	
	return *boundaryEdgePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
