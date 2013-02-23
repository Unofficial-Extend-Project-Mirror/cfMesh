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

#include "triSurf.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurf::calculateFaceGroups() const
{
	nFaceGroups_ = 0;
	
	faceGroupPtr_ = new labelListPMG(triSurface::size(), -1);
	labelListPMG& faceGroup = *faceGroupPtr_;
	
	const labelListList& faceEdges = this->faceEdges();
	const labelListList& edgeFaces = this->edgeFaces();
	
	labelListPMG front;
	
	forAll(faceGroup, fI)
	{
		if( faceGroup[fI] != -1 )
			continue;
		
		front.clear();
		front.append(fI);
		faceGroup[fI] = nFaceGroups_;
		
		while( front.size() != 0 )
		{
			const label fLabel = front.removeLastElement();
			
			const labelList& fEdges = faceEdges[fLabel];
			forAll(fEdges, feI)
			{
				const label eI = fEdges[feI];
				
				if( edgeFaces[eI].size() != 2 )
					continue;
				
				label nei = edgeFaces[eI][0];
				if( nei == fLabel )
					nei = edgeFaces[eI][1];
				
				if( faceGroup[nei] == -1 )
				{
					faceGroup[nei] = nFaceGroups_;
					front.append(nei);
				}
			}
		}
		
		++nFaceGroups_;
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurf::triSurf(const fileName& name)
:
	triSurface(name),
	triSubsets_(),
	nFaceGroups_(-1),
	faceGroupPtr_(NULL)
{
}

triSurf::triSurf
(
	const List<labelledTri>& triangles,
	const geometricSurfacePatchList& patches,
	const pointField& points,
	std::map<word, labelListPMG>& faceSubsets
)
:
	triSurface(triangles, patches, points),
	triSubsets_(),
	nFaceGroups_(-1),
	faceGroupPtr_(NULL)
{
	std::map<word, labelListPMG>::const_iterator iter;
	for(iter=faceSubsets.begin();iter!=faceSubsets.end();++iter)
	{
		triSubsets_.insert
		(
			std::pair<word, labelListPMG>(iter->first, iter->second)
		);
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurf::~triSurf()
{
}
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurf::readFaceSubsets(const fileName& fName)
{
	IFstream file(fName);
	
	//- read the number of subsets
	label nSubsets;
	file >> nSubsets;
	
	// read subsets
	for(label subsetI=0;subsetI<nSubsets;++subsetI)
	{
		word name;
		file >> name;
		labelListPMG set;
		file >> set;
		
		this->addFacetsToSubset(name, set);
	}
}

void triSurf::writeFaceSubsets(const fileName& fName) const
{
	OFstream file(fName);
	
	label nSubsets(0);
	
	//- find the number of subsets
	std::map<word, labelListPMG>::const_iterator iter;
	for(iter = triSubsets_.begin();iter != triSubsets_.end();++iter)
		++nSubsets;
	file << nSubsets << nl;
	
	//- write subsets
	for(iter = triSubsets_.begin();iter != triSubsets_.end();++iter)
	{
		file << iter->first << nl;
		file << iter->second << nl;
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
