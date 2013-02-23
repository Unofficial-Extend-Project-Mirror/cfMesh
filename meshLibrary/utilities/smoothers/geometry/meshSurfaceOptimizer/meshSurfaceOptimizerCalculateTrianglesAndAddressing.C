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
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngine.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::calculateTrianglesAndAddressing() const
{
	if( trianglesPtr_ || pointTrianglesPtr_ )
		FatalErrorIn
		(
			"void meshSurfaceOptimizer::calculateTrianglesAndAddressing() const"
		) << "Addressing is already calculated!" << abort(FatalError);
	
	const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
	const labelList& bp = surfaceEngine_.bp();
	
	triFace triangle; //- helper
	
	pointTrianglesPtr_ = new VRWGraph(surfaceEngine_.boundaryPoints().size());
	VRWGraph& pointTriangles = *pointTrianglesPtr_;
	
	trianglesPtr_ = new LongList<triFace>();
	LongList<triFace>& triangles = *trianglesPtr_;
	
	//- start creating triangles
	forAll(bFaces, bfI)
	{
		const face& bf = bFaces[bfI];
		
		forAll(bf, pI)
		{
            const label nTrias = bf.size() - 2;
			triangle[0] = bp[bf[pI]];
            
            for(label i=0;i<nTrias;++i)
            {
                triangle[1] = bp[bf[(pI+i+1)%bf.size()]];
                triangle[2] = bp[bf[(pI+i+2)%bf.size()]];
			
                triangles.append(triangle);
            }
		}
	}
    
    //- create point-triangles addressing
    labelList nTriaAtPoint(pointTriangles.size(), 0);
    forAll(triangles, triI)
    {
        const triFace& tria = triangles[triI];
        ++nTriaAtPoint[tria[0]];
    }
    
    forAll(nTriaAtPoint, bpI)
        pointTriangles.setRowSize(bpI, nTriaAtPoint[bpI]);
    nTriaAtPoint = 0;
    
    forAll(triangles, triI)
    {
        const triFace& tria = triangles[triI];
        pointTriangles(tria[0], nTriaAtPoint[tria[0]]++) = triI;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
