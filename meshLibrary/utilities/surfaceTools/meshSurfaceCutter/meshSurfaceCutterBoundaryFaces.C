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

#include "meshSurfaceCutter.H"
#include "boundaryFacesGenerator.H"

//#define DEBUG_meshSurfaceCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList meshSurfaceCutter::patchesForPoint(const label vI) const
{
    if( pointTriIndex_[vI].size() )
    {
        labelList pcs(pointTriIndex_[vI].size());
        const List<labelledTri>& fcs = surface_.localFaces();

        const DynList<label>& tria = pointTriIndex_[vI];
        short i(0);
		forAll(tria, tI)
            pcs[i++] = fcs[tria[tI]].region();

        return pcs;
    }
    else
    {
        return labelList();
    }
}

void meshSurfaceCutter::createBoundaryFaces()
{
	VRWGraph pointPatches(nPoints_);
	forAll(pointPatches, pI)
	{
		const labelList pp = patchesForPoint(pI);
		forAll(pp, ppI)
			pointPatches.append(pI, pp[ppI]);
	}
	
	boundaryFacesGenerator
	(
		mesh_,
		nFacesInCell_,
		nPoints_,
		nIntFaces_,
		boundaryCell_,
		pointPatches,
		surface_
	);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
