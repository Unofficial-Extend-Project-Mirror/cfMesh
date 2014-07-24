/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "symmetryPlaneOptimisation.H"
#include "partTetMesh.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"

#include "partTetMeshSimplex.H"
#include "meshUntangler.H"
#include "volumeOptimizer.H"
#include "knuppMetric.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void symmetryPlaneOptimisation::detectSymmetryPlanes()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    symmetryPlanes_.clear();

    typedef std::map<label, std::pair<vector, label> > mapType;
    mapType centreSum, normalSum;

    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].type() == "symmetryPlane" )
        {
            std::pair<vector, label>& cs = centreSum[patchI];
            cs = std::pair<vector, label>(vector::zero, 0);

            std::pair<vector, label>& ns = normalSum[patchI];
            ns = std::pair<vector, label>(vector::zero, 0);

            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for(label faceI=start;faceI<end;++faceI)
            {
                cs.first += faces[faceI].centre(points);
                ns.first += faces[faceI].normal(points);
            }

            cs.second = ns.second = boundaries[patchI].patchSize();
        }
    }

    if( Pstream::parRun() )
    {
        //- sum up all normals and centres of all processors
        //- every symmetry plane patch must be present on all processors
        forAllIter(mapType, centreSum, pIter)
        {
            std::pair<vector, label>& cs = pIter->second;
            reduce(cs.second, sumOp<label>());
            reduce(cs.first, sumOp<vector>());

            std::pair<vector, label>& ns = normalSum[pIter->first];
            reduce(ns.first, sumOp<vector>());
            ns.second = cs.second;
        }
    }

    //- create planes corresponding to each symmetry plane
    forAllConstIter(mapType, centreSum, it)
    {
        const point c = it->second.first / it->second.second;

        const std::pair<vector, label>& ns = normalSum[it->first];
        const point n = ns.first / ns.second;

        const word pName = boundaries[it->first].patchName();
        symmetryPlanes_.insert(std::make_pair(pName, plane(c, n)));
    }
}

void symmetryPlaneOptimisation::pointInPlanes(VRWGraph&) const
{

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
symmetryPlaneOptimisation::symmetryPlaneOptimisation(polyMeshGen& mesh)
:
    mesh_(mesh)
{
    detectSymmetryPlanes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

symmetryPlaneOptimisation::~symmetryPlaneOptimisation()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void symmetryPlaneOptimisation::optimizeSymmetryPlanes()
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
