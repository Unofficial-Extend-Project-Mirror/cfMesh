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

#include "triangulateNonPlanarBoundaryFaces.H"
#include "decomposeCells.H"
#include "faceDecomposition.H"
#include "triangle.H"

//#define DEBUGDecompose

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateNonPlanarBoundaryFaces::findBoundaryCellsToDecompose()
{
    const pointFieldPMG& points = mesh_.points();
    const labelList& owner = mesh_.owner();
    const faceListPMG& faces = mesh_.faces();
    
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    
    label nDecompose(0);
    
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();
        
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];
            
            forAll(f, pI)
                pointPatches_[f[pI]].appendIfNotIn(patchI);
            
            //- do not waste time on triangles
            if( f.size() < 4 ) continue;
        
            if( !faceDecomposition(f, points).isFacePlanar() )
            {
                ++nDecompose;
                decomposeCell_[owner[faceI]] = true;
            }
        }
    }
    
    Info << "Found " << nDecompose << " non-planar boundary faces" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triangulateNonPlanarBoundaryFaces::triangulateNonPlanarBoundaryFaces
(
    polyMeshGen& mesh
)
    :
    mesh_(mesh),
    decomposeCell_(mesh.cells().size(), false),
    pointPatches_(mesh.points().size())
{
}

triangulateNonPlanarBoundaryFaces::~triangulateNonPlanarBoundaryFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateNonPlanarBoundaryFaces::decomposeMesh()
{
    FatalError << "Not implemented" << exit(FatalError);
    
    Info << "Decomposing non-planar boundary faces" << endl;
    findBoundaryCellsToDecompose();
    
    decomposeCells decomposeNonPlanarCells(mesh_);
    decomposeNonPlanarCells.decomposeMesh(decomposeCell_);
    
    Info << "Finished decomposing non-planar boundary faces" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
