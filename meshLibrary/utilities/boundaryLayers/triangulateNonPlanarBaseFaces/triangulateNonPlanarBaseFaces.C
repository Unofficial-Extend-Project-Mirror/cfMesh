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

#include "triangulateNonPlanarBaseFaces.H"
#include "dictionary.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

<<<<<<< .mine=======label voronoiMeshExtractor::sameOrientation_[6] = {3, 1, 2, 2, 3, 0};

label voronoiMeshExtractor::oppositeOrientation_[6] = {2, 3, 1, 0, 0, 1};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::clearOut()
{
    deleteDemandDrivenData(pointEdgesPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(edgeTetsPtr_);
    deleteDemandDrivenData(boundaryEdgePtr_);
}

>>>>>>> .theirs// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::triangulateNonPlanarBaseFaces
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
<<<<<<< .mine    invertedCell_(mesh_.cells().size(), false),
    decomposeFace_(mesh_.faces().size(), false),
    tol_(0.5)
{}
=======    pointEdgesPtr_(NULL),
    edgesPtr_(NULL),
    edgeTetsPtr_(NULL),
    boundaryEdgePtr_(NULL)
{}
>>>>>>> .theirs
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::~triangulateNonPlanarBaseFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateNonPlanarBaseFaces::setRelativeTolerance(const scalar tol)
{
    tol_ = tol;
}

void triangulateNonPlanarBaseFaces::triangulateLayers()
{
    if( findNonPlanarBoundaryFaces() )
    {
        Info << "Decomposing twisted boundary faces" << endl;

        decomposeBoundaryFaces();

        decomposeCellsIntoPyramids();
    }
    else
    {
        Info << "All boundary faces are flat" << endl;
    }
}

void triangulateNonPlanarBaseFaces::readSettings
(
    const dictionary& meshDict,
    triangulateNonPlanarBaseFaces& triangulator
)
{
<<<<<<< .mine    if( meshDict.found("boundaryLayers") )
    {
        const dictionary& layersDict = meshDict.subDict("boundaryLayers");
=======    Info << "Extracting voronoi mesh" << endl;

    //- copy tet points into the mesh
    createPoints();
>>>>>>> .theirs
<<<<<<< .mine        if( layersDict.found("optimisationParameters") )
        {
            const dictionary& optLayerDict =
                layersDict.subDict("optimisationParameters");

            if( optLayerDict.found("relFlatnessTol") )
            {
                const scalar relTol =
                    readScalar(optLayerDict.lookup("relFlatnessTol"));

                triangulator.setRelativeTolerance(relTol);
            }
        }
    }
=======    //- create the mesh
    createPolyMesh();

    polyMeshGenModifier(mesh_).reorderBoundaryFaces();
    polyMeshGenModifier(mesh_).removeUnusedVertices();

    Info << "Mesh has :" << nl
        << mesh_.points().size() << " vertices " << nl
        << mesh_.faces().size() << " faces" << nl
        << mesh_.cells().size() << " cells" << endl;

    Info << "Finished extracting voronoi mesh" << endl;
>>>>>>> .theirs}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
