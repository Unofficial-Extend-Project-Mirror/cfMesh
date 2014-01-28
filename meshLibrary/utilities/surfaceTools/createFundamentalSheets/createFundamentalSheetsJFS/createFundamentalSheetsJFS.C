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

#include "createFundamentalSheetsJFS.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "extrudeLayer.H"

#include "addToRunTimeSelectionTable.H"

#include <omp.h>

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(createFundamentalSheetsJFS, 0);
addToRunTimeSelectionTable
(
    createFundamentalSheets,
    createFundamentalSheetsJFS,
    polyMeshGen
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void createFundamentalSheetsJFS::createInitialSheet()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    const label start = boundaries[0].patchStart();
    const label end
    (
        boundaries[boundaries.size()-1].patchStart() +
        boundaries[boundaries.size()-1].patchSize()
    );

    faceListPMG::subList bFaces(mesh_.faces(), end-start, start);

    const labelList& owner = mesh_.owner();

    LongList<labelPair> extrudeFaces(end-start);

    # ifdef USE_OMP
    # pragma omp parallel for
    # endif
    for(label faceI=start;faceI<end;++faceI)
        extrudeFaces[faceI-start] = labelPair(faceI, owner[faceI]);

    extrudeLayer(mesh_, extrudeFaces);
}

void createFundamentalSheetsJFS::createSheetsAtFeatureEdges()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    const label start = boundaries[0].patchStart();
    const label end
    (
        boundaries[boundaries.size()-1].patchStart() +
        boundaries[boundaries.size()-1].patchSize()
    );

    faceListPMG::subList bFaces(mesh_.faces(), end-start, start);

    labelList patchCell(mesh_.cells().size());

    LongList<labelPair> front;

    # ifdef USE_OMP
    const label nThreads = 2 * omp_get_num_procs();
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for
        # endif
        forAll(patchCell, cellI)
            patchCell[cellI] = -1;

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for
        # endif
        for(label faceI=start;faceI<end;++faceI)
            patchCell[owner[faceI]] = mesh_.faceIsInPatch(faceI);

        //- create the front faces
        LongList<labelPair> localFront;

        # ifdef USE_OMP
        # pragma omp for
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const cell& c = cells[owner[faceI]];
            const label patchI = mesh_.faceIsInPatch(faceI);

            forAll(c, fI)
            {
                if( neighbour[c[fI]] < 0 )
                    continue;

                label nei = owner[c[fI]];
                if( nei == owner[faceI] )
                    nei = neighbour[c[fI]];

                if( patchCell[nei] != patchI )
                    localFront.append(labelPair(c[fI], owner[faceI]));
            }
        }

        label frontStart(-1);
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            frontStart = front.size();
            front.setSize(front.size()+localFront.size());
        }

        //- copy the local front into the global front
        forAll(localFront, lfI)
            front[frontStart+lfI] = localFront[lfI];
    }

    //- extrude the layer
    extrudeLayer(mesh_, front);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, octree, regions for boundary vertices
createFundamentalSheetsJFS::createFundamentalSheetsJFS
(
    polyMeshGen& mesh
)
:
    createFundamentalSheets(mesh)
{
    createInitialSheet();

    createSheetsAtFeatureEdges();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

createFundamentalSheetsJFS::~createFundamentalSheetsJFS()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
