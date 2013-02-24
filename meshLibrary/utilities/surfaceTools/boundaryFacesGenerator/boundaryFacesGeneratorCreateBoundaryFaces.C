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

#include "boundaryFacesGenerator.H"
#include "sortEdgesIntoChains.H"

//#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList boundaryFacesGenerator::patchesForPoint(const label vI) const
{
    labelList pcs(pointRegions_.sizeOfRow(vI));

    forAll(pcs, pI)
        pcs[pI] = pointRegions_(vI, pI);

    return pcs;
}

void boundaryFacesGenerator::determineCellPatches
(
    const label cellI,
    DynList<label>& patches
)
{

    const cellListPMG& polyCells_ = mesh_.cells();
    const faceListPMG& polyFaces_ = mesh_.faces();
    const cell& c = polyCells_[cellI];

    # ifdef DEBUGCutter
    Info << "Cell faces are " << c << endl;
    forAll(c, fI)
        Info << "Face " << c[fI] << " is "
            << polyFaces_[c[fI]] << endl;
    # endif

    const labelList cellPoints = c.labels(polyFaces_);

    # ifdef DEBUGCutter
    Info << "Cell points are " << cellPoints << endl;
    # endif

    //- find which boundary regions intersect the polyhedron
    //- this can be determined by monitoring patch labels
    //- of cell points
    patches.clear();

    forAll(cellPoints, cpI)
    {
        const labelList pp = patchesForPoint(cellPoints[cpI]);
        forAll(pp, ppI)
            patches.appendIfNotIn(pp[ppI]);
    }

    # ifdef DEBUGCutter
    Info << "Patches for cell " << cellI << " are " << patches << endl;
    # endif
}

void boundaryFacesGenerator::createBoundaryChainsForCell
(
    const label cellI,
    labelListList& boundaryChains
)
{
    const faceListPMG& polyFaces_ = mesh_.faces();
    const cell& c = mesh_.cells()[cellI];
    //- find edges which appear once
    //- these edges are on the boundary
    const edgeList cellEdges = c.edges(polyFaces_);

    labelList nAppearances(cellEdges.size(), 0);
    forAll(c, fI)
    {
        const face& f = polyFaces_[c[fI]];

        const edgeList fe = f.edges();

        forAll(fe, eI)
            forAll(cellEdges, i)
            if( fe[eI] == cellEdges[i] )
            {
                ++nAppearances[i];
                break;
            }
    }

    //- create a list of boundary edges
    DynList<edge> boundaryEdges;

    forAll(nAppearances, aI)
        if( nAppearances[aI] == 1 )
        {
            boundaryEdges.append(cellEdges[aI]);
        }
        else if( nAppearances[aI] > 2 )
        {
            const cellListPMG& polyCells_ = mesh_.cells();
            forAll(polyCells_[cellI], fI)
                Info << "Face " << polyCells_[cellI][fI] << " is "
                    << polyFaces_[polyCells_[cellI][fI]] << endl;
            FatalErrorIn
            (
                "void createBoundaryFaces()"
            ) << "Cell " << cellI << " is not topologically closed!"
                << abort(FatalError);
        }

    //- sort edges into chains
    boundaryChains = sortEdgesIntoChains(boundaryEdges).sortedChains();

    # ifdef DEBUGCutter
    Info << "boundaryChains " << boundaryChains << endl;
    # endif
}

void boundaryFacesGenerator::createBoundaryFacesForCell(const label cellI)
{
    //- check if there are any faces in the current cell
    if( nFacesInCell_[cellI] == 0 )
        FatalErrorIn
        (
            "void boundaryFacesGenerator::createBoundaryFacesForCell"
            "(const label)"
        ) << "Cell " << cellI << " has no faces!!"
            << nl << "This generally means that there is not enough cells "
            << "in that region. It can be remedied by refining "
            << "the mesh there." << exit(FatalError);

    //- determine cell patches
    DynList<label> patches;
    determineCellPatches(cellI, patches);

    # ifdef DEBUGCutter
    Info << "Cell patches " << patches << endl;
    # endif

    //- create boundary chains for cell
    labelListList boundaryChains;
    createBoundaryChainsForCell(cellI, boundaryChains);

    //- create storage for faces created from the boundary chain
    //- this list contains faces created from the chain
    //- with respect to the boundary patch they belong to
    //- first list : faces for each chain
    //- second list: faces with respect to the patch
    //- third list : there can be more than one face per patch
    List< List< DynList<face> > > facePoints
    (
        boundaryChains.size()
    );

    forAll(facePoints, chI)
    {
        facePoints[chI].setSize(patches.size()+1);
    }

    //- create faces for the given chains
    forAll(boundaryChains, bcI)
        createFacesForChain
        (
            boundaryChains[bcI],
            patches,
            facePoints[bcI]
        );

    # ifdef DEBUGCutter
    Info << "facePoints " << facePoints << endl;
    # endif

    //- store faces into the list
    forAll(facePoints, chI)
    {
        const List< DynList<face> >& chfcs = facePoints[chI];

        forAll(chfcs, patchI)
        {
            const DynList<face>& patchFaces = chfcs[patchI];

            if( patchI < patches.size() )
            {
                forAll(patchFaces, pfI)
                {
                    boundaryFaces_.appendList(patchFaces[pfI]);
                    boundaryOwners_.append(cellI);
                    boundaryPatches_.append(patches[patchI]);
                }
            }
            else
            {
                # ifdef DEBUGCutter
                Info << "Faces " << patchFaces
                    << " are in default patch" << endl;
                # endif
                //- faces in default patch
                forAll(patchFaces, pfI)
                {
                    hasFacesInDefaultPatch_ = true;
                    boundaryFaces_.appendList(patchFaces[pfI]);
                    boundaryOwners_.append(cellI);
                    boundaryPatches_.append(defaultPatchID_);
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
