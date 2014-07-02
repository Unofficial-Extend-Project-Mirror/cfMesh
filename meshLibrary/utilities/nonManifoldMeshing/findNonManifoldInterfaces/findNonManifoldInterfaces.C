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

#include "findNonManifoldInterfaces.H"
#include "meshOctree.H"
#include "findCellsIntersectingSurface.H"
#include "triSurf.H"
#include "demandDrivenData.H"

#define DEBUGNonManifoldInterfaces

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void findNonManifoldInterfaces::findIntersectedCells()
{
    Info << "Starting finding cells intersected by the surface" << endl;

    findCellsIntersectingSurface cis(mesh_, octree_);

    const VRWGraph& facetsIntersectedCells = cis.facetsIntersectingCells();

    //- find which surface patches intersect a cell
    const triSurf& surface = octree_.surface();
    patchesIntersectingCells_.clear();
    patchesIntersectingCells_.setSize(facetsIntersectedCells.size());

    forAll(facetsIntersectedCells, cellI)
    {
        forAllRow(facetsIntersectedCells, cellI, i)
        {
            const label triI = facetsIntersectedCells(cellI, i);
            patchesIntersectingCells_.appendIfNotIn
            (
                cellI,
                surface[triI].region()
            );
        }
    }

    # ifdef DEBUGNonManifoldInterfaces
    Map<label> regionToID;
    forAll(patchesIntersectingCells_, cellI)
    {
        forAllRow(patchesIntersectingCells_, cellI, i)
        {
            const label regionI = patchesIntersectingCells_(cellI, i);
            Map<label>::iterator iter = regionToID.find(regionI);
            if( iter == regionToID.end() )
            {
                const word sName =
                    "intersected_"+mesh_.boundaries()[regionI].patchName();
                regionToID.insert
                (
                    regionI,
                    mesh_.addCellSubset(sName)
                );

                iter = regionToID.find(regionI);
            }

            mesh_.addCellToSubset(iter(), cellI);
        }
    }
    # endif

    Info << "Finished finding cells intersected by the surface" << endl;
}

void findNonManifoldInterfaces::findNeighbouringGroups()
{
    //- start assigning neighbouring cells of the intersected cells into groups
    neiGroups_.setSize(0);
    neiGroups_.setSize(patchesIntersectingCells_.size());
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    for(label fI=0;fI<mesh_.nInternalFaces();++fI)
    {
        const label own = owner[fI];
        const label nei = neighbour[fI];

        if( patchesIntersectingCells_.sizeOfRow(own) )
        {
            forAllRow(patchesIntersectingCells_, own, i)
            {
                const label patchI = patchesIntersectingCells_(own, i);

                if( !patchesIntersectingCells_.contains(nei, patchI) )
                    neiGroups_.append(nei, patchI);
            }
        }
        else if( patchesIntersectingCells_.sizeOfRow(nei) )
        {
            forAllRow(patchesIntersectingCells_, nei, i)
            {
                const label patchI = patchesIntersectingCells_(nei, i);

                if( !patchesIntersectingCells_.contains(own, patchI) )
                    neiGroups_.append(own, patchI);
            }
        }
    }

    # ifdef DEBUGNonManifoldInterfaces
    Map<label> regionToID;
    forAll(neiGroups_, cellI)
    {
        forAllRow(neiGroups_, cellI, i)
        {
            const label regionI = neiGroups_(cellI, i);
            Map<label>::iterator iter = regionToID.find(regionI);
            if( iter == regionToID.end() )
            {
                const word sName =
                    "neighboursOfPatch_" +
                    mesh_.boundaries()[regionI].patchName();

                regionToID.insert
                (
                    regionI,
                    mesh_.addCellSubset(sName)
                );

                iter = regionToID.find(regionI);
            }

            mesh_.addCellToSubset(iter(), cellI);
        }
    }
    # endif
}

bool findNonManifoldInterfaces::findFaceCandidates()
{
    internalFacePatches_.clear();
    internalFacePatches_.setSize(mesh_.nInternalFaces());

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    forAll(internalFacePatches_, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        forAllRow(patchesIntersectingCells_, own , i)
        {
            const label patchI = patchesIntersectingCells_(own, i);

            if(
                !patchesIntersectingCells_.contains(nei, patchI) &&
                neiGroups_.contains(nei, patchI)
            )
                internalFacePatches_.appendIfNotIn(faceI, patchI);
        }

        forAllRow(patchesIntersectingCells_, nei , i)
        {
            const label patchI = patchesIntersectingCells_(nei, i);

            if(
                !patchesIntersectingCells_.contains(own, patchI) &&
                neiGroups_.contains(own, patchI)
            )
                internalFacePatches_.appendIfNotIn(faceI, patchI);
        }
    }

    # ifdef DEBUGNonManifoldInterfaces
    Map<label> regionToID;
    forAll(internalFacePatches_, faceI)
    {
        forAllRow(internalFacePatches_, faceI, i)
        {
            const label regionI = internalFacePatches_(faceI, i);
            Map<label>::iterator iter = regionToID.find(regionI);
            if( iter == regionToID.end() )
            {
                const word sName =
                    "interfaceFaceToPatch_" +
                    mesh_.boundaries()[regionI].patchName();

                regionToID.insert
                (
                    regionI,
                    mesh_.addFaceSubset(sName)
                );

                iter = regionToID.find(regionI);
            }

            mesh_.addFaceToSubset(iter(), faceI);
        }
    }

    const label problemFacesID = mesh_.addFaceSubset("problematicFaces");
    # endif

    //- check if there exist any proximity problems
    forAll(internalFacePatches_, faceI)
    {
        if( internalFacePatches_.sizeOfRow(faceI) < 2 )
            continue;

        //- check if it is a proximity problem
        Info << "More than one internal interface assigned to face "
            << faceI << endl;

        # ifdef DEBUGNonManifoldInterfaces
        mesh_.addFaceToSubset(problemFacesID, faceI);
        # endif
    }

    return false;
}

bool findNonManifoldInterfaces::extractInterfaces()
{

    //- the procedure ended without any problems
    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and octree
findNonManifoldInterfaces::findNonManifoldInterfaces
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octree_(octree),
    patchesIntersectingCells_(),
    neiGroups_(),
    internalFacePatches_()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

findNonManifoldInterfaces::~findNonManifoldInterfaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void findNonManifoldInterfaces::createNonManifoldInterfaces()
{
    do
    {
        //- find cells intersected by the surface
        findIntersectedCells();

        //- find cells near the intersected cells
        findNeighbouringGroups();
    } while( findFaceCandidates() );

    if( extractInterfaces() )
        WarningIn
        (
            "void findNonManifoldInterfaces::createNonManifoldInterfaces()"
        ) << "Could not generate all internal interfaces" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
