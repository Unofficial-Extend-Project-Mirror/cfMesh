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

#include "removeCellsInSelectedDomains.H"
#include "meshOctree.H"
#include "findCellsIntersectingSurface.H"
#include "triSurf.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const VRWGraph& removeCellsInSelectedDomains::cellsIntersectedBySurfaceFacets()
{
    if( !surfIntersectionPtr_ )
    {
        Info << "Calculating intersections" << endl;
        surfIntersectionPtr_ = new findCellsIntersectingSurface(mesh_, octree_);
    }

    Info << "Intersection pointer " << surfIntersectionPtr_ << endl;

    return surfIntersectionPtr_->facetsIntersectingCells();
}

void removeCellsInSelectedDomains::markSelectedFacets()
{
    const triSurf& surf = octree_.surface();
    const geometricSurfacePatchList& patches = surf.patches();

    boolList usedPatch(patches.size(), false);

    selectedFacets_.setSize(surf.size());
    selectedFacets_ = false;

    forAll(domains_, domainI)
    {
        const DynList<word, 4>& domain = domains_[domainI];

        forAll(domain, subsetI)
        {
            const word& sName = domain[subsetI];

            label index = surf.facetSubsetIndex(sName);

            if( index >= 0 )
            {
                //- it is a subset
                labelListPMG facetsInSubset;
                surf.facetsInSubset(index, facetsInSubset);

                forAll(facetsInSubset, i)
                    selectedFacets_[facetsInSubset[i]] = true;
            }
            else
            {
                //- check if it is a patch
                forAll(patches, i)
                {
                    if( patches[i].name() == sName )
                    {
                        usedPatch[i] = true;
                        break;
                    }
                }
            }
        }
    }

    //- selected facets in used patches
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];

        if( usedPatch[tri.region()] )
            selectedFacets_[triI] = true;
    }
}

void removeCellsInSelectedDomains::findLeavesInsideRegions()
{
    insideLeaves_.setSize(octree_.numberOfLeaves());
    insideLeaves_ = false;

    //- find octree leaves intersected by the selected domains
    boolList intersectedLeaves(insideLeaves_.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(intersectedLeaves, leafI)
    {
        DynList<label> trianglesInLeaf;
        octree_.containedTriangles(leafI, trianglesInLeaf);

        forAll(trianglesInLeaf, i)
        {
            if( selectedFacets_[trianglesInLeaf[i]] )
            {
                intersectedLeaves[leafI] = true;
                break;
            }
        }
    }

    //- find islands of octree leaves which are not intersected by
    //- the selected domains

    //- select the islands which do not have contact with the outer bondary
}

void removeCellsInSelectedDomains::findAndRemoveCells()
{
    //- find the cells which are intersected by the selected domains
    boolList intersectedCells(mesh_.cells().size(), false);

    const VRWGraph& facetsIntersectingCell = cellsIntersectedBySurfaceFacets();

    Info << "1.Here" << endl;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(facetsIntersectingCell, cellI)
    {
        forAllRow(facetsIntersectingCell, cellI, triI)
        {
            if( selectedFacets_[facetsIntersectingCell(cellI, triI)] )
            {
                intersectedCells[cellI] = true;
                break;
            }
        }
    }

    //- TODO: implement the group marking algorithm properly using templates
    //- find islands of cells which are not intersected by the selected domains
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    labelListPMG findCellGroups(intersectedCells.size(), -1);
    label nGroups(0);

    Info << "2. Here" << endl;

    forAll(findCellGroups, cellI)
    {
        if( findCellGroups[cellI] != -1 )
            continue;
        if( intersectedCells[cellI] )
            continue;

        findCellGroups[cellI] = nGroups;

        labelListPMG front;
        front.append(cellI);

        while( front.size() )
        {
            const label cLabel = front.removeLastElement();

            const cell& c = cells[cLabel];

            //- find neighbouring elements
            DynList<label> neis;
            forAll(c, fI)
            {
                const label faceI = c[fI];

                if( owner[faceI] == cellI && neighbour[faceI] >= 0 )
                {
                    neis.append(neighbour[faceI]);
                }
                else
                {
                    neis.append(owner[faceI]);
                }
            }

            //- check which neighbours pass the criteria
            forAll(neis, i)
            {
                const label nei = neis[i];

                if( findCellGroups[nei] != -1 )
                    continue;
                if( intersectedCells[nei] )
                    continue;
            }
        }

        ++nGroups;
    }

    //- check which groups have cells at the outer boundary
    //- which are not intersected by the selected domains
    Info << "3. Here" << endl;

    boolList outerGroup(nGroups, false);
    const PtrList<writePatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            if( intersectedCells[owner[faceI]] )
                continue;

            outerGroup[findCellGroups[owner[faceI]]] = true;
        }
    }

    //- check which islands of cells correspond to octree boxes
    //- marked as internal boxes
    Info << "4. Here" << endl;
    boolList internalCells(intersectedCells.size(), false);

    forAll(findCellGroups, cellI)
    {
        const label groupI = findCellGroups[cellI];

        if( outerGroup[groupI] )
            continue;

        internalCells[cellI] = true;
    }

    //- remove cells inside the selected domains
    Info << "5. Here" << endl;
    polyMeshGenModifier(mesh_).removeCells(internalCells);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

removeCellsInSelectedDomains::removeCellsInSelectedDomains
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octree_(octree),
    domains_(),
    selectedFacets_(),
    surfIntersectionPtr_(NULL),
    insideLeaves_(octree_.numberOfLeaves(), false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

removeCellsInSelectedDomains::~removeCellsInSelectedDomains()
{
    deleteDemandDrivenData(surfIntersectionPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void removeCellsInSelectedDomains::selectCellsInDomain(const wordList& domain)
{
    const label s = domains_.size();

    domains_.setSize(s+1);

    domains_[s] = domain;
}

void removeCellsInSelectedDomains::removeCells()
{
    Info << "Starting removing cells in user-selected domains" << endl;

    markSelectedFacets();

    //findLeavesInsideRegions();

    Info << "Finding and removing cells" << endl;

    findAndRemoveCells();

    Info << "Finished removing cells in user-selected domains" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
