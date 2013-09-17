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

#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const VRWGraph& removeCellsInSelectedDomains::cellsIntersectedBySurfaceFacets()
{
    if( !surfIntersectionPtr_ )
        surfIntersectionPtr_ = new findCellsIntersectingSurface(mesh_, octree_);

    return surfIntersectionPtr_->facetsIntersectingCells();
}

void removeCellsInSelectedDomains::markSelectedFacets()
{
    const triSurf& surf = octree_.surface();
    const geometricSurfacePatchList& patches = surf.patches();

    VRWGraph usedPatch(patches.size(), 0);

    selectedFacets_.setSize(surf.size());
    forAll(selectedFacets_, rowI)
        selectedFacets_.setRowSize(rowI, 0);

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
                    selectedFacets_[facetsInSubset[i]].append(domainI);
            }
            else
            {
                //- check if it is a patch
                forAll(patches, i)
                {
                    if( patches[i].name() == sName )
                    {
                        usedPatch.append(i, domainI);
                    }
                }
            }
        }
    }

    //- selected facets in used patches
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];

        forAllRow(usedPatch, tri.region(), i)
            selectedFacets_[triI].append(usedPatch(tri.region(), i));
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
            if( selectedFacets_[trianglesInLeaf[i]].size() )
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
    LongList<DynList<label, 2> > intersectedCells(mesh_.cells().size());

    const VRWGraph& facetsIntersectingCell = cellsIntersectedBySurfaceFacets();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(facetsIntersectingCell, cellI)
    {
        forAllRow(facetsIntersectingCell, cellI, triI)
        {
            const label fI = facetsIntersectingCell(cellI, triI);

            forAllRow(selectedFacets_, fI, domainI)
                intersectedCells[cellI].append(selectedFacets_(fI, domainI));
        }
    }

    # ifdef DEBUGRemoveDomains
    labelList domainIds(domains_.size());
    forAll(domainIds, i)
    {
        domainIds[i] =
            mesh_.addCellSubset("intersectingDomain_"+help::scalarToText(i));
    }

    forAll(intersectedCells, i)
    {
        const DynList<label, 2>& doms = intersectedCells[i];

        forAll(doms, j)
            mesh_.addCellToSubset(domainIds[doms[j]], i);
    };
    # endif

    //- TODO: implement the group marking algorithm properly using templates
    //- find islands of cells which are not intersected by the selected domains
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    labelListPMG findCellGroups(intersectedCells.size(), -1);
    label nGroups(0);

    forAll(findCellGroups, cellI)
    {
        if( findCellGroups[cellI] != -1 )
            continue;
        if( intersectedCells[cellI].size() )
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

                if( owner[faceI] == cLabel && neighbour[faceI] >= 0 )
                {
                    neis.append(neighbour[faceI]);
                }
                else if( neighbour[faceI] == cLabel )
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
                if( intersectedCells[nei].size() )
                    continue;

                findCellGroups[nei] = nGroups;
                front.append(nei);
            }
        }

        ++nGroups;
    }

    # ifdef DEBUGRemoveDomains
    Info << "Number of groups " << nGroups << endl;
    Info << "2.1 Here" << endl;
    labelList groupToId(nGroups);
    for(label i=0;i<nGroups;++i)
        groupToId[i] = mesh_.addCellSubset("group_"+help::scalarToText(i));

    forAll(findCellGroups, i)
    {
        if( findCellGroups[i] != -1 )
            mesh_.addCellToSubset(groupToId[findCellGroups[i]], i);
    }
    # endif

    //- collect which domains are assigned to facets neighbouring
    //- groups
    List<DynList<label> > neiDomains(nGroups);
    for(label faceI=0;faceI<mesh_.nInternalFaces();++faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        if( intersectedCells[own].size() && (findCellGroups[nei] != -1) )
        {
            forAll(intersectedCells[own], domainI)
                neiDomains[findCellGroups[nei]].appendIfNotIn
                (
                    intersectedCells[own][domainI]
                );
        }
        else if( intersectedCells[nei].size() && (findCellGroups[own] != -1) )
        {
            forAll(intersectedCells[nei], domainI)
                neiDomains[findCellGroups[own]].appendIfNotIn
                (
                    intersectedCells[nei][domainI]
                );
        }
    }

    # ifdef DEBUGRemoveDomains
    Info << "1. Nei domains " << neiDomains << endl;
    # endif

    //- find a common group for all cell neighbours of a group
    for(label faceI=0;faceI<mesh_.nInternalFaces();++faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        if( intersectedCells[own].size() && (findCellGroups[nei] != -1) )
        {
            const label groupI = findCellGroups[nei];
            forAllReverse(neiDomains[groupI], domainI)
            {
                const label neiDomain = neiDomains[groupI][domainI];
                if( !intersectedCells[own].contains(neiDomain) )
                    neiDomains[groupI].removeElement(domainI);
            }
        }
        else if( intersectedCells[nei].size() && (findCellGroups[own] != -1) )
        {
            const label groupI = findCellGroups[own];
            forAllReverse(neiDomains[groupI], domainI)
            {
                const label neiDomain = neiDomains[groupI][domainI];
                if( !intersectedCells[nei].contains(neiDomain) )
                    neiDomains[groupI].removeElement(domainI);
            }
        }
    }

    # ifdef DEBUGRemoveDomains
    Info << "2. Nei domains " << neiDomains << endl;
    # endif

    //- check which islands of cells correspond to octree boxes
    //- marked as internal boxes
    boolList internalCells(intersectedCells.size(), false);

    forAll(findCellGroups, cellI)
    {
        const label groupI = findCellGroups[cellI];

        //- do not remove cells which are not assigned to a domain
        if( groupI < 0 )
            continue;

        //- do not remove cells intersected by the selected surface elements
        if( intersectedCells[cellI].size() )
            continue;

        //- remove domains surrounded by facets in a single user-selected domain
        if( neiDomains[groupI].size() != 1 )
            continue;

        //- mark cell for removal
        internalCells[cellI] = true;
    }

    //- remove cells inside the selected domains
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

    findAndRemoveCells();

    Info << "Finished removing cells in user-selected domains" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
