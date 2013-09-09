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

Class
    removeCellsInSelectedDomains

Description
    Removes cells inside the domains in the multi-domains mesh
    which are specified by the user

SourceFiles
    removeCellsInSelectedDomains.C

\*---------------------------------------------------------------------------*/

#ifndef removeCellsInSelectedDomains_H
#define removeCellsInSelectedDomains_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class findCellsIntersectingSurface;

/*---------------------------------------------------------------------------*\
                Class removeCellsInSelectedDomains Declaration
\*---------------------------------------------------------------------------*/

class removeCellsInSelectedDomains
{
    // Private data
        //- reference to mesh
        polyMeshGen& mesh_;

        //- const reference to meshOctree
        const meshOctree& octree_;

        //- regions which shall not be meshed
        DynList<DynList<word, 4> > domains_;

        //- selected facets
        boolList selectedFacets_;

        //- intersect cells with the surface mesh
        findCellsIntersectingSurface* surfIntersectionPtr_;

        //- octree leaves inside the domains which shall be removed
        boolList insideLeaves_;

    // Private member functions
        //- return a graph containing information which surface facets
        //- intersect each cell in the mesh
        const VRWGraph& cellsIntersectedBySurfaceFacets();

        //- find facets contained in user-selected patches/subsets
        void markSelectedFacets();

        //- find octree boxed inside the given regions
        void findLeavesInsideRegions();

        //- perform removal of cells
        void findAndRemoveCells();

        //- Disallow default bitwise copy construct
        removeCellsInSelectedDomains(const removeCellsInSelectedDomains&);

        //- Disallow default bitwise assignment
        void operator=(const removeCellsInSelectedDomains&);

public:

    // Constructors

        //- Construct from mesh and octree
        removeCellsInSelectedDomains
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

    // Destructor

        ~removeCellsInSelectedDomains();

    // Public member functions
        //- assing a domain by giving a list of patches/subsets bounding it
        //- shall be called prior to removeCells
        void selectCellsInDomain(const wordList&);

        //- perform removal of cells
        void removeCells();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
