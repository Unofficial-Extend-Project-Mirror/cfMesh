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
    meshSurfacePartitioner

Description
    Finds corners and edge nodes at the surface of the volume mesh

SourceFiles
    meshSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfacePartitioner_H
#define meshSurfacePartitioner_H

#include "meshSurfaceEngine.H"
#include "HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class meshSurfacePartitioner
{
    // Private data
        //- reference to mesh surface engine
        const meshSurfaceEngine& meshSurface_;

        //- pointer to face patches
        const labelList& facePatch_;

        //- labels of corner nodes
        labelHashSet corners_;

        //- patches at a corner node
        std::map<label, DynList<label> > cornerPatches_;

        //- labels of edge points
        labelHashSet edgePoints_;

        //- patch-patches addressing
        List<labelHashSet> patchPatches_;

        //- number of edges attached to a surface point
        labelList nEdgesAtPoint_;

        //- labels of feature edges
        labelHashSet featureEdges_;

    // Private member functions
        //- find corners, edge nodes and addressing
        void calculateCornersEdgesAndAddressing();

        //- Disallow default bitwise copy construct
        meshSurfacePartitioner(const meshSurfacePartitioner&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfacePartitioner&);

public:

    // Constructors

        //- Construct from meshSurfaceEngine
        meshSurfacePartitioner(const meshSurfaceEngine&);

        //- Construct from meshSurfaceEngine and face patches
        meshSurfacePartitioner
        (
            const meshSurfaceEngine& meshSurface,
            const labelList& facePatch
        );

    // Destructor

        ~meshSurfacePartitioner();

    // Member Functions
        //- return labels of corner points (from the list of boundary points)
        inline const labelHashSet& corners() const
        {
            return corners_;
        }

        //- return label of patches near corner points
        inline const std::map<label, DynList<label> >& cornerPatches() const
        {
            return cornerPatches_;
        }

        //- return labels of edge points (from the list of boundary points)
        inline const labelHashSet& edgePoints() const
        {
            return edgePoints_;
        }

        //- return the number of feature edges attached to a boundary point
        inline label numberOfFeatureEdgesAtPoint(const label bpI) const
        {
            return nEdgesAtPoint_[bpI];
        }

        //- return patch-patches addressing
        inline const List<labelHashSet>& patchPatches() const
        {
            return patchPatches_;
        }

        //- return labels of boundary edges which are feature edges
        inline const labelHashSet& featureEdges() const
        {
            return featureEdges_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
