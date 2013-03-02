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
    
        //- labels of corner nodes
        labelHashSet corners_;
    
        //- labels of edge nodes
        labelHashSet edgeNodes_;
    
        //- partition-partitions addressing
        List<labelHashSet> partitionPartitions_;
    
        //- number of edges attached to a surface point
        labelList nEdgesAtPoint_;

    // Private member functions
        //- find corners, edge nodes and addressing
        void calculateCornersEdgesAndAddressing();

        //- Disallow default bitwise copy construct
        meshSurfacePartitioner(const meshSurfacePartitioner&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfacePartitioner&);

public:
    
    // Constructors

        //- Construct from mesh
        meshSurfacePartitioner(const meshSurfaceEngine&);

    // Destructor

        ~meshSurfacePartitioner();

    // Member Functions
        //- return labels of corner nodes (from the list of boundary points)
        const labelHashSet& corners() const;
    
        //- return labels of edge nodes (from the list of boundary points)
        const labelHashSet& edgeNodes() const;

        //- return the number of feature edges attached to a boundary point
        label numberOfFeatureEdgesAtPoint(const label bpI) const;

        //- return partition-partitions addressing
        const List<labelHashSet>& partitionPartitions() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
