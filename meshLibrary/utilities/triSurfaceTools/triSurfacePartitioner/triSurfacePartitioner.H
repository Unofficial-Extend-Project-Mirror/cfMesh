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
    triSurfacePartitioner

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfacePartitioner_H
#define triSurfacePartitioner_H

#include "triSurface.H"
#include "DynList.H"
#include "HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class triSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class triSurfacePartitioner
{
    // Private data
        //- reference to triSurface
        const triSurface& surface_;
    
        //- corner nodes
        labelList corners_;
        List<DynList<label> > cornerPatches_;
    
        //- information which partitions share an edge with a given partition
        List<labelHashSet> partitionPartitions_;
    
        //- edge partitions
        labelList edgePartitions_;
    
        //- information which edge partitions share a corner
        //- with a given partition
        List<labelHashSet> edgePartitionEdgePartitions_;
    
        //- partition of the edge between the surface partitions
        std::map<std::pair<label, label>, labelHashSet> partitionsEdgeParts_;
    
        //- corners shared by edge partitions
        std::map<std::pair<label, label>, labelHashSet> edgePartitionsCorners_;

    // Private member functions
        //- calculate partitions
        void calculatePartitionAddressing();
    
        //- find surface corners
        void calculateCornersAndAddressing();
    
        //- calculate partitionPartitions_ addressing
        void calculatePartitionPartitions();
        
        //- calculate edge partitions
        void calculateEdgePartitions();
        
        //- calculate surface partition to edge partitions addressing
        void calculatePartitionsToEdgePartitions();
        
        //- calculate edge partitions to corner addressing
        void calculateEdgePartitionsToCorners();

        //- Disallow default bitwise copy construct
        triSurfacePartitioner(const triSurfacePartitioner&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfacePartitioner&);

public:

    // Constructors

        //- Construct from triSurface
        triSurfacePartitioner(const triSurface& surface);

    // Destructor

        ~triSurfacePartitioner();


    // Member Functions
        //- return corner nodes
        const labelList& corners() const;
        
        //- return corner patches
        const List<DynList<label> >& cornerPatches() const;
        
        //- return partition-partitions addressing
        const List<labelHashSet>& partitionPartitions() const;

        //- return edge partitions. Partitions of edges which are not
        //- the feature edges is set to -1
        const labelList& edgePartitions() const;
        
        //- Edge partition - edge partitions addressing
        const List<labelHashSet>& edgePartitionEdgePartitions() const;
        
        //- Return partitions of feature edges between the given partitions
        void edgePartitionsBetweenPartitions
        (
            const label partition1,
            const label partition2,
            DynList<label>& edgePartitions
        ) const;
        
        //- return corners shared between the given edge partitions
        void cornersBetweenEdgePartitions
        (
            const label edgePartition1,
            const label edgePartition2,
            DynList<label>& corners
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
