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

#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"

//#define DEBUGMorph

# ifdef DEBUGMorph
#include <sstream>
#include "writeMeshEnsight.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfacePartitioner::triSurfacePartitioner
(
	const triSurface& surface
)
:
	surface_(surface),
    corners_(),
	cornerPatches_(),
	partitionPartitions_(surface.patches().size()),
    edgePartitions_(),
    edgePartitionEdgePartitions_(),
    partitionsEdgeParts_(),
    edgePartitionsCorners_()
{
	calculatePartitionAddressing();
}

triSurfacePartitioner::~triSurfacePartitioner()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const labelList& triSurfacePartitioner::corners() const
{
	return corners_;
}

const List<DynList<label> >& triSurfacePartitioner::cornerPatches() const
{
	return cornerPatches_;
}

const List<labelHashSet>& triSurfacePartitioner::partitionPartitions() const
{
	return partitionPartitions_;
}

const labelList& triSurfacePartitioner::edgePartitions() const
{
    return edgePartitions_;
}

const List<labelHashSet>&
triSurfacePartitioner::edgePartitionEdgePartitions() const
{
    return edgePartitionEdgePartitions_;
}

void triSurfacePartitioner::edgePartitionsBetweenPartitions
(
    const label partition1,
    const label partition2,
    DynList<label>& edgePartitions
) const
{
    edgePartitions.clear();
    
    std::pair<label, label> pp
    (
        Foam::min(partition1, partition2),
        Foam::max(partition1, partition2)
    );
    
    std::map<std::pair<label, label>, labelHashSet>::const_iterator it =
        partitionsEdgeParts_.find(pp);
    
    if( it != partitionsEdgeParts_.end() )
    {
        const labelList edgeParts = it->second.toc();
        
        forAll(edgeParts, partI)
            edgePartitions.append(edgeParts[partI]);
    }
}

void triSurfacePartitioner::cornersBetweenEdgePartitions
(
    const label edgePartition1,
    const label edgePartition2,
    DynList<label>& corners
) const
{
    corners.clear();
    
    std::pair<label, label> ep
    (
        Foam::min(edgePartition1, edgePartition2),
        Foam::max(edgePartition1, edgePartition2)
    );
    
    std::map<std::pair<label, label>, labelHashSet>::const_iterator it =
        edgePartitionsCorners_.find(ep);
    
    if( it != edgePartitionsCorners_.end() )
    {
        const labelList corn = it->second.toc();
        
        forAll(corn, cornerI)
            corners.append(corn[cornerI]);
    }
}
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
