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
#include "labelLongList.H"
#include "boolList.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfacePartitioner::calculatePartitionAddressing()
{
    calculateCornersAndAddressing();

    calculatePartitionPartitions();

    calculateEdgePartitions();

    calculatePartitionsToEdgePartitions();

    calculateEdgePartitionsToCorners();
}

void triSurfacePartitioner::calculateCornersAndAddressing()
{
    const VRWGraph& pointFaces = surface_.pointFacets();
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    //- find the number of feature edges connected to each surface node
    List<direction> nEdgesAtNode(surface_.points().size(), direction(0));
    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 )
            continue;

        const label sPatch = surface_[edgeFaces(eI, 0)].region();
        const label ePatch = surface_[edgeFaces(eI, 1)].region();

        if( sPatch != ePatch )
        {
            const edge& e = edges[eI];
            ++nEdgesAtNode[e.start()];
            ++nEdgesAtNode[e.end()];
        }
    }

    //- count the number of feature edges connected to each surface point
    //- corners must have 3 or more edges attached to them
    label nCorners(0);
    forAll(nEdgesAtNode, pI)
    {
        if( nEdgesAtNode[pI] < direction(3) )
            continue;

        ++nCorners;
    }

    corners_.setSize(nCorners);
    cornerPatches_.setSize(nCorners);
    nCorners = 0;

    //- store corner data
    DynList<label> patches;
    forAll(pointFaces, pointI)
    {
        if( nEdgesAtNode[pointI] < direction(3) )
            continue;

        patches.clear();
        forAllRow(pointFaces, pointI, pfI)
            patches.appendIfNotIn(surface_[pointFaces(pointI, pfI)].region());

        corners_[nCorners] = pointI;
        cornerPatches_[nCorners] = patches;
        ++nCorners;
    }
}

void triSurfacePartitioner::calculatePartitionPartitions()
{
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 )
        {
            Warning << "Surface is not a manifold!!" << endl;
            continue;
        }

        const label sPatch = surface_[edgeFaces(eI, 0)].region();
        const label ePatch = surface_[edgeFaces(eI, 1)].region();

        if( sPatch != ePatch )
        {
            partitionPartitions_[sPatch].insert(ePatch);
            partitionPartitions_[ePatch].insert(sPatch);
        }
    }
}

void triSurfacePartitioner::calculateEdgePartitions()
{
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& pointEdges = surface_.pointEdges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    //- make all feature edges
    boolList featureEdge(edgeFaces.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(edgeFaces, eI)
    {
        DynList<label> parts;
        forAllRow(edgeFaces, eI, efI)
            parts.appendIfNotIn(surface_[edgeFaces(eI, efI)].region());

        if( parts.size() > 1 )
            featureEdge[eI] = true;
    }

    //- create a set containing corners for fast searching
    labelHashSet corners;
    forAll(corners_, i)
        corners.insert(corners_[i]);

    edgePartitions_.setSize(edgeFaces. size());
    edgePartitions_ = -1;

    label nPartitions(0);
    forAll(featureEdge, eI)
    {
        if( !featureEdge[eI] )
            continue;

        labelLongList front;
        front.append(eI);
        edgePartitions_[eI] = nPartitions;
        featureEdge[eI] = false;

        while( front.size() )
        {
            const label eLabel = front.removeLastElement();
            const edge& e = edges[eLabel];

            for(label pI=0;pI<2;++pI)
            {
                const label pointI = e[pI];

                if( corners.found(pointI) )
                    continue;

                forAllRow(pointEdges, pointI, peI)
                {
                    const label eJ = pointEdges(pointI, peI);

                    if( featureEdge[eJ] )
                    {
                        edgePartitions_[eJ] = nPartitions;
                        featureEdge[eJ] = false;
                        front.append(eJ);
                    }
                }
            }
        }

        ++nPartitions;
    }

    Info << nPartitions << " edge partitions found!" << endl;

    edgePartitionEdgePartitions_.clear();
    edgePartitionEdgePartitions_.setSize(nPartitions);
}

void triSurfacePartitioner::calculatePartitionsToEdgePartitions()
{
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    forAll(edgeFaces, eI)
    {
        if( edgePartitions_[eI] < 0 )
            continue;

        DynList<label> partitions;
        forAllRow(edgeFaces, eI, efI)
            partitions.appendIfNotIn(surface_[edgeFaces(eI, efI)].region());

        forAll(partitions, i)
        {
            const label partI = partitions[i];
            for(label j=i+1;j<partitions.size();++j)
            {
                const label partJ = partitions[j];

                const std::pair<label, label> pp
                (
                    Foam::min(partI, partJ),
                    Foam::max(partI, partJ)
                );

                partitionsEdgeParts_[pp].insert(edgePartitions_[eI]);
            }
        }
    }
}

void triSurfacePartitioner::calculateEdgePartitionsToCorners()
{
    const VRWGraph& pointEdges = surface_.pointEdges();

    forAll(corners_, cornerI)
    {
        DynList<label> edgePartitionsAtCorner;
        const label pointI = corners_[cornerI];

        forAllRow(pointEdges, pointI, peI)
            edgePartitionsAtCorner.appendIfNotIn
            (
                edgePartitions_[pointEdges(pointI, peI)]
            );

        forAll(edgePartitionsAtCorner, i)
        {
            const label epI = edgePartitionsAtCorner[i];
            if( epI < 0 )
                continue;
            for(label j=i+1;j<edgePartitionsAtCorner.size();++j)
            {
                const label epJ = edgePartitionsAtCorner[j];
                if( epJ < 0 )
                    continue;

                std::pair<label, label> ep
                (
                    Foam::min(epI, epJ),
                    Foam::max(epI, epJ)
                );

                //- create edgepartition - edge partitions addressing
                edgePartitionEdgePartitions_[ep.first].insert(ep.second);
                edgePartitionEdgePartitions_[ep.second].insert(ep.first);

                edgePartitionsCorners_[ep].insert(cornerI);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
