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

#include "triSurfaceDetectMaterials.H"
#include "meshOctree.H"
#include "meshOctreeModifier.H"
#include "meshOctreeInsideOutside.H"
#include "demandDrivenData.H"

//#define DEBUGMaterials

# ifdef DEBUGMaterials
#include "writeOctreeEnsight.H"
#include "helperFunctions.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectMaterials::createPartitions()
{
    facePartition_.setSize(surf_.size());
    facePartition_ = -1;
    nPartitions_ = 0;

    const labelListList& faceEdges = surf_.faceEdges();
    const labelListList& edgeFaces = surf_.edgeFaces();

    forAll(surf_, triI)
    {
        if( facePartition_[triI] != -1 )
            continue;

        facePartition_[triI] = nPartitions_;
        labelListPMG front;
        front.append(triI);

        while( front.size() )
        {
            const label fLabel = front.removeLastElement();

            const labelList& fEdges = faceEdges[fLabel];
            forAll(fEdges, feI)
            {
                const label edgeI = fEdges[feI];

                if( edgeFaces[edgeI].size() != 2 )
                    continue;

                label nei = edgeFaces[edgeI][0];
                if( nei == fLabel )
                    nei = edgeFaces[edgeI][1];

                if( facePartition_[nei] == -1 )
                {
                    front.append(nei);
                    facePartition_[nei] = nPartitions_;
                }
            }
        }

        ++nPartitions_;
    }

    Info << "Found " << nPartitions_ << " surface partitions" << endl;
    # ifdef DEBUGMaterials
    for(label i=0;i<nPartitions_;++i)
    {
        labelList pointMap, faceMap;
        boolList keepTriangle(surf_.size(), false);
        forAll(facePartition_, triI)
            if( facePartition_[triI] == i )
                keepTriangle[triI] = true;

        triSurface subset = surf_.subsetMesh(keepTriangle, pointMap, faceMap);
        fileName fName = "partition";
        fName += help::scalarToText(i);
        fName += ".stl";
        subset.write(fName);
    }
    # endif
}

void triSurfaceDetectMaterials::createOctree()
{
    if( !octreePtr_ )
        octreePtr_ = new meshOctree(surf_);

    const edgeList& edges = surf_.edges();

    if( edges.size() == 0 )
        FatalError << "Surface mesh has no edges!!" << exit(FatalError);

    const pointField& points = surf_.localPoints();

    scalar averageLength(0.0), minLength(VGREAT);
    forAll(edges, eI)
    {
        const scalar l = edges[eI].mag(points);
        minLength = Foam::min(l, minLength);
        averageLength += l;
    }

    averageLength /= edges.size();

    averageLength = (averageLength + minLength) / 2.0;

    const boundBox& bb = octreePtr_->rootBox();
    minLength = bb.max().x() - bb.min().x();

    direction maxLevel(0);
    while( minLength > averageLength )
    {
        ++maxLevel;
        minLength /= 2.0;
    }

    maxLevel = Foam::min(maxLevel, direction(10));

    Info << "Refining octree to level " << label(maxLevel) << endl;

    meshOctreeModifier octreeModifier(*octreePtr_);

    octreeModifier.createListOfLeaves();

    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    //- refine boxes containing triangles to the given refinement level
    label nMarked;
    do
    {
        nMarked = 0;

        List<direction> refineBox(leaves.size(), direction(0));
        forAll(leaves, leafI)
        {
            if( !leaves[leafI]->hasContainedElements() )
                continue;
            if( leaves[leafI]->level() >= maxLevel )
                continue;

            ++nMarked;
            refineBox[leafI] = 1;
        }

        octreeModifier.refineSelectedBoxes(refineBox);
    } while( nMarked && (leaves.size() < 500000) );

    Info << "Number of leaves after initial refinement "
        << leaves.size() << endl;

    //- refine coarse boxes. This is performed by refinining octree boxes
    //- containing more than one partition in its vicinity
    const boundBox& rootBox = octreePtr_->rootBox();
    const labelListList& pointTriangles = surf_.pointFaces();

    DynList<const meshOctreeCube*, 256> neis;
    do
    {
        nMarked = 0;

        List<direction> refineBox(leaves.size(), direction(0));
        forAll(leaves, leafI)
        {
            if( !leaves[leafI]->hasContainedElements() )
                continue;

            const point c = leaves[leafI]->centre(rootBox);
            const scalar s = leaves[leafI]->size(rootBox);

            const boundBox bb
            (
                c - point(s, s, s),
                c + point(s, s, s)
            );

            //- find leaves within the given box
            octreeModifier.findLeavesContainedInBox(bb, neis);

            //- find triangles and points within the boundBox
            labelHashSet containedElmts, containedPts;
            bool hasCorner(false), hasEdge(false);
            forAll(neis, i)
            {
                if( !neis[i]->hasContainedElements() )
                    continue;

                const label ceI = neis[i]->containedElements();
                const VRWGraph& containedTriangles =
                    neis[i]->slotPtr()->containedTriangles_;
                forAllRow(containedTriangles, ceI, j)
                {
                    containedElmts.insert(containedTriangles(ceI, j));

                    const labelledTri& tri = surf_[containedTriangles(ceI, j)];
                    forAll(tri, k)
                        containedPts.insert(tri[k]);
                }
            }

            //- check if there exist any corners or edges in the boundBox
            forAllConstIter(labelHashSet, containedPts, pIter)
            {
                const labelList& pTriangles = pointTriangles[pIter.key()];

                labelHashSet partitions;
                forAll(pTriangles, i)
                    partitions.insert(facePartition_[pTriangles[i]]);

                if( partitions.size() > 2 )
                {
                    hasCorner = true;
                }
                else if( partitions.size() > 1 )
                {
                    hasEdge = true;
                }
            }

            if( hasCorner )
                continue;

            //- find partitions contained in the boundBox
            labelHashSet containedPartitions;
            forAllConstIter(labelHashSet, containedElmts, tIter)
                containedPartitions.insert(facePartition_[tIter.key()]);

            if(
                (hasEdge && (containedPartitions.size() > 2)) ||
                (!hasEdge && (containedPartitions.size() > 1))
            )
            {
                forAll(neis, i)
                    refineBox[neis[i]->cubeLabel()] = 1;
                ++nMarked;
            }
        }

        Info << "Marked " << nMarked << " octree cubes "
            << "for local refinement" << endl;

        if( nMarked )
            octreeModifier.refineSelectedBoxes(refineBox);
    } while( nMarked );

    //- create inside/outside information
    meshOctreeInsideOutside octreeIO(*octreePtr_);

    label nInside(0), nOutside(0), nUnknown(0);
    for(label i=0;i<octreePtr_->numberOfLeaves();++i)
    {
        switch( octreePtr_->returnLeaf(i).cubeType() )
        {
            case meshOctreeCubeBasic::INSIDE:
            {
                ++nInside;
            } break;
            case meshOctreeCubeBasic::OUTSIDE:
            {
                ++nOutside;
            } break;
            case meshOctreeCubeBasic::UNKNOWN:
            {
                ++nUnknown;
            } break;
        };
    }
    Info << "Number of INSIDE leaves " << nInside << endl;
    Info << "Number of OUTSIDE leaves " << nOutside << endl;
    Info << "Number of UNKNOWN leaves " << nUnknown << endl;

    # ifdef DEBUGMaterials
    writeOctreeEnsight(*octreePtr_, "refinedOctree", meshOctreeCubeBasic::DATA);
    # endif
}

void triSurfaceDetectMaterials::refineOctree()
{
    meshOctreeModifier octreeModifier(*octreePtr_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    List<direction> refineBox(leaves.size(), direction(0));

    forAll(partitionMaterials_, partI)
    {
        if(
            (partitionMaterials_.sizeOfRow(partI) == 0) ||
            (partitionMaterials_.sizeOfRow(partI) > 2) ||
            (
                (partitionType_[partI] & meshOctreeCubeBasic::OUTSIDE) &&
                (partitionType_[partI] & meshOctreeCubeBasic::UNKNOWN)
            )
        )
        {
            forAll(leaves, leafI)
            {
                if( !leaves[leafI]->hasContainedElements() )
                    continue;

                const label ceI = leaves[leafI]->containedElements();
                const VRWGraph& containedTriangles =
                    leaves[leafI]->slotPtr()->containedTriangles_;
                forAllRow(containedTriangles, ceI, i)
                {
                    if( facePartition_[containedTriangles(ceI, i)] == partI )
                    {
                        refineBox[leafI] = 1;
                        break;
                    }
                }
            }
        }
    }

    octreeModifier.refineSelectedBoxes(refineBox);
}

void triSurfaceDetectMaterials::findOctreeGroups()
{
    const label nLeaves = octreePtr_->numberOfLeaves();

    octreeGroupForBox_.setSize(nLeaves);
    octreeGroupForBox_ = -1;
    nOctreeGroups_ = 0;

    DynList<label> neighs(24);
    for(label leafI=0;leafI<nLeaves;++leafI)
    {
        if( octreeGroupForBox_[leafI] != -1 )
            continue;
        if( octreePtr_->returnLeaf(leafI).cubeType() & meshOctreeCube::OUTSIDE )
            continue;
        if( octreePtr_->hasContainedTriangles(leafI) )
            continue;

        labelListPMG front;
        octreeGroupForBox_[leafI] = nOctreeGroups_;
        front.append(leafI);

        while( front.size() )
        {
            const label fLabel = front.removeLastElement();
            octreePtr_->findNeighboursForLeaf(fLabel, neighs);

            forAll(neighs, i)
            {
                const label nei = neighs[i];

                if( nei < 0 )
                    continue;
                if( octreeGroupForBox_[nei] != -1 )
                    continue;
                if( octreePtr_->hasContainedTriangles(nei) )
                    continue;

                octreeGroupForBox_[nei] = nOctreeGroups_;
                front.append(nei);
            }
        }

        ++nOctreeGroups_;
    }

    Info << "Found " << nOctreeGroups_ << " groups of octree boxes" << endl;
}

void triSurfaceDetectMaterials::findMaterialsAndWalls()
{
    if( partitionMaterials_.size() == 0 )
    {
        partitionMaterials_.setSizeAndColumnWidth(nPartitions_, nOctreeGroups_);
    }
    else
    {
        for(label partI=0;partI<nPartitions_;++partI)
            partitionMaterials_.setRowSize(partI, 0);
    }
    partitionType_.setSize(nPartitions_);
    partitionType_ = 0;

    const label nLeaves = octreePtr_->numberOfLeaves();
    DynList<label> containedTriangles(64), neighbours(24);
    for(label leafI=0;leafI<nLeaves;++leafI)
    {
        //const meshOctreeCubeBasic& oc = octreePtr_->returnLeaf(leafI);
        if( !octreePtr_->hasContainedTriangles(leafI) )
            continue;

        //- find triangles contained in this cube
        octreePtr_->containedTriangles(leafI, containedTriangles);

        //- find face partitions this cube belongs to
        labelHashSet containedPartitions;
        forAll(containedTriangles, i)
            containedPartitions.insert(facePartition_[containedTriangles[i]]);

        if( containedPartitions.size() > 1 )
            continue;

        const label partI = containedPartitions.toc()[0];
        containedPartitions.clear();

        octreePtr_->findNeighboursForLeaf(leafI, neighbours);

        forAll(neighbours, i)
        {
            const label nei = neighbours[i];
            if( nei < 0 )
                continue;
            if( octreePtr_->hasContainedTriangles(nei) )
                continue;

            const meshOctreeCubeBasic& neiOc = octreePtr_->returnLeaf(nei);

            if( neiOc.cubeType() & meshOctreeCubeBasic::OUTSIDE )
            {
                partitionType_[partI] |= meshOctreeCubeBasic::OUTSIDE;
            }
            else if( neiOc.cubeType() & meshOctreeCubeBasic::INSIDE )
            {
                partitionMaterials_.appendIfNotIn
                (
                    partI,
                    octreeGroupForBox_[nei]
                );

                partitionType_[partI] |= meshOctreeCubeBasic::INSIDE;
            }
            else if( neiOc.cubeType() & meshOctreeCubeBasic::UNKNOWN )
            {
                partitionMaterials_.appendIfNotIn
                (
                    partI,
                    octreeGroupForBox_[nei]
                );

                partitionType_[partI] |= meshOctreeCubeBasic::UNKNOWN;
            }
        }
    }

    //- there shall be no partitions with more than 2 octree groups
    //- partitions having some OUTSIDE neighbours shall have only one group
    labelList newGroupLabels(nOctreeGroups_, -1);

    //- group type is necessary for sampling groups
    List<direction> groupType(nOctreeGroups_, direction(0));
    for(label leafI=0;leafI<nLeaves;++leafI)
    {
        const label groupI = octreeGroupForBox_[leafI];
        if( groupI == -1 )
            continue;

        groupType[groupI] |= octreePtr_->returnLeaf(leafI).cubeType();
    }

    //- sort outer boundaries
    nOctreeGroups_ = 0;
    forAll(partitionMaterials_, partitionI)
    {
        if( !(partitionType_[partitionI] & meshOctreeCubeBasic::OUTSIDE) )
            continue;
        if( partitionType_[partitionI] & meshOctreeCubeBasic::UNKNOWN )
            continue;
        if( partitionMaterials_.sizeOfRow(partitionI) == 0 )
            continue;

        //- this partition shall be in one octree group
        label newGroupLabel(-1);
        forAllRow(partitionMaterials_, partitionI, i)
        {
            const label groupI = partitionMaterials_(partitionI, i);

            if( groupType[groupI] != meshOctreeCubeBasic::INSIDE )
                continue;
            if( newGroupLabels[groupI] == -1 )
                continue;

            newGroupLabel = newGroupLabels[groupI];
        }

        if( newGroupLabel == -1 )
            newGroupLabel = nOctreeGroups_++;

        forAllRow(partitionMaterials_, partitionI, i)
        {
            const label groupI = partitionMaterials_(partitionI, i);
            if( groupType[groupI] == meshOctreeCubeBasic::INSIDE )
                newGroupLabels[groupI] = newGroupLabel;
        }
    }

    //- internal partitions contain one or two octree groups
    //- zero-thickness walls contain one partition
    //- internal surface contains two partitions
    forAll(partitionMaterials_, partitionI)
    {
        if( partitionType_[partitionI] & meshOctreeCubeBasic::OUTSIDE )
            continue;
        if( !(partitionType_[partitionI] & meshOctreeCubeBasic::UNKNOWN) )
            continue;
        if( partitionMaterials_.sizeOfRow(partitionI) == 0 )
            continue;

        //- this partition shall be in one octree group
        label newGroupLabel(-1);
        forAllRow(partitionMaterials_, partitionI, i)
        {
            const label groupI = partitionMaterials_(partitionI, i);

            if( groupType[groupI] != meshOctreeCubeBasic::UNKNOWN )
                continue;
            if( newGroupLabels[groupI] == -1 )
                continue;

            newGroupLabel = newGroupLabels[groupI];
        }

        if( newGroupLabel == -1 )
            newGroupLabel = nOctreeGroups_++;

        forAllRow(partitionMaterials_, partitionI, i)
        {
            const label groupI = partitionMaterials_(partitionI, i);
            if( groupType[groupI] == meshOctreeCubeBasic::UNKNOWN )
                newGroupLabels[groupI] = newGroupLabel;
        }
    }

    //- rename octree groups
    forAll(octreeGroupForBox_, cubeI)
    {
        if( octreeGroupForBox_[cubeI] == -1 )
            continue;

        octreeGroupForBox_[cubeI] = newGroupLabels[octreeGroupForBox_[cubeI]];
    }

    //- renumber partitionMaterials_ according to newGroupLabels
    forAll(partitionMaterials_, partitionI)
    {
        DynList<label> partitions(partitionMaterials_.sizeOfRow(partitionI));
        forAllRow(partitionMaterials_, partitionI, i)
            partitions.appendIfNotIn
            (
                newGroupLabels[partitionMaterials_(partitionI, i)]
            );

        partitionMaterials_.setRow(partitionI, partitions);
    }

    Info << "Finished finding partitions and materials" << endl;
}

bool triSurfaceDetectMaterials::checkMaterials() const
{
    forAll(partitionMaterials_, partitionI)
    {
        if( partitionMaterials_.sizeOfRow(partitionI) == 0 )
            return true;
        if( partitionMaterials_.sizeOfRow(partitionI) > 2 )
            return true;

        if(
            (partitionType_[partitionI] & meshOctreeCubeBasic::OUTSIDE) &&
            (partitionType_[partitionI] & meshOctreeCubeBasic::UNKNOWN)
        )
            return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
