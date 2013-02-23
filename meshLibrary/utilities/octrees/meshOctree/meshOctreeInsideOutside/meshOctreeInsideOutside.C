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

#include "meshOctreeInsideOutside.H"
#include "triSurf.H"
#include "boundBox.H"
#include "labelListPMG.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOctreeInsideOutside::meshOctreeInsideOutside
(
    meshOctree& octree
)
:
    octreeModifier_(octree),
    cubeGroup_(octree.numberOfLeaves(), -1),
    cubesInGroup_(),
    groupType_(),
    boundaryDATACubes_(),
    hasOutsideNeighbour_(octree.numberOfLeaves(), false),
    communicationCubes_(),
    neighbouringGroups_()
{
    initialiseBoxes();

    frontalMarking();

    markOutsideCubes();

    reviseDataBoxes();

    markInsideCubes();

    label nInternal(0), nUnknown(0), nData(0), nOutside(0);

    const label nLeaves = octree.numberOfLeaves();
    for(label leafI=0;leafI<nLeaves;++leafI)
    {
        const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);

        if( oc.cubeType() & meshOctreeCube::INSIDE )
        {
            ++nInternal;
        }
        else if( oc.cubeType() & meshOctreeCube::UNKNOWN )
        {
            ++nUnknown;
        }
        else if( oc.cubeType() & meshOctreeCube::DATA )
        {
            ++nData;
        }
        else if( oc.cubeType() & meshOctreeCube::OUTSIDE )
        {
            ++nOutside;
        }
    }

    if( octree.neiProcs().size() )
    {
        reduce(nInternal, sumOp<label>());
        reduce(nUnknown, sumOp<label>());
        reduce(nData, sumOp<label>());
        reduce(nOutside, sumOp<label>());
    }

    Info << "Number of internal boxes is " << nInternal << endl;
    Info << "Number of outside boxes is " << nOutside << endl;
    Info << "Number of data boxes is " << nData << endl;
    Info << "Number of unknown boxes is " << nUnknown << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOctreeInsideOutside::~meshOctreeInsideOutside()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeInsideOutside::initialiseBoxes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();

    # pragma omp parallel for if( leaves.size() > 1000 )
    forAll(leaves, leafI)
    {
        if( leaves[leafI]->hasContainedElements() )
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::DATA);
        }
        else
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::UNKNOWN);
        }
    }
}

void meshOctreeInsideOutside::frontalMarking()
{
    communicationCubes_.clear();
    neighbouringGroups_.clear();

    labelListPMG frontCubes;
    DynList<label> neighbours(24);

    label nGroup(0), chunkI(0), nChunks, chunkSize;

    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();

    boolList commCubes(leaves.size(), false);

    # pragma omp parallel if( leaves.size() > 1000 ) \
    private(frontCubes, neighbours)
    {
        LongList<std::pair<label, label> > threadCommPairs;

        # pragma omp master
        {
            nChunks = 3 * omp_get_num_threads();
            chunkSize = leaves.size() / nChunks + 1;
        }

        # pragma omp barrier

        while( chunkI < nChunks )
        {
            label minLeaf, maxLeaf;
            # pragma omp critical
            minLeaf = chunkI++ * chunkSize;

            if( minLeaf >= leaves.size() )
                break;

            maxLeaf = Foam::min(leaves.size(), minLeaf + chunkSize);

            for(label leafI=minLeaf;leafI<maxLeaf;++leafI)
            {
                if( leaves[leafI]->hasContainedElements() )
                    continue;
                if( cubeGroup_[leafI] != -1 )
                    continue;

                label groupI;
                # pragma omp critical
                groupI = nGroup++;

                direction cType(meshOctreeCubeBasic::UNKNOWN);
                frontCubes.clear();
                frontCubes.append(leafI);
                cubeGroup_[leafI] = groupI;

                labelListPMG neiDATACubes;

                while( frontCubes.size() )
                {
                    const label fLabel = frontCubes.removeLastElement();
                    octree.findNeighboursForLeaf(fLabel, neighbours);

                    forAll(neighbours, neiI)
                    {
                        const label nei = neighbours[neiI];
                        if( (nei >= minLeaf) && (nei < maxLeaf) )
                        {
                            if( cubeGroup_[nei] != -1 )
                                continue;

                            if( leaves[nei]->hasContainedElements() )
                            {
                                neiDATACubes.append(nei);
                            }
                            else
                            {
                                frontCubes.append(nei);
                                cubeGroup_[nei] = groupI;
                            }
                        }
                        else if( nei == -1 )
                        {
                            cType = meshOctreeCubeBasic::OUTSIDE;
                        }
                        else if( nei == meshOctreeCubeBasic::OTHERPROC )
                        {
                            commCubes[fLabel] = true;
                        }
                        else
                        {
                            if( leaves[nei]->hasContainedElements() )
                            {
                                neiDATACubes.append(nei);
                            }
                            else
                            {
                                threadCommPairs.append
                                (
                                    std::make_pair(fLabel, nei)
                                );
                            }
                        }
                    }
                }

                # pragma omp critical
                {
                    if( groupI >= boundaryDATACubes_.size() )
                        boundaryDATACubes_.setSize(groupI+1);

                    boundaryDATACubes_.setRow(groupI, neiDATACubes);
                    groupType_[groupI] = cType;
                }
            }
        }

        # pragma omp barrier

        # pragma omp master
        neighbouringGroups_.setSize(nGroup);

        # pragma omp barrier

        forAll(threadCommPairs, pairI)
        {
            const std::pair<label, label>& cubesPair = threadCommPairs[pairI];
            const label groupI = cubeGroup_[cubesPair.first];
            const label neiGroup = cubeGroup_[cubesPair.second];

            if( neiGroup >= nGroup )
                FatalError << "neiGroup " << neiGroup << " is >= than "
                    << "nGroups " << nGroup << abort(FatalError);

            if(
                (neiGroup != -1) &&
                !neighbouringGroups_.contains(groupI, neiGroup)
            )
            {
                # pragma omp critical
                neighbouringGroups_.append(groupI, neiGroup);
            }
        }
    }

    //- create cubesInGroup_ addressing
    labelList nCubesInGroup(nGroup, 0);
    forAll(cubeGroup_, leafI)
    {
        if( cubeGroup_[leafI] < 0 )
            continue;

        ++nCubesInGroup[cubeGroup_[leafI]];
    }

    cubesInGroup_.setSizeAndRowSize(nCubesInGroup);

    forAllReverse(cubeGroup_, leafI)
    {
        const label groupI = cubeGroup_[leafI];

        if( groupI < 0 )
            continue;

        cubesInGroup_(groupI, --nCubesInGroup[groupI]) = leafI;
    }

    //- mark cubes at inter-processor boundaries
    forAll(commCubes, leafI)
    {
        if( commCubes[leafI] )
            communicationCubes_.append(leafI);
    }

    # ifdef DEBUGSearch
    label nMarked(0);
    forAll(cubeGroup_, leafI)
    {
        if( cubeGroup_[leafI] != -1 )
            ++nMarked;
    }
    reduce(nMarked, sumOp<label>());
    const label totalLeaves = returnReduce(leaves_.size(), sumOp<label>());
    Info << "Total number of leaves " << totalLeaves << endl;
    Info << "Number of marked leaves " << nMarked << endl;
    # endif
}

void meshOctreeInsideOutside::markOutsideCubes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();

    DynList<label> neighbours(24);
    label nChanged;
    bool keepUpdating;

    do
    {
        keepUpdating = false;

        do
        {
            nChanged = 0;

            //- make sure that groups created by different threads
            //- have the same information
            forAll(neighbouringGroups_, groupI)
            {
                if( groupType_[groupI] & meshOctreeCubeBasic::OUTSIDE )
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_(groupI, i);
                        if( groupType_[neiGroup] & meshOctreeCube::UNKNOWN )
                        {
                            ++nChanged;
                            groupType_[neiGroup] = meshOctreeCube::OUTSIDE;
                        }
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        do
        {
            nChanged = 0;
            LongList<meshOctreeCubeCoordinates> dataToSend;

            //- go through the list of communicationCubes and send the ones
            //- which are marked as outside
            forAll(communicationCubes_, i)
            {
                const label groupI = cubeGroup_[communicationCubes_[i]];

                if( groupI < 0 )
                    continue;

                if( groupType_[groupI] & meshOctreeCube::OUTSIDE )
                    dataToSend.append(*leaves[communicationCubes_[i]]);
            }

            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree.exchangeRequestsWithNeighbourProcessors
            (
                dataToSend,
                receivedCoords
            );

            //- go through the list of received coordinates and check if any
            //- local boxes are their neighbours. If a local neighbour is
            //- a DATA box set the hasOutsideNeighbour_ flag to true. If the
            //- local neighbour is of UNKNOWN type set it to OUTSIDE.
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20)
            forAll(receivedCoords, i)
            {
                octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if( leaves[nei]->hasContainedElements() )
                    {
                        hasOutsideNeighbour_[nei] = true;
                        continue;
                    }
                    if( groupType_[cubeGroup_[nei]] & meshOctreeCube::UNKNOWN )
                    {
                        groupType_[cubeGroup_[nei]] = meshOctreeCube::OUTSIDE;
                        ++nChanged;
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

            reduce(nChanged, sumOp<label>());
        } while( nChanged != 0 );

        reduce(keepUpdating, maxOp<bool>());

    } while( keepUpdating );

    //- set OUTSIDE type to the cubes in OUTSIDE groups
    for
    (
        std::map<label, direction>::const_iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        if( it->first < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::OUTSIDE )
        {
            const label groupI = it->first;

            //- set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_(groupI, i)]->setCubeType
                (
                    meshOctreeCube::OUTSIDE
                );

            //- set true to the collected DATA boxes
            forAllRow(boundaryDATACubes_, groupI, neiI)
                hasOutsideNeighbour_[boundaryDATACubes_(groupI, neiI)] = true;
        }
    }
}

void meshOctreeInsideOutside::reviseDataBoxes()
{
    //- remove DATA flag from boxes which do not have an OUTSIDE neighbour
    //- and are not surrounded with DATA boxes containing different surface
    //- triangles in different patches
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();
    const triSurface& surface = octree.surface();
    DynList<label> neighbours(60);

    boolList checkedPatches(leaves.size(), false);

    label nMarked;

    do
    {
        nMarked = 0;

        LongList<meshOctreeCubeCoordinates> checkCoordinates;
        labelHashSet transferCoordinates;

        # pragma omp parallel for if( leaves.size() > 1000 ) \
        private(neighbours) schedule(dynamic, 20) reduction(+ : nMarked)
        forAll(leaves, leafI)
            if( Pstream::parRun() && hasOutsideNeighbour_[leafI] )
            {
                octree.findAllLeafNeighbours(leafI, neighbours);
                forAll(neighbours, neiI)
                    if( neighbours[neiI] == meshOctreeCubeBasic::OTHERPROC )
                    {
                        # pragma omp critical
                        {
                            if( !transferCoordinates.found(leafI) )
                            {
                                checkCoordinates.append
                                (
                                    leaves[leafI]->coordinates()
                                );
                                transferCoordinates.insert(leafI);
                            }
                        }

                        break;
                    }
            }
            else if(
                (leaves[leafI]->cubeType() & meshOctreeCube::DATA) &&
                !hasOutsideNeighbour_[leafI]
            )
            {
                meshOctreeCube* oc = leaves[leafI];

                # ifdef DEBUGSearch
                Info << "Box " << leafI << " may not be a DATA box" << endl;
                # endif

                DynList<label> patches;
                const VRWGraph& ct =
                    oc->slotPtr()->containedTriangles_;
                const constRow el = ct[oc->containedElements()];
                forAll(el, elI)
                    patches.appendIfNotIn(surface[el[elI]].region());

                if( patches.size() > 1 )
                    continue;

                checkedPatches[leafI] = true;

                //- check if there exist neighbours
                //- which have some DATA neighbours
                octree.findAllLeafNeighbours(leafI, neighbours);
                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if( hasOutsideNeighbour_[nei] )
                    {
                        oc->setCubeType(meshOctreeCube::INSIDE);

                        ++nMarked;
                        break;
                    }
                }
            }

        if( octree.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree.exchangeRequestsWithNeighbourProcessors
            (
                checkCoordinates,
                receivedCoords
            );

            //- check if any of the local neighbours is a data box with
            //- no OUTSIDE neighbours
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20) reduction(+ : nMarked)
            forAll(receivedCoords, i)
            {
                octree.findAllLeafNeighbours(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if
                    (
                        (leaves[nei]->cubeType() & meshOctreeCube::DATA) &&
                        !hasOutsideNeighbour_[nei] && checkedPatches[nei]
                    )
                    {
                        leaves[nei]->setCubeType(meshOctreeCube::INSIDE);

                        ++nMarked;
                    }
                }
            }

            reduce(nMarked, sumOp<label>());
        }
    } while( nMarked != 0 );

    # ifdef DEBUGSearch
    label nOutside(0), nData(0), hasOutNei(0);
    forAll(leaves, leafI)
    {
        const direction cType = leaves[leafI]->cubeType();
        if( cType & meshOctreeCubeBasic::OUTSIDE )
            ++nOutside;
        else if( cType & meshOctreeCubeBasic::DATA )
            ++nData;

        if( hasOutsideNeighbour_[leafI] )
            ++hasOutNei;
    }

    reduce(hasOutNei, sumOp<label>());
    reduce(nData, sumOp<label>());
    reduce(nOutside, sumOp<label>());
    Info << "Number of outside boxes " << nOutside << endl;
    Info << "Number of data boxes " << nData << " real data "
        << hasOutNei << endl;
    returnReduce(1, sumOp<label>());
    //::exit(1);
    # endif
}

void meshOctreeInsideOutside::markInsideCubes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();
    label nChanged;
    bool keepUpdating;
    DynList<label> neighbours;

    //- make INSIDE groups for which it is possible
    for
    (
        std::map<label, direction>::iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        const label groupI = it->first;

        if( groupI < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::UNKNOWN )
        {
            forAllRow(boundaryDATACubes_, groupI, neiI)
            {
                const label cLabel = boundaryDATACubes_(groupI, neiI);
                if(
                    hasOutsideNeighbour_[cLabel] ||
                    (
                        leaves[cLabel]->cubeType() & meshOctreeCube::INSIDE
                    )
                )
                {
                    it->second = meshOctreeCube::INSIDE;
                    break;
                }
            }
        }
    }

    do
    {
        keepUpdating = false;

        //- mark INSIDE groups created by different threads
        do
        {
            nChanged = 0;

            forAll(neighbouringGroups_, groupI)
            {
                if( groupType_[groupI] & meshOctreeCube::INSIDE )
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_(groupI, i);

                        if( groupType_[neiGroup] & meshOctreeCube::UNKNOWN )
                        {
                            ++nChanged;
                            groupType_[neiGroup] = meshOctreeCube::INSIDE;
                        }
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        if( octree.neiProcs().size() == 0 )
            continue;

        //- the code for exchanging data between different processes
        LongList<meshOctreeCubeCoordinates> dataToSend, receivedCoords;

        //- send coordinates of boxes with hasOutsideNeighbour_ flag and
        //- the boxes which have been marked as INSIDE to the neighbouring procs
        forAll(hasOutsideNeighbour_, leafI)
        {
            if(
                hasOutsideNeighbour_[leafI] ||
                (leaves[leafI]->cubeType() & meshOctreeCubeBasic::INSIDE)
            )
            {
                octree.findNeighboursForLeaf(leafI, neighbours);

                forAll(neighbours, neiI)
                    if( neighbours[neiI] == meshOctreeCube::OTHERPROC )
                    {
                        dataToSend.append(leaves[leafI]->coordinates());
                        break;
                    }
            }
        }

        octree.exchangeRequestsWithNeighbourProcessors
        (
            dataToSend,
            receivedCoords
        );

        # pragma omp parallel for if( receivedCoords.size() > 100 ) \
        private(neighbours) schedule(dynamic, 20)
        forAll(receivedCoords, i)
        {
            octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

            forAll(neighbours, neiI)
            {
                const label nei = neighbours[neiI];

                if( nei < 0 )
                    continue;

                const label groupI = cubeGroup_[nei];

                if( groupI < 0 )
                    continue;

                if( groupType_[groupI] & meshOctreeCube::UNKNOWN )
                    groupType_[groupI] = meshOctreeCubeBasic::INSIDE;
            }
        }

        do
        {
            nChanged = 0;
            dataToSend.clear();

            //- go through the list of communicationCubes and send the ones
            //- which are marked as outside
            forAll(communicationCubes_, i)
            {
                if(
                    groupType_[cubeGroup_[communicationCubes_[i]]] &
                    meshOctreeCubeBasic::INSIDE
                )
                    dataToSend.append
                    (
                        leaves[communicationCubes_[i]]->coordinates()
                    );
            }

            receivedCoords.clear();
            octree.exchangeRequestsWithNeighbourProcessors
            (
                dataToSend,
                receivedCoords
            );

            //- go through the list of received coordinates and check if any
            //- local boxes are their neighbours. If a local neighbour is
            //- a DATA box set the hasOutsideNeighbour_ flag to true. If the
            //- local neighbour is of UNKNOWN type set it to OUTSIDE.
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20) reduction(+ : nChanged)
            forAll(receivedCoords, i)
            {
                octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    const label groupI = cubeGroup_[nei];

                    if( groupI < 0 )
                        continue;

                    if( groupType_[groupI] & meshOctreeCube::UNKNOWN )
                    {
                        groupType_[groupI] = meshOctreeCube::INSIDE;
                        ++nChanged;
                    }
                }
            }

            reduce(nChanged, sumOp<label>());

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        reduce(keepUpdating, maxOp<bool>());

    } while( keepUpdating );

    //- set INSIDE type to the cubes in INSIDE groups
    for
    (
        std::map<label, direction>::const_iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        if( it->first < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::INSIDE )
        {
            const label groupI = it->first;

            //- set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_(groupI, i)]->setCubeType
                (
                    meshOctreeCube::INSIDE
                );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
