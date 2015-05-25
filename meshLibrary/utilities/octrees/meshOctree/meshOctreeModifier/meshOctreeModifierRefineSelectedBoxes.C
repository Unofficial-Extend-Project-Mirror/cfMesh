/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "meshOctreeModifier.H"
#include "triSurf.H"
#include "HashSet.H"
#include "helperFunctions.H"
#include "meshOctreeCubeCoordinatesScalar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#include <set>

#include <sys/stat.h>

//#define OCTREETiming
//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::markAdditionalLayers
(
    labelList& refineBox,
    const label nLayers
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for parallel runs to reduce the communication messages
    labelHashSet transferCoordinates;

    DynList<label> neiLeaves;

    for(label i=1;i<=nLayers;++i)
    {
        LongList<meshOctreeCubeCoordinates> processorChecks;

        transferCoordinates.clear();

        labelLongList activeLeaves;
        forAll(leaves, leafI)
            if( refineBox[leafI] == i )
                activeLeaves.append(leafI);

        # ifdef USE_OMP
        # pragma omp parallel for private(neiLeaves) schedule(dynamic, 20)
        # endif
        forAll(activeLeaves, lI)
        {
            const label leafI = activeLeaves[lI];

            const meshOctreeCubeCoordinates& oc = leaves[leafI]->coordinates();

            neiLeaves.clear();
            octree_.findAllLeafNeighbours(oc, neiLeaves);

            forAll(neiLeaves, posI)
            {
                const label neiLabel = neiLeaves[posI];

                if( neiLabel == meshOctreeCubeBasic::OTHERPROC )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        if( !transferCoordinates.found(leafI) )
                        {
                            processorChecks.append(oc);
                            transferCoordinates.insert(leafI);
                        }
                    }

                    continue;
                }

                if( neiLabel < 0 )
                    continue;

                if( !refineBox[neiLabel] )
                    refineBox[neiLabel] = i+1;
            }
        }

        if( octree_.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree_.exchangeRequestsWithNeighbourProcessors
            (
                processorChecks,
                receivedCoords
            );

            //- check consistency with received cube coordinates
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 1000 ) \
            schedule(dynamic, 20) private(neiLeaves)
            # endif
            forAll(receivedCoords, ccI)
            {
                octree_.findAllLeafNeighbours(receivedCoords[ccI], neiLeaves);

                forAll(neiLeaves, posI)
                {
                    if( neiLeaves[posI] < 0 )
                        continue;

                    if( !refineBox[neiLeaves[posI]] )
                        refineBox[neiLeaves[posI]] = i+1;
                }
            }
        }
    }
}

void meshOctreeModifier::markAdditionalLayersOfFaceNeighbours
(
    labelList& refineBox,
    const label nLayers
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for parallel runs to reduce the communication messages
    labelHashSet transferCoordinates;

    DynList<label> neiLeaves;

    for(label i=1;i<=nLayers;++i)
    {
        LongList<meshOctreeCubeCoordinates> processorChecks;

        transferCoordinates.clear();

        labelLongList activeLeaves;
        forAll(leaves, leafI)
            if( refineBox[leafI] == i )
                activeLeaves.append(leafI);

        # ifdef USE_OMP
        # pragma omp parallel for private(neiLeaves) schedule(dynamic, 20)
        # endif
        forAll(activeLeaves, lI)
        {
            const label leafI = activeLeaves[lI];

            const meshOctreeCubeCoordinates& oc = leaves[leafI]->coordinates();

            neiLeaves.clear();
            octree_.findNeighboursForLeaf(oc, neiLeaves);

            forAll(neiLeaves, posI)
            {
                const label neiLabel = neiLeaves[posI];

                if( neiLabel == meshOctreeCubeBasic::OTHERPROC )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        if( !transferCoordinates.found(leafI) )
                        {
                            processorChecks.append(oc);
                            transferCoordinates.insert(leafI);
                        }
                    }

                    continue;
                }

                if( neiLabel < 0 )
                    continue;

                if( !refineBox[neiLabel] )
                    refineBox[neiLabel] = i+1;
            }
        }

        if( octree_.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree_.exchangeRequestsWithNeighbourProcessors
            (
                processorChecks,
                receivedCoords
            );

            //- check consistency with received cube coordinates
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 1000 ) \
            schedule(dynamic, 20) private(neiLeaves)
            # endif
            forAll(receivedCoords, ccI)
            {
                neiLeaves.clear();
                octree_.findNeighboursForLeaf(receivedCoords[ccI], neiLeaves);

                forAll(neiLeaves, posI)
                {
                    if( neiLeaves[posI] < 0 )
                        continue;

                    if( !refineBox[neiLeaves[posI]] )
                        refineBox[neiLeaves[posI]] = i+1;
                }
            }
        }
    }
}

# ifdef DEBUGSearch
void writeLeaves
(
    const fileName& fName,
    const meshOctree& octree,
    const labelList& markedBoxes,
    const label layer
)
{
    labelLongList activeLeaves;

    forAll(markedBoxes, leafI)
        if( markedBoxes[leafI] == layer )
            activeLeaves.append(leafI);

    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << (8 * activeLeaves.size()) << " float\n";
    forAll(activeLeaves, i)
    {
        const label leafI = activeLeaves[i];
        FixedList<point, 8> vertices;
        octree.returnLeaf(leafI).vertices(octree.rootBox(), vertices);

        forAll(vertices, vI)
        {
            const point& p = vertices[vI];

            file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
    }

    //- write lines
    file << "\nPOLYGONS " << (6*activeLeaves.size())
         << " " << 30*activeLeaves.size() << nl;
    forAll(activeLeaves, i)
    {
        const label startNode = 8 * i;
        for(label fI=0;fI<6;++fI)
        {
            file << 4;

            for(label pI=0;pI<4;++pI)
                file << " " << (startNode+meshOctreeCube::faceNodes_[fI][pI]);

            file << nl;
        }
    }

    file << "\n";
}
# endif

label meshOctreeModifier::markAdditionalLayers
(
    labelList& refineBox,
    labelList& nLayers,
    List<direction>& targetRefLevel
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    # ifdef DEBUGSearch
    Info << "Marking additional layers " << endl;
    # endif

    //- sort leaves based on the number of additional layers
    label maxLevel = Foam::max(nLayers);
    reduce(maxLevel, maxOp<label>());

    List<labelLongList> leavesForLayer(maxLevel+1);
    forAll(nLayers, leafI)
    {
        if( nLayers[leafI] < 1 )
            continue;

        leavesForLayer[nLayers[leafI]].append(leafI);
    }

    //- set refinement flag to additional boxes marked separately
    label nMarked(0);

    forAllReverse(leavesForLayer, layerI)
    {
        const labelLongList& activeLeaves = leavesForLayer[layerI];

        if( returnReduce(activeLeaves.size(), sumOp<label>()) == 0 )
            continue;

        //- find the max required refinement level
        direction maxLevel(0);

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            direction localMax(0);

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 50)
            # endif
            forAll(activeLeaves, i)
                localMax = Foam::max(localMax, targetRefLevel[activeLeaves[i]]);

            # ifdef USE_OMP
            # pragma omp critical
            # endif
            maxLevel = Foam::max(localMax, maxLevel);
        }

        label ml = maxLevel;
        reduce(ml, maxOp<label>());
        maxLevel = ml;

        //- mark additional boxes for refinement
        for(direction levelI=maxLevel;levelI>0;--levelI)
        {
            labelList markedBoxes(leaves.size(), 0);

            label counter(0);
            # ifdef USE_OMP
            # pragma omp parallel for reduction(+:counter)
            # endif
            forAll(activeLeaves, lI)
            {
                const label leafI = activeLeaves[lI];

                if( targetRefLevel[leafI] == levelI )
                {
                    markedBoxes[leafI] = 1;
                    ++counter;
                }
            }

            if( returnReduce(counter, sumOp<label>()) == 0 )
                continue;

            //- mark additional cells at this refinement level
            markAdditionalLayersOfFaceNeighbours(markedBoxes, layerI);
            //markAdditionalLayersInRange(markedBoxes, 150.0);

/*
            label minLayer(layerI+1);
            forAll(markedBoxes, leafI)
            {
                if( markedBoxes[leafI] < 2 )
                    continue;

                if( leaves[leafI]->level() < levelI )
                    minLayer = Foam::min(minLayer, markedBoxes[leafI]);
            }

            Info << "min layer " << minLayer << end
*/;

/*            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 20)
            # endif
            forAll(markedBoxes, leafI)
            {
                if( markedBoxes[leafI] & 1)
                {
                    const point c = leaves[leafI]->centre(octree_.rootBox());
                    const scalar r =
                     nLayers[leafI] * leaves[leafI]->size(octree_.rootBox());

                    bool atProcessorBnd(false);
                    octree_.initialCubePtr_->markLeavesInSphere
                    (
                        octree_.rootBox(),
                        c,
                        r,
                        markedBoxes,
                        atProcessorBnd
                    );
                }
            }
*/
            # ifdef DEBUGSearch
            for(label i=1;i<(layerI+1);++i)
            {
                const fileName fName("leaves_"+help::labelToText(i)+".vtk");
                writeLeaves(fName, octree_, markedBoxes, i);
            }

            Info << "LayerI " << layerI << endl;
            # endif

            //- update the main list
            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 100) \
            reduction(+:nMarked)
            # endif
            forAll(markedBoxes, leafI)
            {
                if( markedBoxes[leafI] < 2 )
                    continue;
//                if( markedBoxes[leafI] > minLayer)
//                    continue;
                if( leaves[leafI]->level() >= levelI )
                    continue;

//                label currLayer = markedBoxes[leafI]-1;
//                direction currLevel = leaves[leafI]->level();
//                while( ++currLevel < levelI )
//                {
//                    currLayer *= 2;
//                }
//                --currLayer;

//                if( currLayer > layerI )
//                    continue;

                if( !refineBox[leafI] )
                {
                    refineBox[leafI] |= 1;
                    ++nMarked;
                }
            }
        }
    }

    reduce(nMarked, sumOp<label>());

    return nMarked;
}

void meshOctreeModifier::refineSelectedBoxes
(
    labelList& refineBox,
    const bool hexRefinement
)
{
    # ifdef OCTREETiming
    const scalar startTime = omp_get_wtime();
    # endif

    //- ensure that refinement will produce 1-irregular octree
    do
    {
        ensureCorrectRegularity(refineBox);
    } while( hexRefinement && ensureCorrectRegularitySons(refineBox) );

    # ifdef OCTREETiming
    const scalar regTime = omp_get_wtime();
    Info << "Time for ensuring regularity " << (regTime-startTime) << endl;
    # endif

    const triSurf& surface = octree_.surface();
    const boundBox& rootBox = octree_.rootBox();
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for thread safety
    //- such solutions make me a sad bunny :(
    surface.facetEdges();
    surface.edgeFacets();
    surface.edges();

    # ifdef USE_OMP
    # pragma omp parallel num_threads(octree_.dataSlots_.size())
    # endif
    {
        # ifdef USE_OMP
        meshOctreeSlot* slotPtr = &octree_.dataSlots_[omp_get_thread_num()];
        # else
        meshOctreeSlot* slotPtr = &octree_.dataSlots_[0];
        # endif

        if( !octree_.isQuadtree() )
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(leaves, leafI)
            {
                if( refineBox[leafI] )
                    leaves[leafI]->refineCube(surface, rootBox, slotPtr);
            }
        }
        else
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(leaves, leafI)
            {
                if( refineBox[leafI] )
                    leaves[leafI]->refineCube2D(surface, rootBox, slotPtr);
            }
        }
    }

    createListOfLeaves();

    # ifdef OCTREETiming
    Info << "Time for actual refinement " << (omp_get_wtime()-regTime) << endl;
    # endif
}

void meshOctreeModifier::refineSelectedBoxesAndAdditionalLayers
(
    labelList& refineBox,
    const scalarList& refThickness,
    const List<direction>& targetRefLevel
)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    typedef std::map<direction, LongList<meshOctreeCube*> > lMap;
    lMap leavesMap;

    # ifdef DEBUGSearch
    Info << "Marking additional layers " << endl;
    # endif

    //- sort leaves based on the number of additional layers
    label maxLevel = Foam::max(targetRefLevel);
    reduce(maxLevel, maxOp<label>());

    List<labelLongList> leavesForLevel(maxLevel+1);
    forAll(refineBox, i)
    {
        if( targetRefLevel[i] < 1 )
            continue;

        leavesForLevel[targetRefLevel[i]].append(i);
    }

    forAllReverse(leavesForLevel, levelI)
    {
        const labelLongList& activeLeaves = leavesForLevel[levelI];

        if( returnReduce(activeLeaves.size(), sumOp<label>()) == 0 )
            continue;

        LongList<meshOctreeCube*>& currLeaves = leavesMap[levelI];

        forAll(activeLeaves, i)
            currLeaves.append(leaves[activeLeaves[i]]);
    }

    forAllConstIter(lMap, leavesMap, it)
    {
        const direction levelI = it->first;

        const LongList<meshOctreeCube*>& selectedLeaves = it->second;
        const labelLongList& leafLabels = leavesForLevel[levelI];

        label nMarked;
        do
        {
            nMarked = 0;

            labelList markedLeaves(leaves.size(), 0);
            DynList<label> neiLeaves;

            LongList<meshOctreeCubeCoordinates> procCheck;
            LongList<scalar> procDistance;

            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 10) private(neiLeaves)
            # endif
            forAll(selectedLeaves, i)
            {
                const meshOctreeCube& oc = *selectedLeaves[i];
                const point c = oc.centre(octree_.rootBox());

                const scalar range = refThickness[leafLabels[i]];

                neiLeaves.clear();
                octree_.findLeavesInSphere
                (
                    c,
                    range,
                    neiLeaves
                );

                forAll(neiLeaves, j)
                {
                    const label neiLabel = neiLeaves[j];

                    if( neiLabel < 0 )
                    {
                        if( neiLabel == meshOctreeCube::OTHERPROC )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical
                            # endif
                            {
                                procCheck.append(oc.coordinates());
                                procDistance.append(range);
                            }
                        }

                        continue;
                    }

                    if( leaves[neiLabel]->level() >= levelI )
                        continue;

                    if( leaves[neiLabel]->level() <= oc.level() )
                    {
                        markedLeaves[neiLabel] = 1;
                        nMarked = 1;
                    }
                }
            }

            if( octree_.neiProcs().size() )
            {
                LongList<meshOctreeCubeCoordinates> receivedData;
                LongList<scalar> receivedDistance;

                octree_.exchangeRequestsWithNeighbourProcessors
                (
                    procCheck,
                    procDistance,
                    receivedData,
                    receivedDistance
                );

                returnReduce(1, sumOp<label>());

                # ifdef USE_OMP
                # pragma omp parallel for schedule(dynamic, 5) \
                private(neiLeaves)
                # endif
                forAll(receivedData, i)
                {
                    const meshOctreeCubeCoordinates& cc = receivedData[i];
                    const scalar range = receivedDistance[i];

                    const point c = cc.centre(octree_.rootBox());

                    octree_.findLeavesInSphere
                    (
                        c,
                        range,
                        neiLeaves
                    );

                    forAll(neiLeaves, j)
                    {
                        const label neiLabel = neiLeaves[j];

                        if( neiLabel < 0 )
                            continue;

                        if( leaves[neiLabel]->level() >= levelI )
                            continue;

                        if( leaves[neiLabel]->level() <= cc.level() )
                        {
                            markedLeaves[neiLabel] = 1;
                            nMarked = 1;
                        }
                    }
                }
            }

            reduce(nMarked, sumOp<label>());

            if( nMarked == 0 )
                break;

            refineSelectedBoxes(markedLeaves);
        } while( nMarked );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
