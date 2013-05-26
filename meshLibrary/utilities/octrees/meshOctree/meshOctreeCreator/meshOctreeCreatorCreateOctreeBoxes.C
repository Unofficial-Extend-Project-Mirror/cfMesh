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

#include "meshOctreeCreator.H"
#include "meshOctreeModifier.H"
#include "IOdictionary.H"
#include "boundBox.H"
#include "triSurf.H"
#include "meshOctreeAutomaticRefinement.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGOctree

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

void meshOctreeCreator::setRootCubeSizeAndRefParameters()
{
    meshOctreeModifier octreeModifier(octree_);
    if( octreeModifier.isRootInitialisedAccess() )
        return;
    if( !meshDictPtr_ )
    {
        WarningIn("void meshOctreeCreator::setRootCubeSizeAndRefParameters()")
            << "meshDict is not available" << endl;
        return;
    }

    const scalar maxSize
    (
        scalingFactor_ *
        readScalar(meshDictPtr_->lookup("maxCellSize"))
    );

    octreeModifier.searchRangeAccess() = 0.25 * maxSize;

    const triSurf& surface = octree_.surface();
    boundBox& rootBox = octreeModifier.rootBoxAccess();
    const point c = (rootBox.max() + rootBox.min()) / 2.0;

    scalar size = rootBox.max().x() - rootBox.min().x();
    if( (rootBox.max().y() - rootBox.min().y()) > size )
        size = rootBox.max().y() - rootBox.min().y();
    if( (rootBox.max().z() - rootBox.min().z()) > size )
        size = rootBox.max().z() - rootBox.min().z();

    size += 0.5 * maxSize;

    globalRefLevel_ = 0;
    bool finished;
    do
    {
        finished = false;

        const scalar lSize = size / pow(2, globalRefLevel_);

        if( lSize < (maxSize * (1.0-SMALL)) )
        {
            const scalar bbSize = 0.5 * maxSize * pow(2, globalRefLevel_);
            rootBox.max() = c + point(bbSize, bbSize, bbSize);
            rootBox.min() = c - point(bbSize, bbSize, bbSize);
            finished = true;
        }
        else
        {
            ++globalRefLevel_;
        }
    } while( !finished );

    //- exchange data
    if( Pstream::parRun() )
    {
        label l = globalRefLevel_;
        reduce(l, maxOp<label>());
        globalRefLevel_ = l;
    }

    Info << "Root box " << rootBox << endl;
    Info << "Requested cell size corresponds to octree level "
        << label(globalRefLevel_) << endl;

    octreeModifier.isRootInitialisedAccess() = true;

    surfRefLevel_ = globalRefLevel_;

    //- set other refinement parameters
    //- set boundary ref level
    if( meshDictPtr_->found("boundaryCellSize") )
    {
        direction boundaryRefLevel_ = globalRefLevel_;
        scalar cs(readScalar(meshDictPtr_->lookup("boundaryCellSize")));
        cs *= (scalingFactor_ + SMALL);

        octreeModifier.searchRangeAccess() = cs;

        label addLevel(0);
        do
        {
            finished = false;

            const scalar lSize = maxSize / pow(2, addLevel);

            if( lSize < cs )
            {
                finished = true;
            }
            else
            {
                ++addLevel;
            }
        } while( !finished );

        boundaryRefLevel_ += addLevel;

        //- exchange data
        if( Pstream::parRun() )
        {
            label l = boundaryRefLevel_;
            reduce(l, maxOp<label>());
            boundaryRefLevel_ = l;
        }

        Info << "Requested boundary cell size corresponds to octree level "
            << label(boundaryRefLevel_) << endl;

        surfRefLevel_ = boundaryRefLevel_;
    }

    //- set patch-wise ref levels
    if( meshDictPtr_->found("patchCellSize") )
    {
        patchRefinementList refPatches(meshDictPtr_->lookup("patchCellSize"));

        forAll(refPatches, patchI)
        {
            const label regionI = refPatches[patchI].patchInSurface(surface);

            if( regionI < 0 )
                continue;

            scalar cs = refPatches[patchI].cellSize();
            cs *= (scalingFactor_ + SMALL);

            label addLevel(0);
            do
            {
                finished = false;

                const scalar lSize = maxSize / pow(2, addLevel);

                if( lSize < cs )
                {
                    finished = true;
                }
                else
                {
                    ++addLevel;
                }
            } while( !finished );

            if( Pstream::parRun() )
                reduce(addLevel, maxOp<label>());

            const direction level = globalRefLevel_ + addLevel;
            forAll(surface, triI)
            {
                if(
                    (surface[triI].region() == regionI) &&
                    (surfRefLevel_[triI] < level)
                )
                    surfRefLevel_[triI] = level;
            }
        }
    }

    if( meshDictPtr_->found("subsetCellSize") )
    {
        patchRefinementList refPatches(meshDictPtr_->lookup("subsetCellSize"));

        forAll(refPatches, patchI)
        {
            const word subsetName = refPatches[patchI].patchName();

            const label subsetID = surface.facetSubsetIndex(subsetName);
            if( subsetID < 0 )
            {
                Warning << "Surface subset " << subsetName
                    << " does not exist" << endl;
                continue;
            }

            scalar cs = refPatches[patchI].cellSize();
            cs *= (scalingFactor_+ SMALL);

            label addLevel(0);
            do
            {
                finished = false;

                const scalar lSize = maxSize / pow(2, addLevel);

                if( lSize < cs )
                {
                    finished = true;
                }
                else
                {
                    ++addLevel;
                }
            } while( !finished );

            if( Pstream::parRun() )
                reduce(addLevel, maxOp<label>());

            labelListPMG subsetFaces;
            surface.facetsInSubset(subsetID, subsetFaces);
            const direction level = globalRefLevel_ + addLevel;
            forAll(subsetFaces, tI)
                if( surfRefLevel_[subsetFaces[tI]] < level )
                    surfRefLevel_[subsetFaces[tI]] = level;
        }
    }
}

void meshOctreeCreator::refineInsideAndUnknownBoxes()
{
    const direction cType = meshOctreeCube::INSIDE + meshOctreeCube::UNKNOWN;

    refineBoxes(globalRefLevel_, cType);
}

void meshOctreeCreator::createOctreeBoxes()
{
    //- set root cube size in order to achieve desired maxCellSize
    Info << "Setting cubes size and oher parameters" << endl;
    setRootCubeSizeAndRefParameters();

    //- refine to required boundary resolution
    Info << "Refining boundary" << endl;
    refineBoundary();

    //- perform automatic octree refinement
    if( !Pstream::parRun() )
    {
        Info << "Performing automatic refinement" << endl;
        meshOctreeAutomaticRefinement autoRef(octree_, *meshDictPtr_, false);

        if( hexRefinement_ )
            autoRef.activateHexRefinement();

        autoRef.automaticRefinement();
    }

    //- set up inside-outside information
    createInsideOutsideInformation();

    //- refine INSIDE and UNKNOWN boxes to to the given cell size
    refineInsideAndUnknownBoxes();

    //- refine boxes inside the geometry objects
    refineBoxesContainedInObjects();

    //- make sure that INSIDE and UNKNOWN neighbours of DATA boxes
    //- have the same or higher refinement level
    refineBoxesNearDataBoxes(1);

    //- distribute octree such that each processor has the same number
    //- of leaf boxes which will be used as mesh cells
    if( Pstream::parRun() )
    {
        loadDistribution(true);
    }

    //- delete octree data which is not needed any more
    if( Pstream::parRun() )
    {
        meshOctreeModifier om(octree_);
        om.reduceMemoryConsumption();
    }
}

void meshOctreeCreator::createOctreeWithRefinedBoundary
(
    const direction maxLevel,
    const label nTrianglesInLeaf
)
{
    const triSurf& surface = octree_.surface();
    surface.facetEdges();
    surface.edgeFacets();
    const boundBox& rootBox = octree_.rootBox();
    meshOctreeModifier octreeModifier(octree_);
    List<meshOctreeSlot>& slots = octreeModifier.dataSlotsAccess();

    while( true )
    {
        const LongList<meshOctreeCube*>& leaves =
            octreeModifier.leavesAccess();

        label nMarked(0);

        # ifdef USE_OMP
        # pragma omp parallel reduction(+ : nMarked)
        # endif
        {
            # ifdef USE_OMP
            meshOctreeSlot* slotPtr = &slots[omp_get_thread_num()];
            # else
            meshOctreeSlot* slotPtr = &slots[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 20)
            # endif
            forAll(leaves, leafI)
            {
                meshOctreeCube& oc = *leaves[leafI];

                DynList<label> ct;
                octree_.containedTriangles(leafI, ct);

                if( (oc.level() < maxLevel) && (ct.size() > nTrianglesInLeaf) )
                {
                    oc.refineCube(surface, rootBox, slotPtr);
                    ++nMarked;
                }
            }
        }

        if( nMarked == 0 )
            break;

        octreeModifier.createListOfLeaves();

    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
