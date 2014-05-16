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

Application
    Test of the cartesian mesh

Description
    - creates an octree and creates cartesian mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshOctreeCreator.H"
#include "meshOctreeAutomaticRefinement.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "cartesianMeshExtractor.H"
#include "triSurf.H"
#include "helperFunctions.H"
#include "findCellsIntersectingSurface.H"

#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- find octree leaves connected over faces
class octreeNeighbours
{
    const meshOctree& octree_;

public:

    octreeNeighbours(const meshOctree& octree)
    :
      octree_(octree)
    {}

    label size() const
    {
        return octree_.numberOfLeaves();
    }

    void operator()(const label leafI, DynList<label>& neighbours) const
    {
        octree_.findNeighboursForLeaf(leafI, neighbours);
    }
};

class octreeSelectOperator
{
    const meshOctree& octree_;

public:

    octreeSelectOperator(const meshOctree& octree)
    :
      octree_(octree)
    {}

    bool operator()(const label leafI) const
    {
        if( octree_.returnLeaf(leafI).cubeType() & meshOctreeCube::DATA )
            return false;

        return true;
    }
};

class meshNeighbourOperator
{
    const polyMeshGen& mesh_;

public:

    meshNeighbourOperator(const polyMeshGen& mesh)
    :
        mesh_(mesh)
    {}

    label size() const
    {
        return mesh_.cells().size();
    }

    void operator()(const label cellI, DynList<label>& neighbourCells) const
    {
        neighbourCells.clear();

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const cell& c = mesh_.cells()[cellI];

        forAll(c, fI)
        {
            label nei = owner[c[fI]];

            if( nei == cellI )
                nei = neighbour[c[fI]];

            if( nei >= 0 )
                neighbourCells.append(nei);
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >& neiGroups,
        const labelListType& elementInGroup,
        const DynList<label>& localGroupLabel
    ) const
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        const labelList& owner = mesh_.owner();

        //- send the data to other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelList groupOwner(procBoundaries[patchI].patchSize());
            for(label faceI=0;faceI<size;++faceI)
            {
                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                {
                    groupOwner[faceI] = -1;
                    continue;
                }

                groupOwner[faceI] = localGroupLabel[groupI];
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                groupOwner.byteSize()
            );

            toOtherProc << groupOwner;
        }

        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            forAll(receivedData, faceI)
            {
                if( receivedData[faceI] < 0 )
                    continue;

                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                    continue;

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(receivedData[faceI]);
            }
        }
    }
};

class meshSelectorOperator
{
    const VRWGraph& patchesIntersectingCell_;

public:

    meshSelectorOperator(const VRWGraph& patchesIntersectingCell)
    :
        patchesIntersectingCell_(patchesIntersectingCell)
    {}

    bool operator()(const label cellI) const
    {
        if( patchesIntersectingCell_[cellI].size() )
            return false;

        return true;
    }
};

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surf;
    surf.readSurface(runTime.path()/surfaceFile);

    // construct the octree
    meshOctree mo(surf);
    meshOctreeCreator(mo).createOctreeWithRefinedBoundary(10, 30);

    //- create and read the mesh
    polyMeshGen pmg(runTime);
    pmg.read();

    //- find cells intersecting surface
    findCellsIntersectingSurface fis(pmg, mo);

    meshNeighbourOperator mnop(pmg);
    meshSelectorOperator msop(fis.facetsIntersectingCells());

    labelLongList result;
    help::frontalMarking(result, 0, mnop, msop);
    Info << "Cells in the group " << result.size() << endl;

    const label nGroups = help::groupMarking(result, mnop, msop);

    Info << "Number of groups " << nGroups << endl;
    labelList groupId(nGroups);
    forAll(groupId, i)
        groupId[i] = pmg.addCellSubset("group_"+help::scalarToText(i));

    forAll(result, cellI)
    {
        if( result[cellI] < 0 )
            continue;

        pmg.addCellToSubset(groupId[result[cellI]], cellI);
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
