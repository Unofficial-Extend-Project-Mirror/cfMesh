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

#include "tetMeshExtractor.H"
#include "tetTessellation.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "polyMeshGenModifierAddCellByCell.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshExtractor::selectElements()
{
    const LongList<point>& tetPoints = tessellation_.points();
    const LongList<tessellationElement>& elmts = tessellation_.elmts();

    boolList useBox ( octree_.numberOfLeaves(), false );

    boolList usePatch ( octree_.surface().patches().size(), false );
    forAll ( patchesForBoundaryTet_, bndI )
    {
        usePatch[patchesForBoundaryTet_[bndI]] = true;
    }

    forAll ( useBox, boxI )
    if ( octree_.returnLeaf ( boxI ).cubeType() & meshOctreeCubeBasic::INSIDE )
    {
        useBox[boxI] = true;
    }
    else if (
        octree_.returnLeaf ( boxI ).cubeType() & meshOctreeCubeBasic::DATA
    )
    {
        if ( useBoundaryTets_ )
            useBox[boxI] = true;

        DynList<label> patches;
        octree_.findBoundaryPatchesForLeaf ( boxI,  patches );
        forAll ( patches, patchI )
        if ( usePatch[patches[patchI]] )
            useBox[boxI] = true;
    }

    boolList internalPoint ( tetPoints.size(), false );

    forAll ( tetPoints, pointI )
    {
        const label leafI =
            octree_.findLeafContainingVertex ( tetPoints[pointI] );

        if ( ( leafI != -1 ) && useBox[leafI] )
            internalPoint[pointI] = true;
    }

    forAll ( elmts, elmtI )
    {
        bool keep ( true );
        const tessellationElement& elmt = elmts[elmtI];
        for ( direction i=0;i<DIM1;++i )
            if ( !internalPoint[elmt[i]] )
            {
                keep = false;
                break;
            }

        if ( keep )
            useElement_[elmtI] = true;
    }

    useElement_ = true;
}

void tetMeshExtractor::createPoints()
{
    polyMeshGenModifier meshModifier ( mesh_ );
    pointFieldPMG& points = meshModifier.pointsAccess();

    const LongList<point>& tetPoints = tessellation_.points();
    points.setSize ( tetPoints.size() );

    forAll ( tetPoints, pointI )
    points[pointI] = tetPoints[pointI];
}

void tetMeshExtractor::createPolyMesh()
{
    faceList cellFaces ( 4, face ( 3 ) );

    const LongList<tessellationElement>& elmts = tessellation_.elmts();

    polyMeshGenModifierAddCellByCell meshModifier ( mesh_ );

    forAll ( useElement_, elmtI )
    if ( useElement_[elmtI] )
    {
        const tessellationElement& elmt = elmts[elmtI];

        for ( label i=0;i<4;++i )
        {
            face& f = cellFaces[i];
            triFace tf = elmt.face ( i );
            f[0] = tf[0];
            f[1] = tf[2];
            f[2] = tf[1];
        }

        meshModifier.addCell ( cellFaces );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
