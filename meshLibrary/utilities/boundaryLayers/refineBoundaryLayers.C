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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& refineBoundaryLayers::surfaceEngine() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

refineBoundaryLayers::refineBoundaryLayers(polyMeshGen& mesh)
:
    mesh_(mesh),
    msePtr_(NULL),
    globalNumLayers_(1),
    globalThicknessRatio_(1.0),
    globalMaxThicknessFirstLayer_(-1.0),
    numLayersForPatch_(),
    thicknessRatioForPatch_(),
    maxThicknessForPatch_(),
    discontinuousLayersForPatch_(),
    done_(false),
    nLayersAtBndFace_(),
    splitEdges_(),
    splitEdgesAtPoint_(),
    newVerticesForSplitEdge_(),
    facesFromFace_(),
    newFaces_()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

refineBoundaryLayers::~refineBoundaryLayers()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::setGlobalNumberOfLayers(const label nLayers)
{
    if( nLayers < 2 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalNumberOfLayers(const label)"
        ) << "The specified global number of boundary layers is less than 2"
          << endl;

        return;
    }

    globalNumLayers_ = nLayers;
}

void refineBoundaryLayers::setGlobalThicknessRatio(const scalar thicknessRatio)
{
    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalThicknessRatio(const scalar)"
        ) << "The specified global thickness ratio is less than 1.0" << endl;

        return;
    }

    globalThicknessRatio_ = thicknessRatio;
}

void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer
(
    const scalar maxThickness
)
{
    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const scalar)"
        ) << "The specified global maximum thickness of the first"
          << " boundary layer is negative!!" << endl;

        return;
    }

    globalMaxThicknessFirstLayer_ = maxThickness;
}

void refineBoundaryLayers::setNumberOfLayersForPatch
(
    const word& patchName,
    const label nLayers
)
{
    if( nLayers < 2 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setNumberOfLayersForPatch"
            "(const word&, const label)"
        ) << "The specified number of boundary layers for patch " << patchName
          << " is less than 2" << endl;

        return;
    }

    numLayersForPatch_[patchName] = nLayers;
}

void refineBoundaryLayers::setThicknessRatioForPatch
(
    const word& patchName,
    const scalar thicknessRatio
)
{
    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setThicknessRatioForPatch"
            "(const word&, const scalar)"
        ) << "The specified thickness ratio for patch " << patchName
          << " is less than 1.0" << endl;

        return;
    }

    thicknessRatioForPatch_[patchName] = thicknessRatio;
}

void refineBoundaryLayers::setMaxThicknessOfFirstLayerForPatch
(
    const word& patchName,
    const scalar maxThickness
)
{
    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const word&, const scalar)"
        ) << "The specified maximum thickness of the first boundary layer "
          << "for patch " << patchName << " is negative!!" << endl;

        return;
    }

    maxThicknessForPatch_[patchName] = maxThickness;
}

void refineBoundaryLayers::setInteruptForPatch(const word& patchName)
{
    discontinuousLayersForPatch_.insert(patchName);
}

void refineBoundaryLayers::refineLayers()
{
    Info << "Starting refining boundary layers" << endl;

    if( done_ )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineLayers()"
        ) << "Boundary layers are already refined! Stopping refinement" << endl;

        return;
    }

    analyseLayers();

    if( !findSplitEdges() )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineLayers()"
        ) << "Boundary layers do not exist in the mesh! Cannot refine" << endl;

        return;
    }

    generateNewVertices();

    generateNewFaces();

    generateNewCells();

    done_ = true;

    Info << "Finished refining boundary layers" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
