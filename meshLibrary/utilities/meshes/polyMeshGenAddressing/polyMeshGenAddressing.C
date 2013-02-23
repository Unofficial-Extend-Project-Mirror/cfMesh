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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyMeshGenAddressing, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyMeshGenAddressing::polyMeshGenAddressing(const polyMeshGenCells& mesh)
:
    mesh_(mesh),
    edgesPtr_(NULL),
    ccPtr_(NULL),
    ecPtr_(NULL),
    pcPtr_(NULL),
    efPtr_(NULL),
    pfPtr_(NULL),

    cePtr_(NULL),
    fePtr_(NULL),
    pePtr_(NULL),
    ppPtr_(NULL),
    cpPtr_(NULL),
    cellCentresPtr_(NULL),
    faceCentresPtr_(NULL),
    cellVolumesPtr_(NULL),
    faceAreasPtr_(NULL),
    globalPointLabelPtr_(NULL),
    globalFaceLabelPtr_(NULL),
    globalCellLabelPtr_(NULL),
    globalEdgeLabelPtr_(NULL),
    pProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    pointNeiProcsPtr_(NULL),
    eProcsPtr_(NULL),
    globalToLocalEdgeAddressingPtr_(NULL),
    edgeNeiProcsPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMeshGenAddressing::~polyMeshGenAddressing()
{
    clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
