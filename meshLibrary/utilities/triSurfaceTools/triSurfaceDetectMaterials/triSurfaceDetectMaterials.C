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
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectMaterials::triSurfaceDetectMaterials(const triSurf& surface)
:
    surf_(surface),
    octreePtr_(NULL),
    facePartition_(),
    nPartitions_(),
    octreeGroupForBox_(),
    nOctreeGroups_(),
    partitionMaterials_(),
    partitionType_()
{
    if( Pstream::parRun() )
        FatalError << "Material detection does not run in parallel"
            << exit(FatalError);
    
    createPartitions();
    
    createOctree();
}

triSurfaceDetectMaterials::triSurfaceDetectMaterials(meshOctree& octree)
:
    surf_(octree.surface()),
    octreePtr_(&octree),
    facePartition_(),
    nPartitions_(),
    octreeGroupForBox_(),
    nOctreeGroups_(),
    partitionMaterials_(),
    partitionType_()
{
    createPartitions();
}

triSurfaceDetectMaterials::~triSurfaceDetectMaterials()
{
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectMaterials::detectMaterialsAndInternalWalls()
{
    if( Pstream::parRun() )
        FatalError << "Material detection does not run in parallel"
            << exit(FatalError);
    
    bool finished;
    do
    {
        finished = true;
        
        findOctreeGroups();
        
        findMaterialsAndWalls();
        
        if( checkMaterials() )
        {
            finished = false;
            refineOctree();
        }
    } while( !finished );
}

void triSurfaceDetectMaterials::writeMaterials(const fileName& fName) const
{
    forAll(partitionMaterials_, partI)
    {
        Info << "Partition " << partI << " is in materials "
            << partitionMaterials_[partI] << endl;
        
        if( partitionMaterials_.sizeOfRow(partI) != 1 )
            continue;
        
        if( partitionType_[partI] == meshOctreeCubeBasic::INSIDE )
            Info << "Partition " << partI << " is a baffle" << endl;
        if( partitionType_[partI] == meshOctreeCubeBasic::UNKNOWN )
            Info << "Partition " << partI << " is an enclosed baffle" << endl;
    }
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
