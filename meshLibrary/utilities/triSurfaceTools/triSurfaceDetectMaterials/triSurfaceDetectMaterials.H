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

Class
    triSurfaceDetectMaterials

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfaceDetectMaterials.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDetectMaterials_H
#define triSurfaceDetectMaterials_H

#include "triSurf.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// Forward declaration
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class triSurfaceDetectMaterials Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDetectMaterials
{
    // Private data
        //- reference to triSurf
        const triSurf& surf_;
    
        //- pointer to the octree
        meshOctree* octreePtr_;
    
        //- partition for each surface facet
        labelList facePartition_;
        
        //- number of partitions found in the surface
        label nPartitions_;
    
        //- octree groups
        labelList octreeGroupForBox_;
    
        //- number of octree groups
        label nOctreeGroups_;
    
        //- graph of face partitions belonging to each material
        VRWGraph partitionMaterials_;
        List<direction> partitionType_;

    // Private member functions
        //- create partitioning
        void createPartitions();
        
        //- create the octree
        void createOctree();
        
        //- refine octree locally if needed
        void refineOctree();
        
        //- find octree groups
        void findOctreeGroups();
        
        //- find materials and zero-thickness walls
        void findMaterialsAndWalls();
        
        //- check if the partitions are well behaved
        bool checkMaterials() const;

        //- Disallow default bitwise copy construct
        triSurfaceDetectMaterials(const triSurfaceDetectMaterials&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceDetectMaterials&);

public:

    // Constructors

        //- Construct from triSurface
        triSurfaceDetectMaterials(const triSurf& surface);

        //- Construct from meshOctree
        triSurfaceDetectMaterials(meshOctree& octree);

    // Destructor

        ~triSurfaceDetectMaterials();


    // Member Functions
        
        //- detect surface parts belonging to different materials
        void detectMaterialsAndInternalWalls();

        //- return the number of partitions
        inline label numberOfPartitions() const;

        //- return partitions
        inline const labelList& facePartitions() const;
        
        //- return the number of domains found in the surface
        inline label numberOfDomains() const;
        
        //- return domains that partitions belongs to
        inline const VRWGraph& partitionDomains() const;
        
        //- return partition type
        inline const List<direction>& partitionType() const;

        //- write result to a file
        void writeMaterials(const fileName&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfaceDetectMaterialsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
