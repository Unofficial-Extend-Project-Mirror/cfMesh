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
    tetMeshOptimisation

Description
    Mesh smoothing without any topological changes

SourceFiles
    tetMeshOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshOptimisation_H
#define tetMeshOptimisation_H

#include "DynList.H"
#include "parPartTet.H"
#include "boundBox.H"
#include "labelListPMG.H"
#include "boolList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMesh;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                    Class tetMeshOptimisation Declaration
\*---------------------------------------------------------------------------*/

class tetMeshOptimisation
{
    // Private data
        //- reference to the tet mesh
        partTetMesh& tetMesh_;
    
    // Private member functions needed for parallel runs
        //- make sure that all processors have the same points
        //- marked as negativeNode
        void unifyNegativePoints(boolList& negativeNode) const;
    
        //- exchange data with other processors
        void exchangeData
        (
            std::map<label, DynList<parPartTet> >& m,
            boolList* negativeNodePtr = NULL
        );
    
        //- update buffer layer points
        void updateBufferLayerPoints();
    
        //- make sure that coordinates of moved points remain the same
        //- on all processor containing those points
        void unifyCoordinatesParallel(const boolList* negativeNodePtr = NULL);

public:

    // Constructors

        //- Construct from tet mesh
        tetMeshOptimisation(partTetMesh& mesh);


    // Destructor

        ~tetMeshOptimisation();

    // Member Functions
        //- untangle mesh by using Patrik Knupp's simple metric
        void optimiseUsingKnuppMetric();

        //- smooth using mesh untangler
        void optimiseUsingMeshUntangler();

        //- smooth using volume optimizer
        void optimiseUsingVolumeOptimizer();

        //- smooth boundary using the volume optimizer
        void optimiseBoundaryVolumeOptimizer
        (
            const bool nonShrinking = false
        );
        
        //- smooth boundary using shrinking surface laplace
        void optimiseBoundarySurfaceLaplace();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
