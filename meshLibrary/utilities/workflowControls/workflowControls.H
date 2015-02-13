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

Class
    workflowControls

Description
    A controller fo managing the mesh generation workflow. It allows the user
    to stop the meshing process after selected steps in the workflow,
    and restart from a selected point.

SourceFiles
    workflowControls.C

\*---------------------------------------------------------------------------*/

#ifndef workflowControls_H
#define workflowControls_H

#include "objectRegistry.H"
#include "Time.H"
#include "IOdictionary.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMeshGen;

/*---------------------------------------------------------------------------*\
                    Class workflowControls Declaration
\*---------------------------------------------------------------------------*/

class workflowControls
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- status of the workflow (sets bits to completed steps)
        label status_;

        //- current step in the workflow
        word currentStep_;

    // static private data
        static const std::map<word, label> workflowSteps_;

    // Private member functions
        //- check if it is possible to restart after the current step
        bool restartPossibleAfterCurrentStep() const;

        //- sets the current step to completed
        void setStepCompleted() const;

        //- shall the procedure stop after the current step
        bool exitAfterCurrentStep() const;

    // Static private member functions
        //- populate workflowSteps with values
        static std::map<word, label> populateWorkflowSteps();

public:

    // Constructors

        //- Construct from IOdictionary
        workflowControls(polyMeshGen& mesh);

    // Destructor
        ~workflowControls();

    // Public member functions

        //- set current step
        void setCurrentStep(const word&);

        //- shall the workflow stop after the current step
        bool stopAfterCurrentStep() const;

        //- shall the workflow restart after the current step
        bool restartAfterCurrentStep() const;

        //- set the workflow completed flag
        void workflowCompleted();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
