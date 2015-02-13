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

#include "workflowControls.H"
#include "polyMeshGen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const std::map<word, label> workflowControls::workflowSteps_ =
    populateWorkflowSteps();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool workflowControls::restartPossibleAfterCurrentStep() const
{
    return false;
}

void workflowControls::setStepCompleted() const
{
    mesh_.metaData().add("lastStep", currentStep_);
}

bool workflowControls::exitAfterCurrentStep() const
{
    const IOdictionary& meshDict =
        mesh_.returnTime().lookupObject<IOdictionary>("meshDict");

    const dictionary& workflowControls = meshDict.subDict("workflowControls");

    if( workflowControls.found("stopAfter") )
    {
        const word exitStep(workflowControls.lookup("stopAfter"));

        if( exitStep == currentStep_ )
            return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::map<word, label> workflowControls::populateWorkflowSteps()
{
    std::map<word, label> workflowSteps;
    workflowSteps.insert(std::make_pair(word("templateGeneration"), 1));
    workflowSteps.insert(std::make_pair(word("surfaceTopology"), 2));
    workflowSteps.insert(std::make_pair(word("surfaceProjection"), 4));
    workflowSteps.insert(std::make_pair(word("patchAssignment"), 8));
    workflowSteps.insert(std::make_pair(word("edgeExtraction"), 16));
    workflowSteps.insert(std::make_pair(word("boundaryLayerGeneration"), 32));
    workflowSteps.insert(std::make_pair(word("boundaryLayerRefinement"), 64));

    return workflowSteps;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

workflowControls::workflowControls(polyMeshGen& mesh)
:
    mesh_(mesh),
    status_(0),
    currentStep_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

workflowControls::~workflowControls()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void workflowControls::setCurrentStep(const word& stepName)
{
    //- check if the requested step exists in the database of steps
    std::map<word, label>::const_iterator it = workflowSteps_.find(stepName);
    if( it == workflowSteps_.end() )
    {
        FatalErrorIn
        (
            "void workflowControls::setCurrentStep(const word&)"
        ) << "Step " << stepName << " is not a valid name." << exit(FatalError);
    }

    currentStep_ = stepName;
    status_ |= it->second;
}

bool workflowControls::stopAfterCurrentStep() const
{
    setStepCompleted();

    if( exitAfterCurrentStep() )
    {
        bool writeSuccess(true);

        try
        {
            mesh_.write();
        }
        catch(...)
        {
            writeSuccess = false;
        }

        returnReduce(writeSuccess, minOp<bool>());

        if( !writeSuccess )
            FatalErrorIn
            (
                "bool workflowControls::stopAfterCurrentStep() const"
            ) << "Mes was not written on disk" << exit(FatalError);

        return true;
    }

    return false;
}

bool workflowControls::restartAfterCurrentStep() const
{
    return false;
}

void workflowControls::workflowCompleted()
{
    mesh_.metaData().remove("lastStep");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
