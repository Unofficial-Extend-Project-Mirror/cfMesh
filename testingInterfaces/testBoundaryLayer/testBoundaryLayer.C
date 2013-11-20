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

Application
    Test for boundary layers

Description
    -

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMeshGen.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "polyMeshGenChecks.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

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

    //boundaryLayers bndLayers(pmg);
    //bndLayers.addLayerForPatch("inlet");
    //bndLayers.addLayerForPatch("symmetryplane");
    //bndLayers.createOTopologyLayers();
    //bndLayers.addLayerForAllPatches();

    Info << "Starting bnd layer refinement "
         << runTime.elapsedClockTime() << endl;

    refineBoundaryLayers refLayers(pmg);

    if( meshDict.isDict("boundaryLayers") )
    {
        const dictionary& bndLayers = meshDict.subDict("boundaryLayers");

        //- read global properties
        if( bndLayers.found("nLayers") )
        {
            const label nLayers = readLabel(bndLayers.lookup("nLayers"));
            refLayers.setGlobalNumberOfLayers(nLayers);
        }
        if( bndLayers.found("thicknessRatio") )
        {
            const scalar ratio = readScalar(bndLayers.lookup("thicknessRatio"));
            refLayers.setGlobalThicknessRatio(ratio);
        }
        if( bndLayers.found("maxFirstLayerThickness") )
        {
            const scalar maxFirstThickness =
                readScalar(bndLayers.lookup("maxFirstLayerThickness"));
            refLayers.setGlobalMaxThicknessOfFirstLayer(maxFirstThickness);
        }

        //- patch-based properties
        if( bndLayers.isDict("patchBoundaryLayers") )
        {
            const dictionary& patchBndLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList patchNames = patchBndLayers.toc();

            forAll(patchNames, patchI)
            {
                const word pName = patchNames[patchI];

                if( patchBndLayers.isDict(pName) )
                {
                    const dictionary& patchDict =
                        patchBndLayers.subDict(pName);

                    if( patchDict.found("nLayers") )
                    {
                        const label nLayers =
                            readLabel(patchDict.lookup("nLayers"));
                        refLayers.setNumberOfLayersForPatch(pName, nLayers);
                    }
                    if( patchDict.found("thicknessRatio") )
                    {
                        const scalar ratio =
                            readScalar(patchDict.lookup("thicknessRatio"));
                        refLayers.setThicknessRatioForPatch(pName, ratio);
                    }
                    if( patchDict.found("maxFirstLayerThickness") )
                    {
                        const scalar maxFirstThickness =
                            readScalar
                            (
                                patchDict.lookup("maxFirstLayerThickness")
                            );
                        refLayers.setMaxThicknessOfFirstLayerForPatch
                        (
                            pName,
                            maxFirstThickness
                        );
                    }
                    if( patchDict.found("allowDiscontinuity") )
                    {
                        const bool allowDiscontinuity =
                            readBool(patchDict.lookup("allowDiscontinuity"));

                        if( allowDiscontinuity )
                            refLayers.setInteruptForPatch(pName);
                    }
                }
                else
                {
                    Warning << "Cannot refine layer for patch "
                        << patchNames[patchI] << endl;
                }
            }
        }
    }

    refLayers.refineLayers();

    Info << "Finished with bnd layer refinement "
         << runTime.elapsedClockTime() << endl;

    //polyMeshGenChecks::checkMesh(pmg, true);

    pmg.write();
    //meshOctree* octreePtr = NULL;
    //meshOptimizer(*octreePtr, pmg).preOptimize();

    //writeMeshEnsight(pmg, "meshWithBndLayers");
    //pmg.addressingData().checkMesh(true);

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
