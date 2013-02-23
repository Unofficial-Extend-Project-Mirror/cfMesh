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
    Test of surface morpher

Description
    - morphs mesh surface to make it smooth for mapping

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyMeshGen.H"
#include "surfaceMorpherCells.H"
#include "topologicalCleaner.H"
#include "writeMeshEnsight.H"
#include "polyMeshGenAddressing.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
	
	polyMeshGen pmg(runTime);
	pmg.read();
	
    //- test number of nodes per processor face
    const pointFieldPMG& points = pmg.points();
    const faceListPMG& faces = pmg.faces();
    const PtrList<writeProcessorPatch>& procBoundaries =
        pmg.procBoundaries();
    
    //- check face count
    const label nIntFaces = pmg.nInternalFaces();
    label nBndFaces(0);
    const PtrList<writePatch>& boundaries = pmg.boundaries();
    forAll(boundaries, patchI)
        nBndFaces += boundaries[patchI].patchSize();
    label nProcFaces(0);
    forAll(procBoundaries, patchI)
        nProcFaces += procBoundaries[patchI].patchSize();

    if( (nIntFaces+nBndFaces+nProcFaces) != faces.size() )
        FatalError << "Number of faces is not correct!" << abort(FatalError);
    
    //- check if the centres of the proc patch match
    forAll(procBoundaries, patchI)
    {
        vectorField cents(procBoundaries[patchI].patchSize());
        
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + cents.size();
        for(label faceI=start;faceI<end;++faceI)
            cents[faceI-start] = faces[faceI].centre(points);
        
        OPstream toOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        toOtherProc << cents;
    }
    
    forAll(procBoundaries, patchI)
    {
        vectorField otherCentres;
        IPstream fromOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        fromOtherProc >> otherCentres;
        
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        for(label faceI=start;faceI<end;++faceI)
            if(
                mag(otherCentres[faceI-start] - faces[faceI].centre(points))
                > 1e-6 * mag(otherCentres[faceI-start])
            )
            {
                Serr << Pstream::myProcNo() << "This centre "
                    << faces[faceI].centre(points) << " other centre "
                    << otherCentres[faceI-start] << endl;
                Serr << "Face centres do not match!" << endl;
            }
    }
    
    //- check number of vertices in proc-boundary faces
    forAll(procBoundaries, patchI)
    {
        labelList nVrtPerFace(procBoundaries[patchI].patchSize());
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + nVrtPerFace.size();
        for(label faceI=start;faceI<end;++faceI)
            nVrtPerFace[faceI-start] = faces[faceI].size();
        
        OPstream toOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        toOtherProc << nVrtPerFace;
    }
    
    forAll(procBoundaries, patchI)
    {
        labelList nVrtPerFace;
        IPstream fromOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        fromOtherProc >> nVrtPerFace;
        
        label nTotalNei(0);
        forAll(nVrtPerFace, fI)
            nTotalNei += nVrtPerFace[fI];
        
        label nTotalOwn(0);
        
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + nVrtPerFace.size();
        for(label faceI=start;faceI<end;++faceI)
            nTotalOwn += faces[faceI].size();
        
        if( nTotalNei != nTotalOwn )
            Serr << "Patches do not contain same number of points!" << endl;
        
        for(label faceI=start;faceI<end;++faceI)
            if( nVrtPerFace[faceI-start] != faces[faceI].size() )
            {
                Serr << Pstream::myProcNo() << "Face " << faces[faceI].size()
                    << " other face size " << nVrtPerFace[faceI-start]
                    << endl;
                FatalErrorIn
                (
                    "void surfaceMorpherCells::morphMesh()"
                ) << "Mesh is not correct!" << Pstream::myProcNo()
                    << abort(FatalError);
            }
    }
    
    //- test if vertices of processor boundaries match
    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        label nPointsToSend(0);
        for(label faceI=start;faceI<end;++faceI)
            nPointsToSend += faces[faceI].size();
        
        pointField pointsToSend(nPointsToSend);
        nPointsToSend = 0;
        for(label faceI=start;faceI<end;++faceI)
        {
            const face f = faces[faceI].reverseFace();
            forAll(f, pI)
                pointsToSend[nPointsToSend++] = points[f[pI]];
        }
        
        OPstream toOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        toOtherProc << pointsToSend;
    }
    
    forAll(procBoundaries, patchI)
    {
        pointField pointsToReceive;
        
        IPstream fromOtherProc
        (
            Pstream::nonBlocking,
            procBoundaries[patchI].neiProcNo()
        );
        fromOtherProc >> pointsToReceive;
        
        label counter(0);
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        for(label faceI=start;faceI<end;++faceI)
        {
            const scalar d = Foam::sqrt(mag(faces[faceI].normal(points)));
            const face& f = faces[faceI];
            forAll(f, pI)
                if( mag(pointsToReceive[counter++] - points[f[pI]]) > 0.01*d )
                {
                    Serr << "Received point " << pointsToReceive[counter-1]
                        << "Searched point " << points[f[pI]] << endl;
                    FatalError << "Face " << faceI << " of processor patch "
                        << patchI << " is not correct!" << abort(FatalError);
                }
        }
	}
    
/*    Serr << Pstream::myProcNo() << "Global point labels "
        << pmg.addressingData().globalPointLabel() << endl;
    Serr << Pstream::myProcNo() << "Point procs "
        << pmg.addressingData().pointAtProcs() << endl;
*/
/*    Serr << Pstream::myProcNo() << "Global face labels "
        << pmg.addressingData().globalFaceLabel() << endl;
    Serr << Pstream::myProcNo() << "Global cell labels "
        << pmg.addressingData().globalCellLabel() << endl;
*/
/*    Serr << Pstream::myProcNo() << "Global edge labels "
        << pmg.addressingData().globalEdgeLabel() << endl;
    Serr << Pstream::myProcNo() << "Edge procs "
        << pmg.addressingData().edgeAtProcs() << endl;
*/
    pmg.write();
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
