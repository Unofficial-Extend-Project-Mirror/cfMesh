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
    Reads the specified surface and writes it in the fms format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "triSurfaceChecks.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);

    triSurf surf(inFileName);

    label nFailed(0);

    boundBox bb;
    triSurfaceChecks::calculateBoundingBox(surf, bb);
    Info << "Surface bounding box is " << bb << endl;

    //- calculate manifolds
    const label nManifolds = triSurfaceChecks::checkSurfaceManifolds(surf);
    if( nManifolds > 1 )
    {
        ++nFailed;

        Info << "Surface mesh consists of " << nManifolds
             << " manifolds." << endl;
        Warning << "You cannot mesh geometries consisting of more than"
                << " one domain, and it must not contain baffles." << endl;
    }
    else
    {
        Info << "Surface mesh consists of a single manifold." << endl;
    }

    //- find open boundary edges
    if( triSurfaceChecks::checkForHoles(surf) )
    {
        ++nFailed;

        Info << "Surface mesh has open boundaries!!" << endl;
        Warning << "This indicates that there may be some holes in the surface"
                << " mesh. Holes in the mesh must be smaller than the specified"
                << " cell size at this location. In addition, please avoid"
                << " using automatic refinement (minCellSize)." << endl;
    }
    else
    {
        Info << "No holes found in the surface mesh." << endl;
    }

    //- find non-manifold edges
    if( triSurfaceChecks::checkForNonManifoldEdges(surf) )
    {
        ++nFailed;

        Info << "Surface mesh has non-manifold edges!!" << endl;
        Warning << "This indicates that the surface mesh consists of multiple"
                << " domains and/or baffles. Please make sure that they are not"
                << " in the domain which shall be meshed." << endl;
    }
    else
    {
        Info << "Surface does not have any non-manifold edges." << endl;
    }

    //- check the number of disconnected parts
    if( triSurfaceChecks::checkDisconnectedParts(surf) > 1 )
    {
        ++nFailed;

        Info << "Surface mesh consists of disconnected parts." << endl;
        Warning << "This is not a problem if there exists a region surrounding"
                << " the other ones! In other case, the mesher will generate"
                << " the mesh in the domains with most cells." << endl;
    }
    else
    {
        Info << "Surface mesh consists of a single region." << endl;
    }

    //- find triangles with small angles
    if( triSurfaceChecks::checkAngles(surf, "smallAngles", 1.0) )
    {
        ++nFailed;

        Info << "Surface mesh has some bad-quality triangles." << endl;
        Warning << "This may cause problems to the automatic refinement"
                << " procedure (minCellSize). " << endl;
    }
    else
    {
        Info << "No sliver triangles found." << endl;
    }

    if( nFailed )
    {
        Warning << "Found " << nFailed
                << " checks indicating potential problems." << endl;
        Warning << "This does not mean that you cannot generate"
                << " a valid mesh. " << endl;

        surf.writeSurface(inFileName);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
