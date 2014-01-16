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
    Reads the AVL's surface mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "objectRegistry.H"
#include "triSurf.H"
#include "triFaceList.H"
#include "labelLongList.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if( inFileName.ext() != "flma" )
    {
        Info << "Cannot convert this mesh" << endl;
        return 0;
    }

    label counter;

    IFstream inFile(inFileName);

    inFile >> counter;

    //- read vertices
    pointField points(counter);
    for(label pointI=0;pointI<counter;++pointI)
    {
        point p;
        inFile >> p.x();
        inFile >> p.y();
        inFile >> p.z();

        points[pointI] = p;
    }

    //- read facets
    inFile >> counter;
    geometricSurfacePatchList patches(1);
    patches[0].name() = "patch";
    LongList<labelledTri> triangles(counter);
    forAll(triangles, triI)
    {
        inFile >> counter;

        if( counter != 3 )
        {
            Info << "Facet " << triI << " is not a triangle!!" << endl;
            Warning << "Cannot convert this surface!" << endl;
            return 0;
        }

        for(label j=0;j<3;++j)
            inFile >> triangles[triI][2-j];

        triangles[triI].region() = 0;
    }

    //- read cell types
    inFile >> counter;
    forAll(triangles, triI)
        inFile >> counter;

    //- create the surface mesh
    triSurf ts(triangles, patches, points);

    //- start reading selections
    inFile >> counter;
    for(label selI=0;selI<counter;++selI)
    {
        //- read selection name
        word selName;
        inFile >> selName;

        //- read selection type
        label selType;
        inFile >> selType;

        //- read selection entries
        label size;
        inFile >> size;
        labelLongList entries(size);
        for(label i=0;i<size;++i)
            inFile >> entries[i];

        //- store cell selections
        if( selType == 2 )
        {
            Info << "Adding subset " << selName << endl;
            const label setID = ts.addFacetSubset(selName);

            forAll(entries, i)
                ts.addFacetToSubset(setID, entries[i]);
        }
    }

    //- write the surface

    ts.writeSurface(outFileName);

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
