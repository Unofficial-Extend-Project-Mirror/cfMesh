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
#include "labelListPMG.H"
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
    argList::validArgs.append("subset file name");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);
    fileName subsetFileName(args.args()[3]);
    
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
    triFaceList triangles(counter);
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
    }
    
    triSurface ts(triangles, points);
    ts.write(outFileName);
    
    //- read cell types
    inFile >> counter;
    forAll(triangles, triI)
        inFile >> counter;
    
    //- start reading selections
    std::map<word, labelListPMG> subsets;
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
        labelListPMG entries(size);
        for(label i=0;i<size;++i)
            inFile >> entries[i];
        
        //- store cell selections
        if( selType == 2 )
        {
            Info << "Adding subset " << selName << endl;
            subsets.insert(std::pair<word, labelListPMG>(selName, entries));
        }
    }
    
    OFstream outFile(subsetFileName);
    
    std::map<word, labelListPMG>::const_iterator iter;
    counter = 0;
    for(iter=subsets.begin();iter!=subsets.end();++iter)
        ++counter;
    
    outFile << counter << nl;
    for(iter=subsets.begin();iter!=subsets.end();++iter)
    {
        outFile << iter->first << nl;
        outFile << iter->second << nl;
    }
    
/*    for(iter=subsets.begin();iter!=subsets.end();++iter)
    {
        labelList faceMap, pointMap;
        boolList useTri(triangles.size(), false);
        
        const word& sName = iter->first;
        const labelListPMG& elmts = iter->second;
        
        forAll(elmts, elI)
            useTri[elmts[elI]] = true;
        
        faceMap.setSize(elmts.size());
        pointMap.setSize(elmts.size());
        
        triSurface subSurface = ts.subsetMesh(useTri, pointMap, faceMap);
        
        subSurface.write(sName+".stl");
    }
*/    
    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
