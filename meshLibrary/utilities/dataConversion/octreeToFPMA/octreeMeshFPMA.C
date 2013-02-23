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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "octreeMeshFPMA.H"
#include "IOmanip.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from meshOctree
Foam::octreeMeshFPMA::octreeMeshFPMA
(
    const meshOctree& octree,
    const dictionary& dict
)
:
	addressing_(octree, dict, false)
{
    addressing_.boxType();
    for(label leafI=0;leafI<octree.numberOfLeaves();++leafI)
        addressing_.setBoxType(leafI, meshOctreeAddressing::MESHCELL);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeMeshFPMA::~octreeMeshFPMA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::octreeMeshFPMA::writePoints
(
    Foam::OFstream& fpmaGeometryFile
) const
{
    const pointField& points = addressing_.octreePoints();
	fpmaGeometryFile << points.size() << nl;
	
	forAll(points, pointI)
	{
		const point& p = points[pointI];
		fpmaGeometryFile << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';
	}
	
	fpmaGeometryFile << nl;
}

void octreeMeshFPMA::writeCells
(
    OFstream& fpmaGeometryFile
) const
{
	const VRWGraph& cells = addressing_.leafFaces();
	
	fpmaGeometryFile << cells.size() << nl;
	forAll(cells, cellI)
	{
		fpmaGeometryFile << cells.sizeOfRow(cellI);
		forAllRow(cells, cellI, fI)
			fpmaGeometryFile << ' ' << cells(cellI, fI);
		fpmaGeometryFile << nl;
	}
}

void Foam::octreeMeshFPMA::writeFaces
(
    OFstream& fpmaGeometryFile
) const
{
	const VRWGraph& faces = addressing_.octreeFaces();
	fpmaGeometryFile << faces.size() << nl;
	forAll(faces, faceI)
	{
		const constRow f = faces[faceI];
		
		fpmaGeometryFile << f.size();
		forAllReverse(f, pI)
			fpmaGeometryFile << ' ' << f[pI];
		fpmaGeometryFile << nl;
	}
}

void Foam::octreeMeshFPMA::writeSelections
(
	Foam::OFstream& fpmaGeometryFile
) const
{
	fpmaGeometryFile << 4 << nl;
    //- create outside selection
    fpmaGeometryFile << "Outside" << nl << 2 << nl;
    labelListPMG selBoxes;
    const meshOctree& octree = addressing_.octree();
    for(label i=0;i<octree.numberOfLeaves();++i)
    {
        if( octree.returnLeaf(i).cubeType() & meshOctreeCube::OUTSIDE )
            selBoxes.append(i);
    }
    fpmaGeometryFile << selBoxes.size() << nl;
    forAll(selBoxes, boxI)
        fpmaGeometryFile << selBoxes[boxI] << ' ';
    fpmaGeometryFile << nl;
    
    //- create inside selection
    fpmaGeometryFile << "Inside" << nl << 2 << nl;
    selBoxes.clear();
    for(label i=0;i<octree.numberOfLeaves();++i)
    {
        if( octree.hasContainedTriangles(i) )
            continue;
        if( octree.returnLeaf(i).cubeType() & meshOctreeCube::INSIDE )
            selBoxes.append(i);
    }
    fpmaGeometryFile << selBoxes.size() << nl;
    forAll(selBoxes, boxI)
        fpmaGeometryFile << selBoxes[boxI] << ' ';
    fpmaGeometryFile << nl;
    
    //- create unknown selection
    fpmaGeometryFile << "Unknown" << nl << 2 << nl;
    selBoxes.clear();
    for(label i=0;i<octree.numberOfLeaves();++i)
    {
        if( octree.returnLeaf(i).cubeType() & meshOctreeCube::UNKNOWN )
            selBoxes.append(i);
    }
    fpmaGeometryFile << selBoxes.size() << nl;
    forAll(selBoxes, boxI)
        fpmaGeometryFile << selBoxes[boxI] << ' ';
    fpmaGeometryFile << nl;
    
    //- create DATA selection
    fpmaGeometryFile << "Data" << nl << 2 << nl;
    selBoxes.clear();
    for(label i=0;i<octree.numberOfLeaves();++i)
    {
        if( octree.hasContainedTriangles(i) )
            selBoxes.append(i);
    }
    fpmaGeometryFile << selBoxes.size() << nl;
    forAll(selBoxes, boxI)
        fpmaGeometryFile << selBoxes[boxI] << ' ';
    fpmaGeometryFile << nl;
}


void octreeMeshFPMA::write
(
    OFstream& fpmaGeometryFile
) const
{
	writePoints(fpmaGeometryFile);
	
	writeFaces(fpmaGeometryFile);

	writeCells(fpmaGeometryFile);

	writeSelections(fpmaGeometryFile);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
