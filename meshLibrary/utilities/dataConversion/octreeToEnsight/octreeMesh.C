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

#include "octreeMesh.H"
#include "cellModeller.H"
#include "IOmanip.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from meshOctree
Foam::octreeMesh::octreeMesh(const meshOctree& octree)
:
    points_(3*octree.numberOfLeaves()),
    boxes_(octree.numberOfLeaves()),
    useBox_(octree.numberOfLeaves(), true),
    octree_(octree)
{
    createPointsAndCells();
}

// Construct from meshOctree and cube type
Foam::octreeMesh::octreeMesh
(
    const meshOctree& octree,
    const direction cubeType
)
:
    points_(3*octree.numberOfLeaves()),
    boxes_(octree.numberOfLeaves()),
    useBox_(octree.numberOfLeaves(), false),
    octree_(octree)
{
    label i, nLeaves = octree.numberOfLeaves();
    for(i=0;i<nLeaves;++i)
        if( octree_.returnLeaf(i).cubeType() & cubeType )
            useBox_[i] = true;
    
    createPointsAndCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeMesh::~octreeMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void octreeMesh::createPointsAndCells()
{
    label boxI, nLeaves = octree_.numberOfLeaves();
    
    const cellModel& hex = *cellModeller::lookup("hex");
    const boundBox& rootBox = octree_.rootBox();
    Info << "Root box for writting octree " << rootBox << endl;
    
    label pointI(0), boxJ(0);
    
    for(boxI=0;boxI<nLeaves;++boxI)
        if( useBox_[boxI] )
        {
            labelList cv(8);
            
            const meshOctreeCubeBasic& c = octree_.returnLeaf(boxI);
            FixedList<point, 8> vrt;
            c.vertices(rootBox, vrt);
            
            forAll(vrt, vI)
            {
                points_.newElmt(pointI) = vrt[vI];
                switch( vI )
                {
                    case 0:
                    {
                        cv[4] = pointI;
                        break;
                    }
                    case 1:
                    {
                        cv[5] = pointI;
                        break;
                    }
                    case 2:
                    {
                        cv[0] = pointI;
                        break;
                    }
                    case 3:
                    {
                        cv[1] = pointI;
                        break;
                    }
                    case 4:
                    {
                        cv[7] = pointI;
                        break;
                    }
                    case 5:
                    {
                        cv[6] = pointI;
                        break;
                    }
                    case 6:
                    {
                        cv[3] = pointI;
                        break;
                    }
                    case 7:
                    {
                        cv[2] = pointI;
                        break;
                    }
                };
                
                ++pointI;
            }
            
            boxes_[boxJ++] = cellShape(hex, cv);
        }
    
    points_.setSize(pointI);
    boxes_.setSize(boxJ);
}

void Foam::octreeMesh::writePoints
(
    const Foam::scalarField& pointsComponent,
    Foam::OFstream& ensightGeometryFile
) const
{
    forAll(pointsComponent, pointi)
    {
        ensightGeometryFile << setw(12) << pointsComponent[pointi] << nl;
    }
}

void Foam::octreeMesh::writeBoxes
(
    OFstream& ensightGeometryFile
) const
{
    ensightGeometryFile << "hexa8" << nl << setw(10) << boxes_.size() << nl;

    forAll(boxes_, boxI)
    {
        const cellShape& box = boxes_[boxI];
        
        forAll(box, pI)
            ensightGeometryFile << setw(10) << box[pI] + 1;
        
        ensightGeometryFile << nl;
    }
}


void octreeMesh::write
(
    OFstream& ensightGeometryFile
) const
{
    // Set Format
    ensightGeometryFile.setf
    (
        ios_base::scientific,
        ios_base::floatfield
    );
    ensightGeometryFile.precision(5);

    ensightGeometryFile
        << "FOAM Geometry File " << nl << nl
        << "node id assign" << nl
        << "element id assign" << nl;
    
    ensightGeometryFile
        << "part" << nl
        << setw(10) << 1 << nl
        << "FOAM cells" << nl
        << "coordinates" << nl
        << setw(10) << points_.size()
        << endl;

    for (direction d=0; d<vector::nComponents; d++)
    {
        writePoints(points_.component(d)(), ensightGeometryFile);
    }

    writeBoxes
    (
        ensightGeometryFile
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
