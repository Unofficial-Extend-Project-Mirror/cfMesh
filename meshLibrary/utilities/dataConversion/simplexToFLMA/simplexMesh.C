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

#include "simplexMesh.H"
#include "IOmanip.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::simplexMesh::writePoints(Foam::OFstream& flmaGeometryFile) const
{
    flmaGeometryFile << mesh_.pts().size() << nl;
    const DynList<point, 128>& points = mesh_.pts();
    forAll(points, pointI)
    {
        const point& p = points[pointI];
        flmaGeometryFile << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';
    }
    
    flmaGeometryFile << nl;
}

void simplexMesh::writeCells(OFstream& flmaGeometryFile) const
{
    const DynList<partTet, 128>& tets = mesh_.tets();
    
    flmaGeometryFile << tets.size() << nl;
    
    forAll(tets, tetI)
    {
        const partTet& tet = tets[tetI];
        flmaGeometryFile << 4;
        for(label i=0;i<4;++i)
            flmaGeometryFile << ' ' << tet[i];
        flmaGeometryFile << nl;
    }
}

void Foam::simplexMesh::writeCellTypes(OFstream& flmaGeometryFile) const
{
    const DynList<partTet, 128>& tets = mesh_.tets();
    flmaGeometryFile << nl << tets.size() << nl;
    forAll(tets, tetI)
        flmaGeometryFile << 4 << nl;
    flmaGeometryFile << nl;
}

void Foam::simplexMesh::writeSelections(Foam::OFstream& flmaGeometryFile) const
{
    flmaGeometryFile << 0;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from partTetMeshSimplex
Foam::simplexMesh::simplexMesh(const partTetMeshSimplex& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simplexMesh::~simplexMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simplexMesh::write
(
    OFstream& flmaGeometryFile
) const
{
    writePoints(flmaGeometryFile);

    writeCells(flmaGeometryFile);
    
    writeCellTypes(flmaGeometryFile);

    writeSelections(flmaGeometryFile);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
