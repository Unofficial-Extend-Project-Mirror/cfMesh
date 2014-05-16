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

#include "boolList.H"
#include "triSurface.H"
#include "surfaceIntersectionsOctreeCube.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void surfaceIntersectionsOctreeCube::refineTree(const short maxNEl, const direction maxL)
{
    if( (containedElements_.size() > maxNEl) && (level_ < maxL) )
    {
        # ifdef DEBUGSearch
        Info << "Cube " << cubeBox_ << " is in level " << level_
            << " Has " << containedElements_.size() << " elements" << endl;
        # endif
        //- create subCubes_
        subCubesPtr_ = new FixedList<surfaceIntersectionsOctreeCube*, 8>();
        FixedList<surfaceIntersectionsOctreeCube*, 8>& subCubes_ = *subCubesPtr_;
        const point c = (cubeBox_.min() + cubeBox_.max()) / 2.0;
        {
            boundBox bb(cubeBox_.min(), c);
            subCubes_[0] = new surfaceIntersectionsOctreeCube(surface_, bb, level_+1);
        }
        {
            const point min =
                point(c.x(), cubeBox_.min().y(), cubeBox_.min().z());
            const point max =
                point(cubeBox_.max().x(), c.y(), c.z());
            subCubes_[1] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }
        {
            const point min = 
                point(c.x(), c.y(), cubeBox_.min().z());
            const point max =
                point(cubeBox_.max().x(), cubeBox_.max().y(), c.z());
            subCubes_[2] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }
        {
            const point min =
                point(cubeBox_.min().x(), c.y(), cubeBox_.min().z());
            const point max =
                point(c.x(), cubeBox_.max().y(), c.z());
            subCubes_[3] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }
        {
            const point min =
                point(cubeBox_.min().x(), cubeBox_.min().y(), c.z());
            const point max =
                point(c.x(), c.y(), cubeBox_.max().z());
            subCubes_[4] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }
        {
            const point min =
                point(c.x(), cubeBox_.min().y(), c.z());
            const point max = 
                point(cubeBox_.max().x(), c.y(), cubeBox_.max().z());
            subCubes_[5] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }
        {
            subCubes_[6] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(c, cubeBox_.max()), level_+1);
        }
        {
            const point min =
                point(cubeBox_.min().x(), c.y(), c.z());
            const point max =
                point(c.x(), cubeBox_.max().y(), cubeBox_.max().z());
            subCubes_[7] =
                new surfaceIntersectionsOctreeCube(surface_, boundBox(min, max), level_+1);
        }

        # ifdef DEBUGSearch
        forAll(subCubes_, scI)
            Info << "Refined cube " << scI << " is "
                << subCubes_[scI]->bb() << endl;
        # endif
        //- check if the subCube contain the element
        //- if it does store it into the cubCube
        for(SLList<label>::const_iterator eIter = containedElements_.begin();
            eIter != containedElements_.end();
            ++eIter
        )
        {
            forAll(subCubes_, cubeI)
                if( subCubes_[cubeI]->intersectsTriangle(eIter()) )
                {
                    subCubes_[cubeI]->append(eIter());
                }
        }

        //- clear contained elements because they are store in the subCubes
        containedElements_.clear();

        //- continue with decomposition
        forAll(subCubes_, cubeI)
            subCubes_[cubeI]->refineTree(maxNEl, maxL);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
