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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "tetTessellation.H"
#include "triSurface.H"
#include "DynList.H"

//#define DEBUGTessalation

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
void tetTessellation::createInitialTets()
{
    boundaryPoints_[0] = points_.size();
    points_.append(min_);
    boundaryPoints_[1] = points_.size();
    points_.append(point(max_.x(), min_.y(), min_.z()));
    boundaryPoints_[2] = points_.size();
    points_.append(point(max_.x(), max_.y(), min_.z()));
    boundaryPoints_[3] = points_.size();
    points_.append(point(min_.x(), max_.y(), min_.z()));
    boundaryPoints_[4] = points_.size();
    points_.append(point(min_.x(), min_.y(), max_.z()));
    boundaryPoints_[5] = points_.size();
    points_.append(point(max_.x(), min_.y(), max_.z()));
    boundaryPoints_[6] = points_.size();
    points_.append(max_);
    boundaryPoints_[7] = points_.size();
    points_.append(point(min_.x(), max_.y(), max_.z()));
    
    DynList<label> tetLabels(5);
    tetLabels.append(elmts_.size());
    elmts_.append
    (
        tessellationElement
        (
            boundaryPoints_[0],
            boundaryPoints_[4],
            boundaryPoints_[1],
            boundaryPoints_[3]
        )
    );
    tetLabels.append(elmts_.size());
    elmts_.append
    (
        tessellationElement
        (
            boundaryPoints_[5],
            boundaryPoints_[1],
            boundaryPoints_[4],
            boundaryPoints_[6]
        )
    );
    tetLabels.append(elmts_.size());
    elmts_.append
    (
        tessellationElement
        (
            boundaryPoints_[7],
            boundaryPoints_[4],
            boundaryPoints_[3],
            boundaryPoints_[6]
        )
    );
    tetLabels.append(elmts_.size());
    elmts_.append
    (
        tessellationElement
        (
            boundaryPoints_[2],
            boundaryPoints_[3],
            boundaryPoints_[1],
            boundaryPoints_[6]
        )
    );
    tetLabels.append(elmts_.size());
    elmts_.append
    (
        tessellationElement
        (
            boundaryPoints_[3],
            boundaryPoints_[4],
            boundaryPoints_[1],
            boundaryPoints_[6]
        )
    );
    
    elmts_[tetLabels[0]].setNeighbour(0, tetLabels[4]);
    elmts_[tetLabels[1]].setNeighbour(0, tetLabels[4]);
    elmts_[tetLabels[2]].setNeighbour(0, tetLabels[4]);
    elmts_[tetLabels[3]].setNeighbour(0, tetLabels[4]);
    
    elmts_[tetLabels[4]].setNeighbour(0, tetLabels[1]);
    elmts_[tetLabels[4]].setNeighbour(1, tetLabels[3]);
    elmts_[tetLabels[4]].setNeighbour(2, tetLabels[2]);
    elmts_[tetLabels[4]].setNeighbour(3, tetLabels[0]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct from list of points
tetTessellation::tetTessellation
(
    const LongList<point>& points
)
:
    points_(points),
    elmts_(),
    pointOk_(true),
    nElmts_(),
    max_(),
    min_()
{
    point c;
    if( points.size() > 1)
    {
        max_ = point(-GREAT, -GREAT, -GREAT);
        min_ = point(GREAT, GREAT, GREAT);

        forAll(points, pI)
        {
            const point& p = points[pI];

            if( p.x() > max_.x() ) max_.x() = p.x();
            if( p.y() > max_.y() ) max_.y() = p.y();
            if( p.z() > max_.z() ) max_.z() = p.z();

            if( p.x() < min_.x() ) min_.x() = p.x();
            if( p.y() < min_.y() ) min_.y() = p.y();
            if( p.z() < min_.z() ) min_.z() = p.z();
        }
        c = (max_ + min_) / 2.0;
    }
    else
    {
        c = max_ = min_ = points_[0];
    }

    createInitialTets();

    # ifdef DEBUGTessalation
    Info << "Elements " << elmts_ << endl;
    forAll(elmts_, eI)
        Info << "Volume of element " << eI
            << " is " << elmts_[eI].mag(points_) << endl;
    Info << "Number of points is " << points_.size() << endl;
    Info << "Number of tetrahedra is  " << elmts_.size() << endl;
    # endif
}

// Construct from triSurface (only the first tet is created)
tetTessellation::tetTessellation
(
    const triSurface& surf
)
:
    points_(),
    elmts_(),
    pointOk_(true),
    nElmts_(),
    max_(Foam::max(surf.localPoints())),
    min_(Foam::min(surf.localPoints()))
{
    createInitialTets();

    # ifdef DEBUGTessalation
    Info << "Elements " << elmts_ << endl;
    forAll(elmts_, eI)
        Info << "Volume of element " << eI
            << " is " << elmts_[eI].mag(points_) << endl; 
    # endif
}

// Construct from boundBox (only the first tet is created)
tetTessellation::tetTessellation(const boundBox& bb)
:
    points_(),
    elmts_(),
    pointOk_(true),
    nElmts_(),
    max_(bb.max()),
    min_(bb.min())
{
    createInitialTets();

    # ifdef DEBUGTessalation
    Info << "Elements " << elmts_ << endl;
    forAll(elmts_, eI)
        Info << "Volume of element " << eI
            << " is " << elmts_[eI].mag(points_) << endl; 
    # endif
}

tetTessellation::~tetTessellation()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
