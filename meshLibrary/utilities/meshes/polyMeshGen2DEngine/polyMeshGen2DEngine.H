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

Class
    polyMeshGen2DEngine

Description
    A simple engine which provides topological information of a 2D mesh
    and allows for maintaining consistency

SourceFiles
    polyMeshGen2DEngine.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGen2DEngine_H
#define polyMeshGen2DEngine_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen.H"
#include "boolList.H"
#include "labelList.H"
#include "boundBox.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGen2DEngine
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- bounding box
        boundBox bb_;

        //- which faces are in the x-y plane and at smaller z
        boolList activeFace_;

        //- labels of active faces
        labelList activeFaces_;

        //- which points are in the x-y plane and at smaller z
        boolList activePoint_;

        //- labels of active points
        labelList activePoints_;

        //- labels of offset points of the active points
        labelList offsetPoints_;

    // Private member functions
        //- find active faces
        void findActiveFacesAndPoints();

        //- disallow bitwise assignment
        void operator=(const polyMeshGen2DEngine&);

        //- disallow copy construct
        polyMeshGen2DEngine(const polyMeshGen2DEngine&);

public:

    // Constructors
        //- Construct from the reference to the mesh
        polyMeshGen2DEngine(polyMeshGen& mesh);

    // Destructor
        ~polyMeshGen2DEngine();

    // Member functions
        //- const access to active faces
        inline const boolList& activeFace() const
        {
            return activeFace_;
        }

        //- labels of active faces
        inline const labelList& activeFaceLabels() const
        {
            return activeFaces_;
        }

        //- const access to active points
        inline const boolList& activePoint() const
        {
            return activePoint_;
        }

        //- labels of active points
        inline const labelList& activePointLabels() const
        {
            return activePoints_;
        }

        //- labels of offset points
        inline const labelList& offsetPointLabels() const
        {
            return offsetPoints_;
        }

        //- correct x and y coordinates of offset points
        //- and unify the z coordinate
        void correctPoints();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
