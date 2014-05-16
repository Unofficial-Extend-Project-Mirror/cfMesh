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

Class
    createFacesFromChain

Description
    Creates boundary faces from a chain of open edges

SourceFiles
    createFacesFromChain.C

\*---------------------------------------------------------------------------*/

#ifndef createFacesFromChain_H
#define createFacesFromChain_H

#include "DynList.H"
#include "labelList.H"
#include "face.H"
#include "boolList.H"

namespace Foam
{

class VRWGraph;

/*---------------------------------------------------------------------------*\
                Class createFacesFromChain Declaration
\*---------------------------------------------------------------------------*/

class createFacesFromChain
{
    // Members
        labelList chainPoints_;

        List<DynList<label> > regionsForPointAtPosition_;

        DynList<face> createdFaces_;
        DynList<label> faceRegions_;

    // Private member functions
        void findPointsBelongingToTheFace
        (
            const label currPos,
            boolList& includePoints,
            boolList& endPoints
        ) const;

        void shrinkTheChain
        (
            const label currPos,
            const boolList& includePoints,
            const boolList& endPoints
        );

    public:

    // Constructors
        //- construct from chain points and a graph of regions
        //- to which each point belongs to
        createFacesFromChain
        (
            const labelList& chPoints,
            const VRWGraph& pointRegions
        );

        //- construct from chain points and regions
        //- to which each point belongs to
        createFacesFromChain
        (
            const labelList& chPoints,
            const List<DynList<label> >& pointRegions
        );

    // Destructor

        ~createFacesFromChain();

    // Member functions
        //- create faces which do not require an additional corner
        void createFacesWithoutACorner();

        //- create faces with a corner
        void createFacesWithACorner(const label);

        //- return created faces and the corresponding patches
        const DynList<face>& createdFaces() const;
        const DynList<label>& faceRegion() const;

        //- a list of points which have not yet been resolved
        const labelList& unresolvedPoints() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
