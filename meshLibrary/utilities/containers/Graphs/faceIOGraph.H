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

Class
    faceIOGraph

Description
    A graph of faces which support automated output

SourceFiles
    faceIOGraph.C

\*---------------------------------------------------------------------------*/

#ifndef faceIOGraph_H
#define faceIOGraph_H

#include "VRWGraph.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceIOGraph Declaration
\*---------------------------------------------------------------------------*/

class faceIOGraph
:
    public regIOobject,
    public VRWGraph
{

public:

    //- Runtime type information
		TypeName("faceIOGraph");

    // Constructors

        //- Construct from IOobject
        faceIOGraph(const IOobject&);

        //- Construct from IOobject and size of faceIOGraph
        faceIOGraph(const IOobject&, const label);

        //- Construct from IOobject and a VRWGraph
        faceIOGraph(const IOobject&, const VRWGraph&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

		void operator=(const faceIOGraph&);
        void operator=(const VRWGraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
