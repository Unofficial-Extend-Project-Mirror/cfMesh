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
    boundaryFacesGenerator

Description
    Generates boundary faces from open edges of boundary cells

SourceFiles
    boundaryFacesGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryFacesGenerator_H
#define boundaryFacesGenerator_H

#include "polyMeshGenModifier.H"
#include "triSurface.H"
#include "DynList.H"
#include "labelListPMG.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class boundaryFacesGenerator Declaration
\*---------------------------------------------------------------------------*/

class boundaryFacesGenerator
{
    // private data
        //- triangulated surface
        const triSurface& surface_;

        //- mesh data
        polyMeshGen& mesh_;
    
        //- helpers
        List<direction>& nFacesInCell_;
        const boolList& boundaryCell_;
        VRWGraph& pointRegions_;
    
        label& nPoints_;
        const label nIntFaces_;

        //- boundary faces
        VRWGraph boundaryFaces_;
        labelListPMG boundaryOwners_;
        labelListPMG boundaryPatches_;
        bool hasFacesInDefaultPatch_;
        label defaultPatchID_;
    
        //- surface corners
        DynList<label> surfaceCorners_;
        DynList< DynList<label> > cornersPatches_;

    // private member functions
        //- find patches for the given vertex
        labelList patchesForPoint(const label) const;

        //- find patches for a given cell
        void determineCellPatches(const label, DynList<label>&);
        
        //- find boundary chains for the given cell
        void createBoundaryChainsForCell(const label, labelListList&);
        
        void createFacesForChain
        (
            const labelList& chainVertices,
            const DynList<label>& patches,
            List<DynList<face> >& facesForChain
        );
        void createFacesForChainWithCorner
        (
            const label cornerLabel,
            const labelList& bChain,
            const DynList<label>& patches,
            List<DynList<face> >& facesForChain
        );
        void createFacesForChainWithoutCorners
        (
            const labelList& bChain,
            const DynList<label>& patches,
            List<DynList<face> >& facesForChain
        );
        void createBoundaryFacesForCell(const label);

        void createBoundaryFaces();
        
        //- store boundary faces into the mesh
        void storeBoundaryFaces();
        
        //- find surface corners
        void createSurfaceCornersAndPatches();

        //- disallows bitwise construct
        void operator=(const boundaryFacesGenerator&);

        //- copy constructor
        boundaryFacesGenerator(const boundaryFacesGenerator&);

public:

    // Constructors

        //- construct components
        boundaryFacesGenerator
        (
            polyMeshGen& mesh,
            List<direction>& nFacesInCell,
            label& nPoints,
            const label nIntFaces,
            const boolList& boundaryCell,
            VRWGraph& pointPatches,
            const triSurface& surface
        );

        //- Destructor
        ~boundaryFacesGenerator();

    // Member functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
