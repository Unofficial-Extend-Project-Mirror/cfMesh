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

Namespace
    Foam::edgeExtraction

Description
    Helper functions used for capturing feature edges
    and transferring of surface patches on the volume mesh

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef edgeExtraction_H
#define edgeExtraction_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



namespace Foam
{

/*---------------------------------------------------------------------------*\
                Namespace EdgeExtraction functions Declaration
\*---------------------------------------------------------------------------*/

namespace EdgeExtraction
{
    //- project face centres on the nearest location at the surface mesh
    //- and assign the patch to the patch of the surface element
    void distributeBoundaryFaces
    (
        const meshSurfaceEngine&,
        const meshOctree&,
        labelListPMG&
    );

    //- find the nearest points on the surface of the volume mesh
    //- to the corners on the triangulated surface
    void findCornerCandidates
    (
        const meshSurfaceEngine&,
        const meshOctree&
        LongList<labelPair>&
    );



} // End namespace EdgeExtraction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
