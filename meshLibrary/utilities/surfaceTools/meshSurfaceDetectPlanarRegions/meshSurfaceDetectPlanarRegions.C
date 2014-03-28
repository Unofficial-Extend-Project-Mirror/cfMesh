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

#include "meshSurfaceDetectPlanarRegions.H"
#include "meshSurfaceEngine.H"
#include "labelledPoint.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace meshSurfaceDetectPlanarRegionsHelper
{

class meshSurfaceDetectPlanarRegionsNeiOp
{
    // Private data
        //- const reference to meshSurfaceEngine
        const meshSurfaceEngine& mse_;

        //- angle tolerance
        const scalar cosTol_;

        //- face normals
        vectorField faceNormals_;

public:

    // Constructors
        //- Construct from meshSurfaceEngine
        meshSurfaceDetectPlanarRegionsNeiOp
        (
            const meshSurfaceEngine& mse,
            const scalar cosTol
        )
        :
            mse_(mse),
            cosTol_(cosTol),
            faceNormals_()
        {
            //- create face-faces addressing outside of a parallel region
            mse_.faceFaces();

            //- calculate face normals
            const pointFieldPMG& pts = mse_.points();
            const faceList::subList& bFaces = mse_.boundaryFaces();
            faceNormals_.setSize(bFaces.size());

            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 50)
            # endif
            forAll(bFaces, bfI)
            {
                const face& bf = bFaces[bfI];

                vector& fn = faceNormals_[bfI];
                fn = bf.normal(pts);
                fn /= mag(fn);
            }
        }

    // Public member functions
        //- return the size (number of boundary faces)
        inline label size() const
        {
            return mse_.boundaryFaces().size();
        }

        //- find neighbours of a boundary face
        void operator()(const label bfI, DynList<label>& neiFaces) const
        {
            neiFaces.clear();
            const VRWGraph& faceFaces = mse_.faceFaces();

            const vector& fn = faceNormals_[bfI];
            forAllRow(faceFaces, bfI, ffI)
            {
                const vector& on = faceNormals_[faceFaces(bfI, ffI)];

                if( (on & fn) > cosTol_ )
                    neiFaces.append(faceFaces(bfI, ffI));
            }
        }

        //- unify region in case of an MPI parallel run
        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynList<label> >& neiGroups,
            const labelListType& elementInGroup,
            const DynList<label>& localGroupLabel
        ) const
        {
        /*
            const Map<label>& globalToLocal =
                mse_.globalToLocalBndEdgeAddressing();
            const Map<label>& otherFaceAtProc = mse_.otherEdgeFaceAtProc();

            const DynList<label>& neiProcs = mse_.beNeiProcs();

            std::map<label, DynList<labelledPoint> > exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label beI = it();

                forAllRow(bpEdges, bpI, i)
                {
                    const label beI = bpEdges(bpI, i);

                    if( !isFeatureEdge_[beI] )
                        continue;

                    const label groupI = elementInGroup[beI];

                    forAllRow(bpAtProcs, bpI, ppI)
                    {
                        const label neiProc = bpAtProcs(bpI, ppI);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append
                        (
                            labelPair(it.key(), localGroupLabel[groupI])
                        );
                    }
                }
            }

            LongList<labelPair> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelPair& lp = receivedData[i];
                const label groupI = elementInGroup[globalToLocal[lp.first()]];

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(lp.second());
            }
            */
        }
};

class meshSurfaceDetectPlanarRegionsSelOp
{
public:

    // Constructors
        //- Null construct
        meshSurfaceDetectPlanarRegionsSelOp()
        {}

    // Public member functions
        //- a dummy operator
        bool operator()(const label) const
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshSurfaceDetectPlanarRegionsHelper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceDetectPlanarRegions::findPlanarRegions()
{
    const scalar cosTol = Foam::cos(angleTol_);

    meshSurfaceDetectPlanarRegionsHelper::meshSurfaceDetectPlanarRegionsNeiOp nop
    (
        meshSurface_,
        cosTol
    );

    meshSurfaceDetectPlanarRegionsHelper::meshSurfaceDetectPlanarRegionsSelOp sop;

    nPlanarRegions_ = help::groupMarking(faceInPlanarRegion_, nop, sop);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceDetectPlanarRegions::meshSurfaceDetectPlanarRegions
(
    const meshSurfaceEngine& meshSurface,
    const scalar angleTol
)
:
    meshSurface_(meshSurface),
    angleTol_(angleTol),
    faceInPlanarRegion_(),
    nPlanarRegions_(0)
{
    findPlanarRegions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceDetectPlanarRegions::~meshSurfaceDetectPlanarRegions()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
