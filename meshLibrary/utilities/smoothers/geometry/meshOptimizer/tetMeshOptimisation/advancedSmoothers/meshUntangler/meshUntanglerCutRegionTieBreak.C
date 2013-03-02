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

#include "demandDrivenData.H"
#include "meshUntangler.H"
#include "sortEdgesIntoChains.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::cutRegion::tieBreak(const DynList<label, 8>& f)
{
    // There must be only one, singly connected region which is cut off
    // from the feasible region. This may not be the case if one operates in the
    // range of computer tolerances. In order to resolve the tie it is necessary
    // to find a single region which will be cut off from the feasible region.
    // This will be done by finding the node with the most negative distance,
    // and then start marking vertices connected to that vertex via an edge.
    
    # ifdef DEBUGSmooth
    Info << "Starting tie break" << endl;
    # endif
    
    //- delete pointer data
    deleteDemandDrivenData(cPtsPtr_);
    deleteDemandDrivenData(cEdgesPtr_);
    deleteDemandDrivenData(cFacesPtr_);
    
    //- remove coincident vertices
    //removeCoincidentVertices();
    
    const DynList<edge, 128>& edges = *edgesPtr_;

    DynList<edge> faceEdges(f.size());
    forAll(f, eI)
        faceEdges.append(edges[f[eI]]);
        
    labelListList fvertices = sortEdgesIntoChains(faceEdges).sortedChains();
    if( fvertices.size() != 1 )
    {
        valid_ = false;
        return;
        
        Info << "Face vertices " << fvertices << endl;
        FatalErrorIn
        (
            "void meshUntangler::cutRegion::tieBreak(const face& f)"
        ) << "Number of created faces is not 1 but "
            << fvertices.size() << abort(FatalError);
    }
    
    const labelList& fv = fvertices[0];
    
    DynList<label, 64> vertexRegion;
    vertexRegion.setSize(fv.size());
    vertexRegion = 0;
    
    label region(1);
    forAll(fv, vI)
        if( !vertexTypes_[fv[vI]] && !vertexRegion[vI] )
        {
            vertexRegion[vI] = region;
            
            label fcI = fv.fcIndex(vI);
            label rcI = fv.rcIndex(vI);
            bool found;
            do
            {
                found = false;
                if( !vertexTypes_[fv[fcI]] )
                {
                    vertexRegion[fcI] = region;
                    fcI = fv.fcIndex(fcI);
                    found = true;
                }
                
                if( !vertexTypes_[fv[rcI]] )
                {
                    vertexRegion[rcI] = region;
                    rcI = fv.rcIndex(rcI);
                    found = true;
                }
            } while( found );
            
            ++region;
        }
        
    # ifdef DEBUGSmooth
    Info << "Tolerance " << tol_ << endl;
    Info << "Number of regions " << region-1 << endl;
    Info << "Vertex regions " << vertexRegion << endl;
    # endif
    
    if( region > 2 )
    {
        //- there are more than two regions which need to be cut off
        # ifdef DEBUGSmooth
        forAll(fv, vI)
            Info << "Distance for vertex " << fv[vI] << " is "
                << vertexDistance_[fv[vI]] << endl;
        # endif
        
        //- there should be only one cut-off region
        //- there this region will be determined by the most negative
        //- distance from plane
        scalar minDist(VGREAT);
        label minRegion(-1);
        forAll(fv, vI)
            if( vertexRegion[vI] && (vertexDistance_[fv[vI]] < minDist) )
            {
                minDist = vertexDistance_[fv[vI]];
                minRegion = vertexRegion[vI];
            }
            
        forAll(vertexRegion, vI)
            if( vertexRegion[vI] && (vertexRegion[vI] != minRegion) )
            {
                vertexTypes_[fv[vI]] |= INPLANE;
            }
    }
    else
    {
        forAll(fv, vI)
            if(
                (vertexTypes_[fv[vI]] & INPLANE) &&
                !vertexRegion[fv.rcIndex(vI)] &&
                !vertexRegion[fv.fcIndex(vI)]
            )
            {
                vertexTypes_[fv[vI]] ^= INPLANE;
                vertexTypes_[fv[vI]] |= KEEP;
                
                # ifdef DEBUGSmooth
                Info << "Node " << vI << " was INPLANE" << endl;
                Info << "New type " << label(vertexTypes_[fv[vI]]) << endl;
                # endif
            }
    }
    
    # ifdef DEBUGSmooth
    forAll(fv, vI)
        Info << "Vertex type for vertex " << fv[vI] << " is "
            << label(vertexTypes_[fv[vI]]) << endl;
    # endif
    
    //- create new points
    const DynList<point, 64>& points = *pointsPtr_;
    cPtsPtr_ = new DynList<point, 64>();
    newVertexLabel_ = -1;
    origNumVertices_ = 0;
    forAll(points, pI)
        if( vertexTypes_[pI] )
        {
            cPtsPtr_->append(points[pI]);
            newVertexLabel_[pI] = origNumVertices_++;
        }
        
    //- find new edges and continue creating faces
    findNewEdges();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
