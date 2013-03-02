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

#include "voronoiMeshExtractor.H"
#include "polyMeshGenModifierAddCellByCell.H"
#include "tessellationElement.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

//#define DEBUGVoronoi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createPoints()
{
    const LongList<point>& tetPoints = tetCreator_.tetPoints();
    const LongList<partTet>& tets = tetCreator_.tets();
    
    pointFieldPMG& points = mesh_.points();
    points.setSize(tets.size());
    
    # ifdef DEBUGVoronoi
    Info << "Number of tets " << tets.size() << endl;
    # endif
    
    forAll(tets, tetI)
    {
        points[tetI] = tets[tetI].centroid(tetPoints);
        
        # ifdef DEBUGVoronoi
        Info << "Point " << tetI << " has coordinates "
            << points[tetI] << endl;
        Info << "Tet of origin " << tetI << " has nodes "
            << tets[tetI] << endl;
        # endif
    }
}

void voronoiMeshExtractor::createPolyMesh()
{
    const VRWGraph& pointEdges = this->pointEdges();
    const VRWGraph& edgeTets = this->edgeTets();
    const boolList& boundaryEdge = this->boundaryEdge();
    const LongList<edge>& edges = this->edges();
    const LongList<partTet>& tets = tetCreator_.tets();
    
    polyMeshGenModifierAddCellByCell meshModifier(mesh_);
    
    forAll(pointEdges, pointI)
    {
        bool create(true);
        forAllRow(pointEdges, pointI, eI)
            if( boundaryEdge[pointEdges(pointI, eI)] )
            {
                create = false;
                break;
            }
        
        if( !create || (pointEdges.sizeOfRow(pointI) == 0) )
            continue;
        
        faceList cellFaces(pointEdges.sizeOfRow(pointI));
        
        forAllRow(pointEdges, pointI, faceI)
        {
            const label edgeI = pointEdges(pointI, faceI);
            
            //- check if the face orientation needs to be changed
            bool flip(false);
            const partTet& tet = tets[edgeTets(edgeI, 0)];
            const partTet& tet1 = tets[edgeTets(edgeI, 1)];
            tessellationElement tEl(tet[0], tet[1], tet[2], tet[3]);
            for(label i=0;i<4;++i)
            {
                const triFace tf = tEl.face(i);
                label nShared(0);
                for(label j=0;j<3;++j)
                    if( tet1.whichPosition(tf[j]) != -1 )
                        ++nShared;
                
                if( nShared == 3 )
                {
                    const edge& e = edges[edgeI];
                    
                    const edgeList tEdges = tf.edges();
                    forAll(tEdges, teI)
                        if( tEdges[teI] == e )
                        {
                            if( pointI == tEdges[teI].start() )
                                flip = true;
                            
                            break;
                        }
                    
                    break;
                }
            }
            
            face& f = cellFaces[faceI];
            f.setSize(edgeTets.sizeOfRow(edgeI));
            
            //- fill the faces with the node labels
            forAll(f, pI)
                f[pI] = edgeTets(edgeI, pI);
            
            if( flip )
                f = f.reverseFace();
        }
        
        meshModifier.addCell(cellFaces);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
