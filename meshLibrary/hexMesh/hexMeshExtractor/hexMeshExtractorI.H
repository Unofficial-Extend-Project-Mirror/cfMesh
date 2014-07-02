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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void hexMeshExtractor::createHexCells::createFaceCentreHexes
(
    const label pointI
)
{
    //- find the maximum level of octree cubes at this point
    bool create(true);
    direction maxLevel(0);
    label pos(-1), startNode(-1);
    forAllRow(nodeLeaves_, pointI, nlI)
    {
        const label leafI = nodeLeaves_(pointI, nlI);
        
        if( (leafI < 0) || (centreNodeLabel_[leafI] < 0) )
        {
            create = false;
            break;
        }
        
        if( maxLevel < octree_.returnLeaf(leafI).level() )
        {
            maxLevel = octree_.returnLeaf(leafI).level();
            startNode = 7 - nlI;
            pos = nlI;
        }
    }
    
    if( !create )
        return;
    
    //- create cells from vertices
    FixedList<label, 8> cellPoints;
    const meshOctreeCubeCoordinates cc =
        octree_.returnLeaf(nodeLeaves_(pointI, pos));
    
    forAllRow(nodeLeaves_, pointI, nlI)
    {
        const label leafI = nodeLeaves_(pointI, nlI);
        
        if( octree_.returnLeaf(leafI).level() == maxLevel )
        {
            //- store the centre node
            cellPoints[nlI] = centreNodeLabel_[leafI];
        }
        else
        {
            const meshOctreeCubeCoordinates& oc =
                octree_.returnLeaf(leafI).coordinates();
            
            //- store the correct sub vertex
            const FixedList<meshOctreeCubeCoordinates, 8>& vlPos =
                octree_.positionsOfLeavesAtNodes()[startNode];
            const meshOctreeCubeCoordinates sc = cc + vlPos[nlI];
            
            label lpos(-1);
                        
            for(label j=0;j<8;j++)
            {
                if( sc == oc.refineForPosition(j) )
                {
                    lpos = j;
                    break;
                }
            }
            
            cellPoints[nlI] = subVertices_(leafI, lpos);
            
            if( cellPoints[nlI] < 0 )
                cellPoints[nlI] = centreNodeLabel_[leafI];
        }
    }
    
    # ifdef DEBUGHex
    Info << "Points for cell are " << cellPoints << endl;
    # endif
    
    //- create a cell and add it to the mesh
    addHexCell(cellPoints);
}

inline void hexMeshExtractor::createHexCells::createHexesAtOctreePoints
(
    const label pointI
)
{
    if( octreeVertexType_[pointI] != FACECENTRE )
        return;
    
    //- find the maximum refinement level of boxes at this point
    direction maxLevel(0);
    forAllRow(nodeLeaves_, pointI, nlI)
    {
        const label leafI = nodeLeaves_(pointI, nlI);
        
        if( leafI < 0 )
            continue;
        
        maxLevel = Foam::max(maxLevel, octree_.returnLeaf(leafI).level());
    }
    
    //- find the direction with all four boxes at low refinement level
    for(label fI=0;fI<6;++fI)
    {
        const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];
        
        //- check if all 4 leaves are at low ref level
        bool create(true);
        for(label i=0;i<4;++i)
        {
            const label leafI = nodeLeaves_(pointI, fNodes[i]);
            
            if( (leafI < 0) || (centreNodeLabel_[leafI] < 0) )
            {
                create = false;
                break;
            }
            else if( octree_.returnLeaf(leafI).level() == maxLevel )
            {
                create = false;
                break;
            }
        }
        
        if( create )
        {
            # ifdef DEBUGHex
            Info << "Creating FACECENTRE hexes at level "
                << label(maxLevel)-1 << endl;
            Info << "fI " << fI << endl;
            # endif
            const label ofI = meshOctreeCubeCoordinates::oppositeFace_[fI];
            const label* ofNodes =
                meshOctreeCubeCoordinates::faceNodes_[ofI];
            
            # ifdef DEBUGHex
            Info << "ofI " << ofI << endl;
            # endif
            
            FixedList<label, 4> centres, prev, curr, next;
            for(label i=0;i<4;++i)
            {
                const label leafI = nodeLeaves_(pointI, fNodes[i]);
                
                centres[i] = centreNodeLabel_[leafI];
                
                label pos(-1);
                for(label j=0;j<4;++j)
                {
                    if( nodeLabels_(leafI, ofNodes[j]) == pointI )
                    {
                        pos = j;
                        break;
                    }
                }
                
                prev[i] = subVertices_(leafI, ofNodes[(pos+3)%4]);
                curr[i] = subVertices_(leafI, ofNodes[pos]);
                next[i] = subVertices_(leafI, ofNodes[(pos+1)%4]);
            }
            
            # ifdef DEBUGHex
            Info << "Centre points " << centres << endl;
            Info << "Curr subvertices " << curr << endl;
            Info << "Prev subvertices " << prev << endl;
            Info << "Next subvertices " << next << endl;
            # endif
                
            //- create four hex cells, one from subVertices, and the second
            //- one from subvertices and centre nodes
            FixedList<label, 8> cellPoints;
            
            //- first cell
            cellPoints[0] = curr[0];
            cellPoints[1] = curr[1];
            cellPoints[2] = curr[3];
            cellPoints[3] = curr[2];
            cellPoints[4] = prev[0];
            cellPoints[5] = next[1];
            cellPoints[6] = next[3];
            cellPoints[7] = prev[2];
            
            addHexCell(cellPoints);
            
            //- second cell
            cellPoints[0] = next[0];
            cellPoints[1] = prev[1];
            cellPoints[2] = curr[0];
            cellPoints[3] = curr[1];
            cellPoints[4] = centres[0];
            cellPoints[5] = centres[1];
            cellPoints[6] = prev[0];
            cellPoints[7] = next[1];
            
            addHexCell(cellPoints);
            
            //- third cell
            cellPoints[0] = curr[3];
            cellPoints[1] = curr[2];
            cellPoints[2] = prev[3];
            cellPoints[3] = next[2];
            cellPoints[4] = next[3];
            cellPoints[5] = prev[2];
            cellPoints[6] = centres[3];
            cellPoints[7] = centres[2];
            
            addHexCell(cellPoints);
            
            //- fourth cell
            cellPoints[0] = prev[0];
            cellPoints[1] = next[1];
            cellPoints[2] = next[3];
            cellPoints[3] = prev[2];
            cellPoints[4] = centres[0];
            cellPoints[5] = centres[1];
            cellPoints[6] = centres[3];
            cellPoints[7] = centres[2];
            
            addHexCell(cellPoints);
            
            break;
        }
    }
}

inline void hexMeshExtractor::createHexCells::createEdgeHexes
(
    const label pointI
)
{
    if( octreeVertexType_[pointI] != MIXED )
        return;
    
    //- find the maximum refinement level of boxes at this point
    direction maxLevel(0);
    forAllRow(nodeLeaves_, pointI, nlI)
    {
        const label leafI = nodeLeaves_(pointI, nlI);
        
        if( leafI < 0 )
            continue;
        
        maxLevel = Foam::max(maxLevel, octree_.returnLeaf(leafI).level());
    }
    
    //- find the direction with all four boxes at low refinement level
    for(label fI=0;fI<6;++fI)
    {
        const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];
        
        //- check if all 4 leaves are at low ref level
        bool create(true);
        for(label i=0;i<4;++i)
        {
            const label leafI = nodeLeaves_(pointI, fNodes[i]);
            
            if( (leafI < 0) || (centreNodeLabel_[leafI] < 0) )
            {
                create = false;
                break;
            }
            else if( octree_.returnLeaf(leafI).level() == maxLevel )
            {
                create = false;
                break;
            }
        }
        
        if( create )
        {
            //- check if there exist any duplicate cubes
            if( nodeLeaves_(pointI, fNodes[0]) ==
                nodeLeaves_(pointI, fNodes[1])
            )
                continue;
            if( nodeLeaves_(pointI, fNodes[0]) ==
                nodeLeaves_(pointI, fNodes[3])
            )
                continue;
            
            FixedList<label, 8> cellPoints;
            //- create a cell from subVertices an centre nodes
            for(label i=0;i<4;++i)
            {
                const label leafI = nodeLeaves_(pointI, fNodes[mapping_[i]]);
                
                cellPoints[i] = subVertices_(leafI, 7-fNodes[mapping_[i]]);
                if( cellPoints[i] < 0 )
                {
                    create = false;
                    break;
                }
                cellPoints[i+4] = centreNodeLabel_[leafI];
            }
            
            if( !create )
                continue;
            
            # ifdef DEBUGHex
            Info << "Edge hex cell has points " << cellPoints << endl;
            # endif
            
            addHexCell(cellPoints);
        }
    }
}

inline void hexMeshExtractor::createHexCells::addHexCell
(
    const FixedList<label, 8>& cellPoints
)
{
    faceList hexCell(6, face(4));
    forAll(hexCell, fI)
    {
        face& f = hexCell[fI];
        const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];

        forAll(f, pI)
            f[pI] = cellPoints[fNodes[pI]];
    }
    
    # ifdef DEBUGHex
    Info << "Adding cell " << hexCell << endl;
    # endif
    
    meshModifier_.addCell(hexCell);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
