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

#include "hexMeshExtractor.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "polyMeshGenModifierAddCellByCell.H"
#include "labelListPMG.H"

//#define DEBUGHex

# ifdef DEBUGHex
#include "writeMeshFPMA.H"
#include "writeMeshFLMA.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void hexMeshExtractor::classifyOctreePoints()
{
    const VRWGraph& nodeLabels = octreeAddressing_.nodeLabels();
    const FRWGraph<label, 8>& nl = octreeAddressing_.nodeLeaves();
    const meshOctree& octree = octreeAddressing_.octree();
    
    octreeVertexType_.setSize(nl.size());
    octreeVertexType_ = NONE;
    
    forAll(nl, pointI)
    {
        if( octreeVertexType_[pointI] )
            continue;
        
        direction maxLevel(0), minLevel(255);
        
        bool boundaryNode(false);
        forAllRow(nl, pointI, nlI)
        {
            const label leafI = nl(pointI, nlI);
            
            if( leafI < 0 )
            {
                boundaryNode = true;
                continue;
            }
            
            maxLevel = Foam::max(maxLevel, octree.returnLeaf(leafI).level());
            minLevel = Foam::min(minLevel, octree.returnLeaf(leafI).level());
        }
        
        if( boundaryNode )
        {
            octreeVertexType_[pointI] = BOUNDARY;
            continue;
        }
        
        if( maxLevel == minLevel )
        {
            //- all octree leaves are at the same level. Point is a corner
            octreeVertexType_[pointI] = CORNER;
        }
        else
        {
            # ifdef DEBUGHex
            Info << "\nClassifying point " << pointI << endl;
            forAllRow(nl, pointI, nlI)
            {
                const label leafI = nl(pointI, nlI);
                
                Info << "Leaf at position " << nlI
                    << " is " << nl(pointI, nlI) << endl;
                
                if( leafI >= 0 )
                    Info << "Leaf coordinates "
                        << octree.returnLeaf(leafI) << endl;
            }
            #endif
            
            //- check if the octree point is at a face centre of a father
            //- cube of the boxes found in the directions of a face
            //- there must exist 4 different boxes at minLevel
            //- which have the same father
            for(label fI=0;fI<6;++fI)
            {
                const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];
                
                //- check if the vertex is a centre over a single face or edge
                if( 
                    (nl(pointI, fNodes[0]) == nl(pointI, fNodes[1])) ||
                    (nl(pointI, fNodes[1]) == nl(pointI, fNodes[2])) ||
                    (nl(pointI, fNodes[2]) == nl(pointI, fNodes[0])) ||
                    (nl(pointI, fNodes[0]) == nl(pointI, fNodes[3]))
                )
                    continue;
                
                //- check if all cubes are at minLevel and the nodes
                //- of the opposite face shall be at maxLevel
                bool check(false);
                for(label i=0;i<4;++i)
                {
                    const meshOctreeCubeBasic& oca =
                        octree.returnLeaf(nl(pointI, fNodes[i]));
                    
                    if( oca.level() != minLevel )
                    {
                        check = true;
                        break;
                    }
                    
                    const meshOctreeCubeBasic& ocb =
                        octree.returnLeaf(nl(pointI, 7-fNodes[i]));
                    
                    if( ocb.level() != maxLevel )
                    {
                        check = true;
                        break;
                    }
                }
                
                if( check )
                    continue;
                
                //- check if the cubes have the same father box
                check = true;
                const meshOctreeCubeBasic& oc =
                        octree.returnLeaf(nl(pointI, fNodes[0]));
                const meshOctreeCubeCoordinates father = oc.reduceLevelBy(1);
                for(label i=1;i<4;++i)
                {
                    const meshOctreeCubeBasic& oca =
                        octree.returnLeaf(nl(pointI, fNodes[i]));
                    
                    if( oca.reduceLevelBy(1) != father )
                    {
                        check = false;
                        break;
                    }
                }
                
                if( check )
                {
                    # ifdef DEBUGHex
                    Info << "Point " << pointI << " is a face centre" << endl;
                    # endif
                    
                    octreeVertexType_[pointI] = FACECENTRE;
                    
                    //- mark opposite vertices at all 4 faces as NOSUBVERTICES
                    const label ofI =
                        meshOctreeCubeCoordinates::oppositeFace_[fI];
                    const label* ofNodes =
                        meshOctreeCubeCoordinates::faceNodes_[ofI];
                    for(label i=0;i<4;++i)
                    {
                        const label leafJ = nl(pointI, fNodes[i]);
                        
                        label pos(-1);
                        for(label j=0;j<4;++j)
                            if( nodeLabels(leafJ, ofNodes[j]) == pointI )
                            {
                                pos = j;
                                break;
                            }
                            
                        if( pos == -1 )
                            FatalErrorIn
                            (
                                "void hexMeshExtractor::classifyOctreePoints()"
                            ) << "Cannot find vertex" << abort(FatalError);

                        
                        const label pJ = nodeLabels(leafJ, ofNodes[(pos+2)%4]);
                        octreeVertexType_[pJ] = NOSUBVERTICES;
                    }
                }
            }
                
            if( !octreeVertexType_[pointI] )
            {
                # ifdef DEBUGHex
                Info << "Point " << pointI << " is of type MIXED" << endl;
                # endif
                
                octreeVertexType_[pointI] = MIXED;
            }
        }
    }
    
    # ifdef DEBUGHex
    Info << "Writting debug mesh" << endl;
    Info << "Number of octree leaves is " << octree.numberOfLeaves() << endl;
    polyMeshGen pmg(mesh_.returnTime());
    
    const List<direction>& boxType = octreeAddressing_.boxType();
    const pointField& octreePoints = octreeAddressing_.octreePoints();
    const VRWGraph& octreeFaces = octreeAddressing_.octreeFaces();
    const VRWGraph& boxFaces = octreeAddressing_.leafFaces();
    const labelListPMG& neighbour = octreeAddressing_.octreeFaceNeighbour();
    
    pointFieldPMG& points = polyMeshGenModifier(pmg).pointsAccess();
    points.setSize(octreePoints.size());
    forAll(octreePoints, pI)
        points[pI] = octreePoints[pI];
    
    polyMeshGenModifierAddCellByCell* bModPtr =
        new polyMeshGenModifierAddCellByCell(pmg);
    forAll(boxType, leafI)
    {
        if( !(boxType[leafI] & meshOctreeAddressing::MESHCELL) )
            continue;
        
        faceList c(boxFaces.sizeOfRow(leafI));
        
        forAllRow(boxFaces, leafI, fI)
        {
            const label faceI = boxFaces(leafI, fI);
            
            face f(octreeFaces.sizeOfRow(faceI));
            forAll(f, pI)
                f[pI] = octreeFaces(faceI, pI);
            
            if( neighbour[faceI] == leafI )
                f = f.reverseFace();
            
            c[fI].transfer(f);
        }
        
        bModPtr->addCell(c);
    }
    
    deleteDemandDrivenData(bModPtr);
    
    Info << "Mesh has " << pmg.faces().size() << " faces" << endl;
    Info << "Mesh has " << pmg.cells().size() << " cells" << endl;
    
    polyMeshGenModifier(pmg).reorderBoundaryFaces();
    
    const label cornerID = pmg.addPointSubset("CORNER");
    const label faceID = pmg.addPointSubset("FACECENTRE");
    const label edgeID = pmg.addPointSubset("NOSUBVERTICES");
    const label bndVrt = pmg.addPointSubset("BOUNDARY");
    
    forAll(octreeVertexType_, pI)
    {
        switch( octreeVertexType_[pI] )
        {
            case CORNER:
            {
                pmg.addPointToSubset(cornerID, pI);
            } break;
            case FACECENTRE:
            {
                pmg.addPointToSubset(faceID, pI);
            } break;
            case NOSUBVERTICES:
            {
                pmg.addPointToSubset(edgeID, pI);
            } break;
            case BOUNDARY:
            {
                pmg.addPointToSubset(bndVrt, pI);
            } break;
        };
    }
    
    writeMeshFPMA(pmg, "markedPoints");
    # endif
}
    
void hexMeshExtractor::createPoints()
{
    clearOut();

    const List<direction>& boxType = octreeAddressing_.boxType();
    const meshOctree& octree = octreeAddressing_.octree();
    const VRWGraph& nodeLabels = octreeAddressing_.nodeLabels();
    const FRWGraph<label, 8>& nodeLeaves = octreeAddressing_.nodeLeaves();
    centreNodeLabelPtr_ = new labelListPMG(boxType.size(), -1);
    labelListPMG& centreNode = *centreNodeLabelPtr_;
    subVerticesPtr_ = new VRWGraph(boxType.size());
    VRWGraph& subVertices = *subVerticesPtr_;

    const boundBox& rootBox = octree.rootBox();
    label nPoints(0);

    //- create vertices and subvertices
    forAll(centreNode, leafI)
    {
        if( boxType[leafI] & meshOctreeAddressing::MESHCELL )
        {
            //- create a vertex in the centre of the octree cube
            centreNode[leafI] = nPoints++;

            //- create subvertices near required corners
            forAllRow(nodeLabels, leafI, nI)
            {
                const label nodeI = nodeLabels(leafI, nI);
                
                //- find the maximum level of octree cubes at this node
                direction maxLevel(0);
                
                forAllRow(nodeLeaves, nodeI, i)
                {
                    const label leafLabel = nodeLeaves(nodeI, i);
                    
                    if( leafLabel < 0 )
                        continue;
                    
                    const meshOctreeCubeBasic& oc =
                        octree.returnLeaf(leafLabel);
                    
                    maxLevel = Foam::max(maxLevel, oc.level());
                }
                    
                if( maxLevel != octree.returnLeaf(leafI).level() )
                {
                    if( subVertices.sizeOfRow(leafI) == 0 )
                    {
                        subVertices.setRowSize(leafI, 8);
                        for(label i=0;i<8;++i)
                            subVertices(leafI, i) = -1;
                    }
                    
                    if( octreeVertexType_[nodeI] & NOSUBVERTICES )
                        continue;
                    
                    subVertices(leafI, nI) = nPoints++;
                }
            }
        }
    }
    
    //- set the number of points
    pointFieldPMG& points = mesh_.points();
    points.setSize(nPoints);
    
    //- create points
    forAll(centreNode, leafI)
    {
        if( centreNode[leafI] != -1 )
        {
            //- create centre node
            points[centreNode[leafI]] =
                octree.returnLeaf(leafI).centre(rootBox);
            
            forAllRow(subVertices, leafI, i)
            {
                const label svI = subVertices(leafI, i);
                
                if( svI < 0 )
                    continue;
                
                //- create the subvertex
                const meshOctreeCubeCoordinates cc =
                    octree.returnLeaf(leafI).refineForPosition(i);
                points[svI] = cc.centre(rootBox);
                
                if( octreeVertexType_[nodeLabels(leafI, i)] & FACECENTRE )
                {
                    const label* pFaces =
                        meshOctreeCubeCoordinates::nodeFaces_[i];
                    for(label j=0;j<3;++j)
                    {
                        DynList<label> neighs;
                        octree.findNeighboursInDirection
                        (
                            leafI,
                            pFaces[j],
                            neighs
                        );
                        
                        if( neighs.size() == 4 )
                        {
                            const scalar disp = 0.5 * cc.size(rootBox);
                            
                            switch( pFaces[j] )
                            {
                                case 0:
                                {
                                    points[svI].x() -= disp;
                                } break;
                                case 1:
                                {
                                    points[svI].x() += disp;
                                } break;
                                case 2:
                                {
                                    points[svI].y() -= disp;
                                } break;
                                case 3:
                                {
                                    points[svI].y() += disp;
                                } break;
                                case 4:
                                {
                                    points[svI].z() -= disp;
                                } break;
                                case 5:
                                {
                                    points[svI].z() += disp;
                                } break;
                            };
                        }
                    }
                }
            }
        }
    }
    
    # ifdef DEBUGHex
    polyMeshGen pmg(mesh_.returnTime());
    
    const FixedList<Vector<label>, 8>& octantVectors = octree.octantVectors();
    polyMeshGenModifier meshModifier(pmg);
    const label centID = pmg.addCellSubset("centrePoints");
    label setID = pmg.addCellSubset("subVertices");
    forAll(subVertices, leafI)
    {
        if( centreNode[leafI] < 0 )
            continue;
        
        const scalar s = 0.1 * octree.returnLeaf(leafI).size(rootBox);
        const point cent = octree.returnLeaf(leafI).centre(rootBox);
        label nPoints = pmg.points().size();
        label nFaces = pmg.faces().size();
        label nCells = pmg.cells().size();
        
        forAll(octree.octantVectors(), ovI)
        {
            point p;
            for(direction i=0;i<Vector<label>::nComponents;++i)
                p[i] = cent[i] + s * octantVectors[ovI][i];
            
            meshModifier.pointsAccess().append(p);
        }
        
        cell c(6);
        
        forAll(c, fI)
        {
            face f(4);
            const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];
            
            forAll(f, pI)
                f[pI] = nPoints + fNodes[pI];
            
            meshModifier.facesAccess().append(f);
            c[fI] = nFaces + fI;
        }
        
        meshModifier.cellsAccess().append(c);
        pmg.addCellToSubset(centID, nCells);
        
        forAllRow(subVertices, leafI, svI)
        {
            const label pI = subVertices(leafI, svI);
            
            if( pI < 0 )
                continue;
            
            // add vertices
            nPoints = pmg.points().size();
            nFaces = pmg.faces().size();
            nCells = pmg.cells().size();
            
            const meshOctreeCubeCoordinates cc =
                octree.returnLeaf(leafI).refineForPosition(svI);
            const scalar s = 0.1 * cc.size(rootBox);
            const point cent = cc.centre(rootBox);
            
            forAll(octree.octantVectors(), ovI)
            {
                point p;
                for(direction i=0;i<Vector<label>::nComponents;++i)
                    p[i] = cent[i] + s * octantVectors[ovI][i];
                
                meshModifier.pointsAccess().append(p);
            }
            
            cell c(6);
            
            forAll(c, fI)
            {
                face f(4);
                const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[fI];
                
                forAll(f, pI)
                    f[pI] = nPoints + fNodes[pI];
                
                meshModifier.facesAccess().append(f);
                c[fI] = nFaces + fI;
            }
            
            meshModifier.cellsAccess().append(c);
            pmg.addCellToSubset(setID, nCells);
        }
    }
    meshModifier.reorderBoundaryFaces();
    
    writeMeshFLMA(pmg, "subVerticesAndCentres");
    # endif
}

void hexMeshExtractor::createHexMesh()
{
    createHexCells hexCells
    (
        mesh_,
        octreeAddressing_.octree(),
        *centreNodeLabelPtr_,
        *subVerticesPtr_,
        octreeVertexType_,
        octreeAddressing_.nodeLabels(),
        octreeAddressing_.nodeLeaves()
    );
    
    hexCells.generateCells();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

hexMeshExtractor::createHexCells::createHexCells
(
    polyMeshGen& mesh,
    const meshOctree& octree,
    const labelListPMG& centreNodeLabel,
    const VRWGraph& subVertices,
    const List<direction>& octreeVertexType,
    const VRWGraph& nodeLabels,
    const FRWGraph<label, 8>& nodeLeaves
)
:
    meshModifier_(mesh),
    octree_(octree),
    centreNodeLabel_(centreNodeLabel),
    subVertices_(subVertices),
    octreeVertexType_(octreeVertexType),
    nodeLabels_(nodeLabels),
    nodeLeaves_(nodeLeaves),
    mapping_()
{
    mapping_[0] = 0;
    mapping_[1] = 1;
    mapping_[2] = 3;
    mapping_[3] = 2;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
hexMeshExtractor::createHexCells::~createHexCells()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void hexMeshExtractor::createHexCells::generateCells()
{
    forAll(nodeLeaves_, nodeI)
    {
        createFaceCentreHexes(nodeI);
        
        createHexesAtOctreePoints(nodeI);
        
        createEdgeHexes(nodeI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
