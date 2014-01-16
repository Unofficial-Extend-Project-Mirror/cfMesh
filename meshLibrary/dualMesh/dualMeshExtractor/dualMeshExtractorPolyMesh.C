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

#include "dualMeshExtractor.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "polyMeshGenModifierAddCellByCell.H"
#include "labelLongList.H"

//#define DEBUGDual

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualMeshExtractor::createPoints()
{
    clearOut();
    
    const List<direction>& boxType = octreeAddressing_.boxType();
    const meshOctree& octree = octreeAddressing_.octree();
    centreNodeLabelPtr_ = new labelLongList(boxType.size(), -1);
    labelLongList& centreNode = *centreNodeLabelPtr_;
    
    const boundBox& rootBox = octree.rootBox();
    label nPoints(0);
    
    forAll(centreNode, boxI)
        if( boxType[boxI] & meshOctreeAddressing::MESHCELL )
            centreNode[boxI] = nPoints++;
        
    pointFieldPMG& points = mesh_.points();
    points.setSize(nPoints);
        
    forAll(centreNode, boxI)
        if( centreNode[boxI] != -1 )
        {
            points[centreNode[boxI]] =
                octree.returnLeaf(boxI).centre(rootBox);
        }
        
}

void dualMeshExtractor::createPolyMesh()
{
    Info << "Creating polyMesh from octree" << endl;
    
    const meshOctree& octree = octreeAddressing_.octree();
    const labelLongList& centreNode = *centreNodeLabelPtr_;
    const FRWGraph<label, 8>& nodeLeaves = octreeAddressing_.nodeLeaves();
    
    polyMeshGenModifierAddCellByCell meshModifier(mesh_);
    
    forAll(nodeLeaves, nodeI)
    {
        bool create(true);
        
        forAllRow(nodeLeaves, nodeI, nlI)
        {
            const label leafI = nodeLeaves(nodeI, nlI);
            
            if( (leafI == -1) || (centreNode[leafI] == -1) )
            {
                create = false;
                break;
            }
        }
        
        if( !create )
            continue;
        
        Map<direction> nodeLevel;
        forAllRow(nodeLeaves, nodeI, nlI)
            nodeLevel.insert
            (
                centreNode[nodeLeaves(nodeI, nlI)],
                octree.returnLeaf(nodeLeaves(nodeI, nlI)).level()
            );
        
        faceList cFaces(12);
        direction fI(0);
        
        for(label i=0;i<6;++i)
        {
            DynList<label> f(4);
            
            for(label j=0;j<4;++j)
            {
                f.appendIfNotIn
                (
                    centreNode[nodeLeaves(nodeI, faceFlip_[i][j])]
                );
            }
            
            if( f.size() == 4 )
            {
                cFaces.newElmt(fI++) = face(f);
            }
            else if( f.size() == 3 )
            {
                f.shrink();
                cFaces.newElmt(fI++) = face(f);
            }
        }
        
        cFaces.setSize(fI);
        
        # ifdef DEBUGDual
        Info << "pLeaves " << pLeaves << endl;
        forAll(pLeaves, plI)
            Info << "Centre " << plI << " is "
                << centreNode[pLeaves[plI]] << endl;
        Info << "Cell faces are " << cFaces << nl << endl;
        forAll(cFaces, cfI)
            Info << "Face " << cfI << " normal is "
                << cFaces[cfI].normal(mesh_.points()) << endl;
        # endif
            
        List<faceList> tmpF;
        decomposeCreatedPoly dcp(cFaces, nodeLevel);
        dcp.decomposeCell(tmpF);
        forAll(tmpF, cI)
        {
            # ifdef DEBUGDual
            Info << cI << ". Cell " << cellI << " has faces "
                << tmpF[cI] << endl;
            # endif
            
            meshModifier.addCell(tmpF[cI]);
        }
    }

    Info << "Finished creating polyMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
