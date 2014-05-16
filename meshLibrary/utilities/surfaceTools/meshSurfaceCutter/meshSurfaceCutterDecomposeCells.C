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

#include "meshSurfaceCutter.H"
#include "helperFunctions.H"

// #define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void meshSurfaceCutter::checkCellTopology()
{
    //- check if there exist a cluster of cells which are
    //- stored as a single cell
    //- this may happen when there exist two or more internal faces originating
    //- from a single poly face
    bool decompose(false);

    # ifdef DEBUGCutter
    Info << "Problematic topology " << problematicTopology_ << endl;
    # endif

    forAll(problematicTopology_, cI)
        if( problematicTopology_[cI] )
        {
            decompose = true;
            break;
        }

    if( decompose )
    {
        SLList<cell> newCells;
        forAll(problematicTopology_, cI)
            if( problematicTopology_[cI] )
            {
                const cell& c = polyCells_[cI];

                # ifdef DEBUGCutter
                Info << "Cell " << cI << " may consist of two or more cells."
                    << nl << "cell faces are " << c << endl;
                # endif

                //- create cell addressing
                const edgeList edges = c.edges(polyFaces_);
                List< DynList<label> > edgeFaces(edges.size());

                forAll(c, fI)
                {
                    const edgeList fe = polyFaces_[c[fI]].edges();

                    forAll(edges, eI)
                        forAll(fe, feI)
                        if( edges[eI] == fe[feI] )
                        {
                            edgeFaces[eI].append(fI);
                            break;
                        }
                }

                forAll(edgeFaces, eI)
                    if( edgeFaces[eI].size() != 2 )
                    {
                        Info << "edgeFaces " << edgeFaces << endl;

                        FatalErrorIn
                        (
                            "void meshSurfaceCutter::checkCellTopology()"
                        ) << "Cell " << cI << " is not closed"
                            << abort(FatalError);
                    }

                labelListList faceFaces(c.nFaces());
                forAll(c, fI)
                    faceFaces[fI].setSize(polyFaces_[c[fI]].size());
                labelList count(c.size(), 0);

                forAll(edgeFaces, eI)
                {
                    const DynList<label>& ef = edgeFaces[eI];
                    faceFaces[ef[0]][count[ef[0]]++] = ef[1];
                    faceFaces[ef[1]][count[ef[1]]++] = ef[0];
                }
                
                //- decide which faces belong to which cell
                labelList cellFace(c.nFaces(), -1);

                bool finished;

                short i(0);
                cellFace[0] = i;

                do
                {
                    do
                    {
                        finished = true;
                        forAll(c, fI)
                            if( cellFace[fI] != -1 )
                                forAll(faceFaces[fI], fJ)
                                {
                                    if( cellFace[faceFaces[fI][fJ]] == -1 )
                                    {
                                        finished = false;
                                        cellFace[faceFaces[fI][fJ]] =
                                            cellFace[fI];
                                    }
                                    else if
                                    (
                                        cellFace[fI] !=
                                        cellFace[faceFaces[fI][fJ]]
                                    )
                                    {
                                        FatalErrorIn
                                        (
                                            "void meshSurfaceCutter::checkCellTopology()"
                                        ) << "Cannot handle this!!"
                                            << abort(FatalError);
                                    }
                                }
                    } while( !finished );

                    forAll(cellFace, cfI)
                        if( cellFace[cfI] == -1 )
                        {
                            i++;
                            cellFace[cfI] = i;
                            finished = false;
                            break;
                        }

                } while( !finished );

                cellList decCells(i+1, c);
                labelList nF(decCells.size(), 0);

                //- store faces into correct cells
                forAll(cellFace, cfI)
                    decCells[cellFace[cfI]].newElmt
                    (nF[cellFace[cfI]]++) = c[cfI];

                //- store cells into the list
                forAll(nF, cI)
                {
                    decCells[cI].setSize(nF[cI]);
                    # ifdef DEBUGCutter
                    Info << "Cell part " << cI << " is "
                        << decCells[cI] << endl;
                    # endif
                    newCells.append(decCells[cI]);
                }
            }
            else
            {
                newCells.append(polyCells_[cI]);
            }
        
        //- store cells into the list
        label i(0);
        polyCells_.setSize(newCells.size());
        for(SLList<cell>::const_iterator cIter = newCells.begin();
            cIter != newCells.end();
            ++cIter
        )
        polyCells_[i++] = cIter();
    }

    problematicTopology_.clear();
    /*
    Info << "Checking for faces which share more than two vertices" << endl;
    bool finished;
    do
    {
        finished = true;
        //- check if there exist faces which share more than two vertices
        List< DynList<label> > pFaces(nPoints_);
        
        forAll(polyFaces_, fI)
        {
            const face& f = polyFaces_[fI];
            
            forAll(f, pI)
                pFaces[f[pI]].append(fI);
        }
        
        facesFromFace_.clear();
        facesFromFace_.setSize(polyFaces_.size());
        
        for(label fI=nIntFaces_;fI<polyFaces_.size();fI++)
        {
            const edgeList edges = polyFaces_[fI].edges();
            
            const face& f = polyFaces_[fI];
            
            forAll(f, pI)
            {
                const DynList<label>& pf = pFaces[f[pI]];
                
                forAll(pf, fJ)
                    if( (pf[fJ] < fI) && facesFromFace_[pf[fJ]].size() < 2 )
                    {
                        const edgeList edg = polyFaces_[pf[fJ]].edges();
                        
                        DynList<edge> ce;
                        
                        forAll(edges, eI)
                            forAll(edg, eJ)
                            if( edges[eI] == edg[eJ] )
                            {
                                ce.append(edges[eI]);
                                break;
                            }
                        
                        if( ce.size() > 1 )
                        {
                            # ifdef DEBUGCutter
                            Info << "Faces " << fI << " and "
                                << pf[fJ] << " share "
                                << ce << " edges!" << endl;
                            # endif

                            forAll(ce, ceI)
                                for(label ceJ=ceI+1;ceJ<ce.size();ceJ++)
                                {
                                    const label cv = ce[ceI].commonVertex(ce[ceJ]);
                                    if( cv != -1 )
                                    {
                                        if( isFaceConvex(fI) )
                                        {
                                            decomposeFaceIntoTriangles(fI, cv);
                                        }
                                        else if( concaveVertex(fI) == cv )
                                        {
                                            decomposeFaceIntoTriangles(fI);
                                        }
                                        else
                                        {
                                            # ifdef DEBUGCutter
                                            Info << "Decomposing face " << fI
                                                << " manually" << endl;
                                            # endif

                                            decomposeFace(fI,cv);

                                            # ifdef DEBUGCutter
                                            Info << "Face " << fI
                                                << " is decomposed into "
                                                << facesFromFace_[fI] << endl;
                                            # endif
                                        }
                                    }
                                }
                        }
                    }
                
                if( facesFromFace_[fI].size() > 1 )
                {
                    finished = false;
                    break;
                }
            }
        }
        
        faceList newFaces(2*polyFaces_.size());
        label counter(0);
        labelListList childs(polyFaces_.size());
        
        forAll(polyFaces_, faceI)
        {
            if(facesFromFace_[faceI].size() > 1 )
            {
                childs[faceI].setSize(facesFromFace_[faceI].size());
                forAll(facesFromFace_[faceI], fJ)
                {
                    childs[faceI][fJ] = counter;
                    newFaces.newElmt(counter++) = facesFromFace_[faceI][fJ];
                }
            }
            else
            {
                childs[faceI].setSize(1);
                childs[faceI][0] = counter;
                newFaces.newElmt(counter++) = polyFaces_[faceI];
            }
        }
        
        facesFromFace_.clear();
        
        newFaces.setSize(counter);
        
        polyFaces_ = newFaces;
        newFaces.clear();
        
        //- change cells
        forAll(polyCells_, cellI)
        {
            cell& c = polyCells_[cellI];
            
            cell copy(c);
            short i(0);
            
            forAll(c, fI)
                forAll(childs[c[fI]], chI)
                copy.newElmt(i++) = childs[c[fI]][chI];
            
            copy.setSize(i);
            
            c = copy;
        }
        
        //- change boundary information
        for(short patchI=1;patchI<nFacesInPatch_.size();patchI++)
        {
            const labelList& chE = childs[patchStart_[patchI]];
            nFacesInPatch_[patchI-1] = chE[chE.size()-1] - patchStart_[patchI-1];
            patchStart_[patchI] = patchStart_[patchI-1] + nFacesInPatch_[patchI-1];
        }
        
        nFacesInPatch_[nFacesInPatch_.size()-1] =
            polyFaces_.size() - patchStart_[patchStart_.size()-1];

    } while( !finished );
        */
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
