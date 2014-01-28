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
#include "helperFunctions.H"

//#define DEBUGDual
//#define DEBUGClosedness

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void dualMeshExtractor::decomposeCreatedPoly::decomposeFaces()
{
    //- check if face needs to be decomposed
    forAll(cellFaces_, cfI)
    {
        const face& f = cellFaces_[cfI];

        if( f.size() == 3 )
        {
            decomposedFaces_.append(f);
            splitFace_.append(false);
        }
        else
        {
            List<direction> nodeLevel(4);
            forAll(f, pI)
                nodeLevel[pI] = nodeLevel_.find(f[pI])();

            # ifdef DEBUGDual
            forAll(nodeLevel, nI)
                Info << "level for node " << f[nI] << " is "
                    << label(nodeLevel[nI]) << endl;
            # endif

            direction minLevel(nodeLevel[0]);
            forAll(nodeLevel, nI)
                if( minLevel > nodeLevel[nI] )
                    minLevel = nodeLevel[nI];

            //- find positions of min level nodes
            DynList<label> posMin;
            forAll(nodeLevel, nI)
                if( nodeLevel[nI] == minLevel )
                    posMin.append(nI);

            switch( posMin.size() )
            {
                case 4:
                {
                    //- store face
                    decomposedFaces_.append(f);
                    splitFace_.append(false);
                } break;
                case 3:
                {
                    # ifdef DEBUGDual
                    Info << "2. Decomposing " << f << endl;
                    # endif
                    //- decompose face into 2 triangles
                    direction posHigh(0);
                    forAll(nodeLevel, nI)
                        if( nodeLevel[nI] != minLevel )
                        {
                            posHigh = nI;
                            break;
                        }

                    splitEdge_.append
                    (
                        edge(f[posHigh], f[(posHigh+2)%4])
                    );
                    splitFaceNode_[cfI][posHigh] = true;
                    splitFaceNode_[cfI][(posHigh+2)%4] = true;

                    face tri(3);
                    tri[0] = f[posHigh];
                    tri[1] = f[f.fcIndex(posHigh)];
                    tri[2] = f[(posHigh+2)%4];
                    # ifdef DEBUGDual
                    Info << "2. tri " << tri << endl;
                    # endif
                    decomposedFaces_.append(tri);
                    splitFace_.append(true);

                    tri[0] = f[posHigh];
                    tri[1] = f[(posHigh+2)%4];
                    tri[2] = f[f.rcIndex(posHigh)];
                    # ifdef DEBUGDual
                    Info << "2. tri " << tri << endl;
                    # endif
                    decomposedFaces_.append(tri);
                    splitFace_.append(true);
                } break;
                case 2:
                {
                    if(
                        (f.fcIndex(posMin[0]) == posMin[1])
                    || (f.rcIndex(posMin[0]) == posMin[1])
                    )
                    {
                        decomposedFaces_.append(f);
                        splitFace_.append(false);
                    }
                    else
                    {
                        # ifdef DEBUGDual
                        Info << "3. Decomposing " << f << endl;
                        # endif
                        face tri(3);
                        direction pm = posMin[0];

                        splitEdge_.append
                        (
                            edge(f[f.fcIndex(pm)], f[f.rcIndex(pm)])
                        );
                        splitFaceNode_[cfI][f.fcIndex(pm)] = true;
                        splitFaceNode_[cfI][f.rcIndex(pm)] = true;

                        tri[0] = f[pm];
                        tri[1] = f[f.fcIndex(pm)];
                        tri[2] = f[f.rcIndex(pm)];
                        # ifdef DEBUGDual
                        Info << "3. tri " << tri << endl;
                        # endif
                        decomposedFaces_.append(tri);
                        splitFace_.append(true);

                        pm = posMin[1];
                        tri[0] = f[pm];
                        tri[1] = f[f.fcIndex(pm)];
                        tri[2] = f[f.rcIndex(pm)];
                        # ifdef DEBUGDual
                        Info << "3. tri " << tri << endl;
                        # endif
                        decomposedFaces_.append(tri);
                        splitFace_.append(true);
                    }
                } break;
                case 1:
                {
                    # ifdef DEBUGDual
                    Info << "4. Decomposing " << f << endl;
                    # endif
                    //- decompose face into 2 triangles
                    const direction pm = posMin[0];

                    splitEdge_.append
                    (
                        edge(f[f.fcIndex(pm)], f[f.rcIndex(pm)])
                    );
                    splitFaceNode_[cfI][f.fcIndex(pm)] = true;
                    splitFaceNode_[cfI][f.rcIndex(pm)] = true;

                    face tri(3);
                    tri[0] = f[pm];
                    tri[1] = f[f.fcIndex(pm)];
                    tri[2] = f[f.rcIndex(pm)];
                    # ifdef DEBUGDual
                    Info << "4. tri " << tri << endl;
                    # endif
                    decomposedFaces_.append(tri);
                    splitFace_.append(true);

                    tri[0] = f[f.fcIndex(pm)];
                    tri[1] = f[(pm+2)%4];
                    tri[2] = f[f.rcIndex(pm)];
                    # ifdef DEBUGDual
                    Info << "4. tri " << tri << endl;
                    # endif
                    decomposedFaces_.append(tri);
                    splitFace_.append(true);
                } break;
            };
        }
    }
}

void dualMeshExtractor::decomposeCreatedPoly::calculateHelperAddressing()
{
    dfe_.setSize(decomposedFaces_.size());
    forAll(decomposedFaces_, faceI)
    {
        dfe_[faceI].setSize(decomposedFaces_[faceI].size());
        dfe_[faceI] = -1;
    }

    direction edgeI(0);
    forAll(decomposedFaces_, faceI)
    {
        const edgeList edg = decomposedFaces_[faceI].edges();
        forAll(edg, eI)
        {
            const label pos = de_.containsAtPosition(edg[eI]);

            if( pos == -1 )
            {
                de_.append(edg[eI]);
                dfe_[faceI][eI] = edgeI;
                def_[edgeI][0] = faceI;
                ++edgeI;
            }
            else
            {
                dfe_[faceI][eI] = pos;
                def_[pos][1] = faceI;
            }
        }
    }

    se_.setSize(edgeI);
    se_ = false;

    forAll(de_, eI)
        if( splitEdge_.contains(de_[eI]) )
            se_[eI] = true;

    # ifdef DEBUGDual
    Info << "de_ " << de_ << endl;
    Info << "def_ " << def_ << endl;
    Info << "dfe_ " << dfe_ << endl;
    Info << "se " << se_ << endl;
    # endif
}

void dualMeshExtractor::decomposeCreatedPoly::selectAdditionalSplitEdges()
{
    DynList<label> splitNode(8);
    forAll(splitEdge_, eI)
    {
        splitNode.appendIfNotIn(splitEdge_[eI].start());
        splitNode.appendIfNotIn(splitEdge_[eI].end());
    }

    forAll(se_, eI)
    {
        if( se_[eI] )
            continue;

        const edgeList edges0 = decomposedFaces_[def_[eI][0]].edges();
        const edgeList edges1 = decomposedFaces_[def_[eI][1]].edges();

        const edge& e = de_[eI];

        # ifdef DEBUGDual
        Info << "Checking edge " << eI << " with nodes " << e << endl;
        # endif

        forAll(e, vI)
        {
            if( !splitNode.contains(e[vI]) )
                continue;

            DynList<edge> connectedEdges(2);
            forAll(edges0, eJ)
                if( (edges0[eJ] != e) && (edges0[eJ].otherVertex(e[vI]) != -1) )
                {
                    connectedEdges.append(edges0[eJ]);
                    break;
                }

            forAll(edges1, eJ)
                if( (edges1[eJ] != e) && (edges1[eJ].otherVertex(e[vI]) != -1) )
                {
                    connectedEdges.append(edges1[eJ]);
                    break;
                }

            # ifdef DEBUGDual
            Info << "Connected edges " << connectedEdges << endl;
            # endif

            forAll(cellFaces_, cfI)
            {
                const face& f = cellFaces_[cfI];
                const edgeList fEdges = f.edges();

                label e0(-1), e1(-1), cv(-1);
                forAll(fEdges, feI)
                {
                    if( fEdges[feI] == connectedEdges[0] )
                    {
                        e0 = feI;
                    }
                    else if( fEdges[feI] == connectedEdges[1] )
                    {
                        e1 = feI;
                    }

                    if( f[feI] == e[vI] )
                        cv = feI;
                }

                # ifdef DEBUGDual
                Info << "e0 " << e0 << endl;
                Info << "e1 " << e1 << endl;
                Info << "cv " << cv << endl;
                # endif

                if( (e0 != -1) && (e1 != -1) && splitFaceNode_[cfI][cv] )
                {
                    se_[eI] = true;
                    break;
                }
            }
        }
    }
}

void dualMeshExtractor::decomposeCreatedPoly::selectFacesForCell
(
    List<faceList>& cFaces
)
{
    if( splitEdge_.size() == 0 )
    {
        //- cell has not been split
        cFaces.setSize(1);
        cFaces[0] = cellFaces_;
    }
    else
    {
        direction cellI(0);

        cFaces.setSize(5);
        boolList storedFace(decomposedFaces_.size(), false);

        forAll(storedFace, faceI)
        {
            if( storedFace[faceI] )
                continue;

            DynList<label> front;
            front.append(faceI);

            # ifdef DEBUGDual
            Info << "Starting front with face " << faceI << endl;
            # endif

            cFaces.newElmt(cellI).setSize(5);
            direction cfI(0);

            do
            {
                DynList<label> newFront;

                forAll(front, fI)
                {
                    const label fLabel = front[fI];
                    if( storedFace[fLabel] )
                        continue;
                    cFaces[cellI].newElmt(cfI++) = decomposedFaces_[fLabel];
                    storedFace[fLabel] = true;

                    # ifdef DEBUGDual
                    Info << "Storing face " << decomposedFaces_[fLabel]
                        << " into cell " << label(cellI) << endl;
                    # endif

                    forAll(dfe_[fLabel], eI)
                    {
                        const label eLabel = dfe_[fLabel][eI];
                        if( !se_[eLabel] )
                        {
                            if( !storedFace[def_[eLabel][1]] )
                            {
                                newFront.append(def_[eLabel][1]);
                                # ifdef DEBUGDual
                                Info << "Adding face " << def_[eLabel][1]
                                    << " into the front" << endl;
                                # endif
                            }
                            else if( !storedFace[def_[eLabel][0]] )
                            {
                                newFront.append(def_[eLabel][0]);
                                # ifdef DEBUGDual
                                Info << "Adding face " << def_[eLabel][0]
                                    << " into the front" << endl;
                                # endif
                            }
                        }
                    }
                }

                front = newFront;
            } while( front.size() );

            cFaces[cellI].setSize(cfI);
            ++cellI;
        }

        cFaces.setSize(cellI);
        # ifdef DEBUGDual
        Info << "Decomposed faces " << decomposedFaces_ << endl;
        Info << "Cell faces " << cFaces << endl;
        forAll(storedFace, fI)
            if( !storedFace[fI] )
                FatalErrorIn
                (
                    "void dualMeshExtractor::decomposeCreatedPoly::"
                    "selectFacesForCell(List<faceList>&)"
                ) << "Face " << fI << " is not stored!" << abort(FatalError);
        # endif

        if( cellI < 2 )
        {
            //selectFacesForCellComplex(cFaces);

            if( cFaces.size() < 2 )
                FatalErrorIn
                (
                    "void dualMeshExtractor::decomposeCreatedPoly::"
                    "selectFacesForCell(List<faceList>&)"
                ) << "Cell is not decomposed, but should be!"
                    << abort(FatalError);
        }
    }
}

void dualMeshExtractor::decomposeCreatedPoly::createMissingFaces
(
    List<faceList>& cFaces
)
{
    if( splitEdge_.size() == 0 )
        return;

    //- missing faces will be created from open edges shared by two
    //- adjacent cells
    # ifdef DEBUGDual
    Info << "Creating missing faces for cells " << cFaces << endl;
    # endif
    List<direction> nFacesInCell(cFaces.size());
    forAll(cFaces, cI)
        nFacesInCell[cI] = cFaces[cI].size();
    List<DynList<edge> > openEdges(cFaces.size());

    forAll(cFaces, cI)
    {
        help::findOpenEdges(cFaces[cI], openEdges[cI]);
    }

    # ifdef DEBUGDual
    Info << "Open edges are " << openEdges << endl;
    # endif

    forAll(openEdges, cellI)
    {
        const DynList<edge>& oe = openEdges[cellI];
        for(label cellJ=(cellI+1);cellJ<openEdges.size();++cellJ)
        {
            const DynList<edge>& oes = openEdges[cellJ];

            DynList<edge> commonEdges(4);
            forAll(oes, eI)
                if( oe.contains(oes[eI]) )
                    commonEdges.append(oes[eI]);

            # ifdef DEBUGDual
            Info << "Common edges are " << commonEdges << endl;
            # endif

            if( commonEdges.size() > 1 )
            {
                //- zip open chains
                help::zipOpenChain(commonEdges);

                //- create a face
                face newBf(help::sortEdgeChain(commonEdges));
                # ifdef DEBUGDual
                Info << "Created face for cells " << label(cellI) << " and "
                    << label(cellJ) << " is " << newBf << endl;
                # endif

                //- store face into new cells
                cFaces[cellI].newElmt(nFacesInCell[cellI]++) = newBf;
                cFaces[cellJ].newElmt(nFacesInCell[cellJ]++) =
                    newBf.reverseFace();
            }
        }
    }

    forAll(cFaces, cI)
    {
        cFaces[cI].setSize(nFacesInCell[cI]);
        help::findOpenEdges(cFaces[cI], openEdges[cI]);
    }

    # ifdef DEBUGDual
    Info << "Open edges making internal cell are " << openEdges << endl;
    # endif
    faceList midCell(4);
    direction fI(0);

    forAll(openEdges, cI)
        if( openEdges[cI].size() > 1 )
        {
            //- zip open chain
            help::zipOpenChain(openEdges[cI]);

            //- create new face
            face f(help::sortEdgeChain(openEdges[cI]));
            midCell[fI++] = f;
            cFaces[cI].newElmt(nFacesInCell[cI]++) = f.reverseFace();
            # ifdef DEBUGDual
            Info << "1. Adding face of internal cell " << f << endl;
            # endif
        }

    if( fI > 3 )
    {
        midCell.setSize(fI);

        cFaces.setSize(nFacesInCell.size()+1);
        cFaces[nFacesInCell.size()] = midCell;
    }

    forAll(nFacesInCell, cI)
        cFaces[cI].setSize(nFacesInCell[cI]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dualMeshExtractor::decomposeCreatedPoly::decomposeCreatedPoly
(
    const faceList& cf,
    const Map<direction>& nl
)
:
    cellFaces_(cf),
    nodeLevel_(nl),
    splitFaceNode_(cf.size()),
    splitEdge_(6),
    splitFace_(2*cf.size()),
    decomposedFaces_(2*cf.size()),
    dfe_(12),
    def_(24, labelList(2)),
    de_(24),
    se_(24)
{
    forAll(cf, fI)
    {
        splitFaceNode_[fI].setSize(cf[fI].size());
        splitFaceNode_[fI] = true;
    }
}

dualMeshExtractor::decomposeCreatedPoly::~decomposeCreatedPoly()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualMeshExtractor::decomposeCreatedPoly::decomposeCell
(
    List<faceList>& decCells
)
{
    //decCells.setSize(1);
    //decCells[0] = cellFaces_;

    # ifdef DEBUGDual
    Info << nl << "Starting decomposing cell " << cellFaces_ << endl;
    # endif

    decomposeFaces();

    if( splitEdge_.size() == 0 )
    {
        decCells.setSize(1);
        decCells[0] = cellFaces_;
        return;
    }

    calculateHelperAddressing();

    selectAdditionalSplitEdges();

    selectFacesForCell(decCells);

    createMissingFaces(decCells);

    # ifdef DEBUGClosedness
    //- check for topological closedness
    forAll(decCells, cI)
    {
        const faceList& cf = decCells[cI];

        DynList<edge> cEdges(12);
        DynList<direction> nAppearances(12);

        forAll(cf, fI)
        {
            const edgeList edges = cf[fI].edges();

            forAll(edges, eI)
            {
                const label pos = cEdges.containsAtPosition(edges[eI]);

                if( pos == -1 )
                {
                    cEdges.append(edges[eI]);
                    nAppearances.append(1);
                }
                else
                {
                    nAppearances[pos]++;
                }
            }
        }

        forAll(nAppearances, eI)
            if( nAppearances[eI] != 2 )
            {
                Info << "Original cell is " << cellFaces_ << endl;
                Info << "cell faces are " << decCells << endl;
                Info << "Edges for cell " << cI << " are "
                    << cEdges << endl;
                FatalErrorIn
                (
                    "void dualMeshExtractor::decomposeCreatedPoly::"
                    "createMissingFaces(List<faceList>&)"
                ) << "Edge " << eI << " is open!!" << abort(FatalError);
            }
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
