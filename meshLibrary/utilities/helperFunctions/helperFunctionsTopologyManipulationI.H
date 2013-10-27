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

#include "error.H"
#include "helperFunctionsTopologyManipulation.H"
#include "edgeList.H"
#include "pointField.H"
#include "boolList.H"
#include "cellList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

template<class labelListType1, class labelListType2>
bool areFacesEqual(const labelListType1& f1, const labelListType2& f2)
{
    //- their size mut be equal
    if( f1.size() != f2.size() )
        return false;

    //- find the starting point the the second face
    label start(-1);
    const label s = f2.size();
    bool equalOrientation(false);

    forAll(f2, pI)
    {
        if( f1[0] == f2[pI] )
        {
            start = pI;

            if( f1[1] == f2[(pI+1)%s] )
            {
                //- faces have the same orientation
                equalOrientation = true;
            }
            else if( f1[1] != f2[(s-1+pI)%s] )
            {
                //- faces shall have opposite orientation, but do not match
                return false;
            }
        }
    }

    if( start < 0 )
        return false;

    if( equalOrientation )
    {
        //- same orientation, check if all points match
        for(label pI=1;pI<s;++pI)
        {
            if( f1[pI] != f2[(pI+start)%s] )
                return false;
        }
    }
    else
    {
        //- faces with opposite orientation, check if all points match
        for(label pI=1;pI<s;++pI)
        {
            if( f1[pI] != f2[(start+s-pI)%s] )
                return false;
        }
    }

    //- faces are equal
    return true;
}

template<class labelListType>
labelListType reverseFace(const labelListType& f)
{
    labelListType rf;
    rf.setSize(f.size());

    rf[0] = f[0];

    const label size = f.size();
    for(label i=1;i<size;++i)
        rf[f.size()-i] = f[i];

    return rf;
}

inline face mergeTwoFaces(const face& f1, const face& f2)
{
    const edgeList e1 = f1.edges();
    const edgeList e2 = f2.edges();

    boolList ce1(e1.size(), false), ce2(e2.size(), false);

    forAll(e1, eI)
        forAll(e2, eJ)
        if( e1[eI] == e2[eJ] )
        {
            ce1[eI] = true;
            ce2[eJ] = true;
            break;
        }

    DynList<edge> fEdges;
    forAll(ce1, eI)
        if( !ce1[eI] ) fEdges.append(e1[eI]);
    forAll(ce2, eI)
        if( !ce2[eI] ) fEdges.append(e2[eI]);

    return face(sortEdgeChain(fEdges));
}

inline edgeList modifyFacesToShareOneEdge(face& f1, face& f2)
{
    const edgeList edg1 = f1.edges();
    const edgeList edg2 = f2.edges();

    boolList she1(f1.size(), false);
    boolList she2(f2.size(), false);

    edgeList sharedEdges(edg1);
    label nSharedEdges(0);

    forAll(edg1, eI)
        forAll(edg2, eJ)
            if( edg1[eI] == edg2[eJ] )
            {
                sharedEdges.newElmt(nSharedEdges++) = edg1[eI];

                she1[eI] = true;
                she2[eJ] = true;
                break;
            }

    face newF(f1);
    label i(0);
    forAll(f1, pI)
        if( !(she1[pI] && she1[(pI-1+f1.size())%f1.size()]) )
            newF[i++] = f1[pI];

    newF.setSize(i);
    f1 = newF;

    newF = f2;
    i = 0;
    forAll(f2, pI)
        if( !(she2[pI] && she2[(pI-1+f2.size())%f2.size()]) )
            newF[i++] = f2[pI];

    newF.setSize(i);
    f2 = newF;

    sharedEdges.setSize(nSharedEdges);
    return sharedEdges;
}

inline face createFaceFromRemovedPart(const face& fOrig, const face& fCut)
{
    if( fCut.size() == 0 )
        return fOrig;

    const edgeList eOrig = fOrig.edges();
    const edgeList eCut = fCut.edges();

    boolList usedEdge(eOrig.size(), false);

    forAll(eOrig, eI)
        forAll(eCut, eJ)
            if( eOrig[eI] == eCut[eJ] )
            {
                usedEdge[eI] = true;
                break;
            }

    face f(fOrig);
    direction i(0);

    forAll(fOrig, pI)
        if( !(usedEdge[pI] && usedEdge[(pI-1+fOrig.size())%fOrig.size()]) )
        {
            f[i++] = fOrig[pI];
        }

    f.setSize(i);

    return f;
}

inline face removeEdgesFromFace
(
    const face& fOrig,
    const DynList<edge>& removeEdges
)
{
    boolList foundEdge(fOrig.size(), false);

    forAll(removeEdges, reI)
        forAll(fOrig, eI)
            if( removeEdges[reI] == fOrig.faceEdge(eI) )
            {
                foundEdge[eI] = true;
                break;
            }

    face newF(fOrig.size());
    label i(0);

    forAll(fOrig, pI)
        if( !(foundEdge[pI] && foundEdge[fOrig.rcIndex(pI)]) )
            newF[i++] = fOrig[pI];

    newF.setSize(i);

    return newF;
}

inline void findOpenEdges(const faceList& cellFaces, DynList<edge>& openEdges)
{
    DynList<edge> cellEdges(12);
    DynList<label> nAppearances(12);

    forAll(cellFaces, fI)
    {
        const edgeList edges = cellFaces[fI].edges();

        forAll(edges, eI)
        {
            const label pos = cellEdges.containsAtPosition(edges[eI]);

            if( pos == -1 )
            {
                cellEdges.append(edges[eI]);
                nAppearances.append(1);
            }
            else
            {
                nAppearances[pos]++;
            }
        }
    }

    openEdges.setSize(12);
    openEdges.clear();

    forAll(nAppearances, eI)
        if( nAppearances[eI] == 1 )
        {
            openEdges.append(cellEdges[eI]);
        }
        else if( nAppearances[eI] > 2 )
        {
            FatalErrorIn
            (
                "void findOpenEdges(const faceList& cellFaces,"
                "DynList<edge>& openEdges)"
            ) << "More than two faces in " << cellFaces
                << " share edge " << cellEdges[eI] << abort(FatalError);
        }
}

inline bool shareAnEdge(const face& f1, const face& f2)
{
    const edgeList edges1 = f1.edges();
    const edgeList edges2 = f2.edges();

    forAll(edges1, eI)
        forAll(edges2, eJ)
        if( edges1[eI] == edges2[eJ] )
            return true;

    return false;
}

inline edge sharedEdge(const face& f1, const face& f2)
{
    const edgeList edges1 = f1.edges();
    const edgeList edges2 = f2.edges();

    forAll(edges1, eI)
        forAll(edges2, eJ)
        if( edges1[eI] == edges2[eJ] )
            return edges1[eI];

    return edge(-1, -1);
}

inline label positionOfEdgeInFace(const edge& e, const face& f)
{
    const edgeList edges = f.edges();

    forAll(edges, eI)
        if( edges[eI] == e )
            return eI;

    return -1;
}

inline label sharedVertex(const face& f1, const face& f2)
{
    forAll(f1, pI)
        forAll(f2, pJ)
            if( f1[pI] == f2[pJ] )
                return f1[pI];

    return -1;
}

inline bool shareAVertex(const face& f1, const face& f2)
{
    forAll(f1, pI)
        forAll(f2, pJ)
        if( f1[pI] == f2[pJ] )
            return true;

    return false;
}

inline label sharedVertex(const faceList& fcs)
{
    const face& f = fcs[0];

    forAll(f, pI)
    {
        bool allFound(true);

        for(label i=1;i<fcs.size();i++)
        {
            bool found(false);
            forAll(fcs[i], pJ)
                if( f[pI] == fcs[i][pJ] )
                    found = true;

            if( !found )
            {
                allFound = false;
                break;
            }
        }

        if( allFound )
            return f[pI];
    }

    return -1;
}

inline bool areElementsInChain(const boolList& sel)
{
    boolList selInChain(sel.size(), false);

    forAll(sel, eI)
    {
        if( sel[eI] )
        {
            selInChain[eI] = true;
            bool found;
            do
            {
                found = false;
                forAll(selInChain, eJ)
                    if(
                        !selInChain[eJ] && sel[eJ] &&
                        (
                            selInChain[sel.fcIndex(eJ)] ||
                            selInChain[sel.rcIndex(eJ)]
                        )
                    )
                    {
                        found = true;
                        selInChain[eJ] = true;
                    }
            } while( found );

            break;
        }
    }

    forAll(sel, eI)
    {
        if( sel[eI] && !selInChain[eI] )
            return false;
    }

    return true;
}

inline void zipOpenChain(DynList<edge>& bEdges)
{
    //- close the chain if open
    DynList<label> chainVertices(5);
    List<label> nAppearances(5);
    forAll(bEdges, eI)
    {
        const edge& e = bEdges[eI];
        forAll(e, pI)
        {
            const label pos = chainVertices.containsAtPosition(e[pI]);

            if( pos == -1 )
            {
                nAppearances[chainVertices.size()] = 1;
                chainVertices.append(e[pI]);
            }
            else
            {
                ++nAppearances[pos];
            }
        }
    }

    bool closed(true);
    DynList<label> openVertices(2);
    forAll(chainVertices, pI)
        if( nAppearances[pI] == 1 )
        {
            closed = false;
            openVertices.append(chainVertices[pI]);
        }

    if( !closed && (openVertices.size() == 2) )
    {
        forAll(bEdges, eI)
            if( bEdges[eI].end() == openVertices[0] )
            {
                bEdges.append(edge(openVertices[0], openVertices[1]));
                break;
            }
            else if( bEdges[eI].end() == openVertices[1] )
            {
                bEdges.append(edge(openVertices[1], openVertices[0]));
                break;
            }
    }
    else if( !closed )
    {
        FatalErrorIn
        (
            "void dualMeshExtractor::decomposeCreatedPoly::"
            "createMissingFaces(List<faceList>& cFaces)"
        ) << "Chain has " << openVertices << " open vertices"
            << abort(FatalError);
    }
}

inline labelList sortEdgeChain(const DynList<edge>& bEdges)
{
    boolList sorted(bEdges.size(), false);

    DynList<edge> sortedEdges;
    sortedEdges.append(bEdges[0]);
    sorted[0] = true;
    direction i(0);

    bool finished;
    do
    {
        finished = true;

        forAll(bEdges, eI)
            if( !sorted[eI] )
            {
                if( sortedEdges[i].end() == bEdges[eI].start() )
                {
                    sorted[eI] = true;
                    finished = false;
                    sortedEdges.append(bEdges[eI]);
                    ++i;
                }
                else if( sortedEdges[i].end() == bEdges[eI].end() )
                {
                    FatalErrorIn
                    (
                        "labelList sortEdgeChain("
                        "const DynList<edge>& bEdges)"
                    ) << "Chain is not oriented correctly!"
                        << abort(FatalError);
                }
            }
    } while( !finished );

    labelList sortPoints(bEdges.size());
    forAll(sortedEdges, eI)
        sortPoints[eI] = sortedEdges[eI].start();

    return sortPoints;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
