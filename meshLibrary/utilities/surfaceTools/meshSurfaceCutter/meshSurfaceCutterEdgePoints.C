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
#include "trianglePlaneIntersections.H"
#include "helperFunctions.H"

#define DEBUGCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCutter::findBoundaryEdgePoints
(
    const label faceI,
    const face& newF,
    const boolList& outwardPoint,
    const boolList& inwardPoint
)
{
    DynList<point>& newPoints = *newPointsPtr_;

    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& meshPoints = mesh_.points();

    # ifdef DEBUGCutter
    Info << "Creating edge points for face " << faceI << endl;
    Info << "Face consists of vertices " << faces[faceI] << endl;
    const pointField fvrt = faces[faceI].points(meshPoints);
    Info << "Vertices " << fvrt << endl;
    # endif

    label inward(-1);
    short start(0);
    forAll(newF, pI)
        if( inwardPoint[pI] )
        {
            inward = pointTriIndex_[newF[pI]][0];
            start = pI;
            break;
        }

    faceListList& facesFromFace_ = *facesFromFacePtr_;

    if( inward > -1 )
    {
        const face& f = faces[faceI];

        vector normal = f.normal(meshPoints);
        normal /= mag(normal);
        #ifdef DEBUGCutter
        Info << "Face normal " << normal << endl;
        # endif

        const pointField& points = surface_.points();
        const LongList<labelledTri>& faces = surface_.facets();
        const VRWGraph& fe = surface_.facetEdges();
        const VRWGraph& ef = surface_.edgeFacets();

        bool finished;

        short faceJ(0);
        short facePointI(0);

        faceList& newFaces = facesFromFace_[faceI];
        newFaces.setSize(1);
        newFaces[0] = face(15);

        forAll(newF, npJ)
        {
            //- look for points on the boundary edges
            //- if the current point is an outwardPoint and the next point
            //- is an inward point
            const short npI = (npJ+start) % newF.size();
            const short next = (npI+1) % newF.size();

            //- store the current vertex
            newFaces[faceJ].newElmt(facePointI++) = newF[npI];

            if(
                outwardPoint[npI] && inwardPoint[next]
            )
            {
                label lastTri(-1);
                label searchTri = pointTriIndex_[newF[npI]][0];
                const label inwardNext = pointTriIndex_[newF[next]][0];

                if(
                    searchTri == inwardNext
                )
                {
                    continue;
                }

                finished = false;

                labelList intersectedNeighbours(2);
                pointField intersectionPoints(2);
                List< DynList<label> > iPointPatches(2);

                while( !finished )
                {
                    intersectedNeighbours = -1;

                    short counter(0);

                    trianglePlaneIntersections tpi
                    (
                        normal,
                        meshPoints[f[0]],
                        surface_,
                        searchTri
                    );

                    const boolList& intersectedPoints = tpi.intersectedPoints();
                    const boolList& intersectedEdges = tpi.intersectedEdges();
                    const pointField& edgePoints = tpi.edgePoints();

                    # ifdef DEBUGCutter
                    Info << "SEarch tri " << searchTri << "   "
                        << faces[searchTri] << endl;
                    Info << "Vertices in the triangle "
                        << faces[searchTri].points(points) << endl;
                    Info << "intersected edges " << intersectedEdges << endl;
                    Info << "Intersections " << edgePoints << endl;
                    Info << "intersected vertices "
                        << intersectedPoints << endl;
                    # endif

                    forAll(intersectedEdges, eI)
                        if( intersectedPoints[eI] )
                        {
                            # ifdef DEBUGCutter
                            Info << "Plane intersects vertex " << endl;
                            # endif

#                           include "meshSurfaceCutterEdgesVrtImpl.H"
                        }
                        else if( intersectedEdges[eI] )
                        {
                            const label curEdge = fe(searchTri, eI);

                            label neighbourTri = ef(curEdge, 0);
                            if( neighbourTri == searchTri )
                                neighbourTri = ef(curEdge, 1);

                            # ifdef DEBUGCutter
                            Info << "Neighbour tri " << neighbourTri << "  "
                                << faces[neighbourTri] << endl;
                            # endif

                            //- check if the vertex is in two patches
                            DynList<label> patches;
                            patches.append(searchTri);
                            if(
                                faces[searchTri].region() !=
                                faces[neighbourTri].region()
                            )
                                patches.append(neighbourTri);

                            //- store the information
                            iPointPatches[counter] = patches;
                            intersectedNeighbours[counter] = neighbourTri;
                            intersectionPoints[counter++] = edgePoints[eI];
                        }

                    # ifdef DEBUGCutter
                    Info << "Intersected neighbours "
                        << intersectedNeighbours << endl;
                    Info << "intersectionPoints " << intersectionPoints << endl;
                    Info << "iPointPatches " << endl;
                    forAll(iPointPatches, ipI)
                        Info << iPointPatches[ipI] << endl;
                    Info << "lastTri " << lastTri << endl;
                    # endif

                    //- move to the next triangle
                    if( lastTri == -1 )
                    {
                        //- determine rotation
                        lastTri = searchTri;

                        bool in0 =
                            help::pointInsideFace
                            (intersectionPoints[0], f, normal, meshPoints);
                        bool in1 =
                            help::pointInsideFace
                            (intersectionPoints[1], f, normal, meshPoints);

                        if( in0 && !in1 )
                        {
                            searchTri = intersectedNeighbours[0];
                        }
                        else if( !in0 && in1 )
                        {
                            searchTri = intersectedNeighbours[1];
                        }
                        else if( in0 && in1 )
                        {
                            //- both vertices are within the face
                            //- this means that one of them lies on the
                            //- edge of the boundary triangle
                            if(
                                tpi.triInfluence() ==
                                trianglePlaneIntersections::TWO_VERTICES
                            )
                               FatalErrorIn
                                (
                                    "void meshSurfaceCutter::"
                                    "findBoundaryEdgePoints()"
                                ) << " Cannot close face " << faceI
                                    << " because your data is higly degenerate"
                                    << ". Please smooth your triangulated"
                                    << " surface and try again!"
                                    << abort(FatalError);

                            const scalar d0 =
                                mag
                                (
                                    intersectionPoints[0] -
                                    newPoints[newF[npI]]
                                );

                            const scalar d1 =
                                mag
                                (
                                    intersectionPoints[1] -
                                    newPoints[newF[npI]]
                                );

                            if( (d1 - d0) > SMALL )
                            {
                                searchTri = intersectedNeighbours[1];
                            }
                            else if( (d0 - d1) > SMALL )
                            {
                                searchTri = intersectedNeighbours[0];
                            }


                            if( mag(d0 - d1) < SMALL )
                            {
                                Info << "intersectionPoints "
                                    << intersectionPoints << endl;
                                Info << "Face vertex "
                                    << newPoints[newF[npI]] << endl;
                                FatalErrorIn
                                (
                                    "void meshSurfaceCutter::"
                                    "findBoundaryEdgePoints()"
                                ) << "Cannot find the vertex on the edge!!"
                                    << " Cannot close face " << faceI
                                    << abort(FatalError);
                            }
                        }
                        else
                        {
                            //- both vertices are outside the face
                            //- this means that one of them lies on the
                            //- edge of the boundary triangle
                            # ifdef DEBUGCutter
                            Info << "Intersected vertices "
                                << intersectedPoints << endl;
                            Info << "Intersected edges "
                                << intersectedEdges << endl;
                            Info << "Intersections " << edgePoints << endl;
                            Info << "Intersection points "
                                << intersectionPoints << endl;
                            Info << "Intersected neighbours "
                                << intersectedNeighbours << endl;
                            Info << "lastTri " << lastTri << endl;
                            Info << "searchTri " << searchTri << endl;
                            Info << "inwardNext " << inwardNext << endl;
                            Info << "inward " << inward << endl;
                            # endif

                            bool foundDirection(false);

                            forAll(intersectedNeighbours, nI)
                                if(
                                    intersectedNeighbours[nI] == inwardNext
                                )
                                {
                                    searchTri = inwardNext;
                                    foundDirection = true;
                                    break;
                                }
                                else if(
                                    intersectedNeighbours[nI] == inward
                                )
                                {
                                    searchTri = inward;
                                    foundDirection = true;
                                    break;
                                }

                            if( !foundDirection )
                            {
                                const scalar d0 =
                                    mag
                                    (
                                        intersectionPoints[0] -
                                        newPoints[newF[npI]]
                                    );

                                const scalar d1 =
                                    mag
                                    (
                                        intersectionPoints[1] -
                                        newPoints[newF[npI]]
                                    );

                                const scalar el =
                                    mag
                                    (
                                        intersectionPoints[1] -
                                        intersectionPoints[0]
                                    );

                                if( d0 > 0.1*el && d1 > 0.1*el )
                                {
                                    //- intersection vertices are on the
                                    //- edge. Stop serching further.
                                    searchTri = inwardNext;
                                    foundDirection = true;
                                }
                                if( (d1 - d0) > SMALL )
                                {
                                    searchTri = intersectedNeighbours[0];
                                    foundDirection = true;
                                }
                                else if( (d0 - d1) > SMALL )
                                {
                                    searchTri = intersectedNeighbours[1];
                                    foundDirection = true;
                                }
                            }

                            if( !foundDirection )
                            {
                                Info << "Search tri " << searchTri << endl;
                                Info << "intersected neighbours "
                                    << intersectedNeighbours << endl;
                                FatalErrorIn
                                (
                                    "void meshSurfaceCutter::"
                                    "findBoundaryEdgePoints()"
                                ) << "Cannot find the vertex on the edge!!"
                                    << " Cannot close face " << faceI
                                    << abort(FatalError);
                            }
                        }
                    }
                    else if( lastTri == intersectedNeighbours[0] )
                    {
                        lastTri = searchTri;
                        searchTri = intersectedNeighbours[1];
                    }
                    else if( lastTri == intersectedNeighbours[1] )
                    {
                        lastTri = searchTri;
                        searchTri = intersectedNeighbours[0];
                    }
                    else
                    {
                        FatalErrorIn
                         (
                             "void meshSurfaceCutter::findBoundaryEdgePoints()"
                         ) << "Cannot find the next triangle!!"
                             << abort(FatalError);
                    }

                    # ifdef DEBUGCutter
                    Info << "New search tri is " << searchTri << endl;
                    Info << "iPointPatches[0].size() "
                        << iPointPatches[0].size() << endl;
                    Info << "iPointPatches[1].size() "
                        << iPointPatches[1].size() << endl;
                    # endif

                    if( searchTri == intersectedNeighbours[0] )
                    {
                        # ifdef DEBUGCutter
                        Info << "1. Here " << endl;
                        Info << "Point in face "
                            << help::pointInsideFace
                            (intersectionPoints[0], f, normal, meshPoints)
                            << endl;
                        # endif
                        if(
                            iPointPatches[0].size() > 1 &&
                            help::pointInsideFace
                            (intersectionPoints[0], f, normal, meshPoints)
                        )
                        {
                            # ifdef DEBUGCutter
                            Info << "Storing vertex "
                                << intersectionPoints[0]
                                << " with the number " << nPoints_ << endl;
                            # endif
                            //- store the vertex
                            newFaces[faceJ].newElmt(facePointI++) =
                                nPoints_;

                            if( nPoints_ >= pointTriIndex_.size() )
                            {
                                const label oldSize = pointTriIndex_.size();
                                pointTriIndex_.setSize(2*oldSize);
                            }

                            newPoints.append(intersectionPoints[0]);
                            pointTriIndex_[nPoints_] = iPointPatches[0];
                            nPoints_++;
                        }
                    }
                    else if( searchTri == intersectedNeighbours[1] )
                    {
                        # ifdef DEBUGCutter
                        Info << "2. Here " << endl;
                        Info << "Point in face "
                            << help::pointInsideFace
                            (intersectionPoints[1], f, normal, meshPoints)
                            << endl;
                        # endif
                        if(
                            iPointPatches[1].size() > 1 &&
                            help::pointInsideFace
                            (intersectionPoints[1], f, normal, meshPoints)
                        )
                        {
                            # ifdef DEBUGCutter
                            Info << "Storing vertex "
                                << intersectionPoints[1]
                                << " with the number " << nPoints_ << endl;
                            # endif
                            //- store the vertex
                            newFaces[faceJ].newElmt(facePointI++) =
                                nPoints_;
                            newPoints.append(intersectionPoints[1]);
                            if( nPoints_ >= pointTriIndex_.size() )
                            {
                                const label oldSize = pointTriIndex_.size();
                                pointTriIndex_.setSize(2*oldSize);
                            }

                            pointTriIndex_[nPoints_] = iPointPatches[1];
                            nPoints_++;
                        }
                    }

                    if( searchTri == inward )
                    {
                        //- face is closed
                        finished = true;
                        # ifdef DEBUGCutter
                        Info << "facePointI " << facePointI << endl;
                        # endif
                        newFaces[faceJ].setSize(facePointI);
                        faceJ++;
                        newFaces.newElmt(faceJ).setSize(15);
                        facePointI = 0;

                        //- set new inward point
                        inward = inwardNext;
                    }
                    else if( searchTri == inwardNext )
                    {
                        finished = true;
                    }

                    if(
                        tpi.triInfluence() ==
                        trianglePlaneIntersections::NONE
                    )
                    {
                        FatalErrorIn
                        (
                            "void meshSurfaceCutter::findEdgePoints"
                        ) << "Triangle is not intersected by the plane!"
                            << abort(FatalError);
                    }
                }
            }
        }

        if( faceJ == 0 )
        {
            newFaces[faceJ++].setSize(facePointI);
        }

        newFaces.setSize(faceJ);

        # ifdef DEBUGCutter
        Info << "New faces for face " << faceI << " are " << newFaces << endl;
        # endif
    }
    else
    {
        # ifdef DEBUGCutter
        Info << "Storing face for face " << faceI << endl;
        # endif
        facesFromFace_[faceI].setSize(1);
        facesFromFace_[faceI][0] = newF;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
