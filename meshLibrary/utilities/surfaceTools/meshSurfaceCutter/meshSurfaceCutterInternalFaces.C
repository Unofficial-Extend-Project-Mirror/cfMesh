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

#include "meshSurfaceCutter.H"
#include "surfaceIntersectionsOctree.H"
#include "polyMeshGenAddressing.H"

#define DEBUGCutter

# ifdef DEBUGCutter
#include "triSurfaceSearch.H"
#include "octreeDataTriSurface.H"
#include "octree.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
void meshSurfaceCutter::findInternalPointsAndIntersections()
{
	//- create poly edges
	const polyMeshGenAddressing& addressingData = mesh_.addressingData();
	const edgeList& edges = addressingData.edges();
	const VRWGraph& faceEdges = addressingData.faceEdges();

    # ifdef DEBUGCutter
    Info << "Edges " << edges << endl;
    Info << "faceEdges " << faceEdges << endl;
    # endif

	//- create the octree
	surfaceIntersectionsOctree so(surface_, 100, 10);
	
	const pointFieldPMG& points = mesh_.points();

	if( !newPointsPtr_ )
		newPointsPtr_ = new DynList<point>();
	DynList<point>& newPoints = *newPointsPtr_;
	newPoints.setSize(points.size());
	
	if( !newPointLabelPtr_ )
		newPointLabelPtr_ = new labelList();
    labelList& newPointLabels = *newPointLabelPtr_;
	newPointLabels.setSize(points.size());
	newPointLabels = -1;
	
    const label np = points.size();
    nPoints_ = 0;
    Serr << "Creating internal vertices:";
    label index(0);

    forAll(points, pI)
    {
        if( so.isPointInside(points[pI]) )
        {
            # ifdef DEBUGCutter
            Info << "Point " << pI << "  " << points[pI]
                << " is inside. New label is " << nPoints_ << endl;
            # endif
            newPointLabels[pI] = nPoints_;
            newPoints.append(points[pI]);
            nPoints_++;
        }
            
        if( index++ >= np / 10 )
        {
            Serr << ".";
            index = 0;
        }
    }
    
    Info << endl;

    # ifdef DEBUGCutter
    Info << "Internal points " << newPoints << endl;
    # endif

    Serr << "Number of internal points " << nPoints_ << endl;

    //- some edges of the Voronoi diagram intersect the boundary surface
    //- it is therefore necessary to find those intersections as
    //- they will be the boundary points
    //- Every edge of the Voronoi diagram appears in three faces so it is
    //- necessary to store the information about these points
    //- into a list (edgePoint)
    //- pointTriIndex list holds information about the index of the intersected
    //- boundary triangle for every mesh point
    Serr << "Cutting edges:";

    pointTriIndex_.setSize(points.size());

	if( !edgePointPtr_ )
		edgePointPtr_ = new edgeList();
    edgeList& edgePoint = *edgePointPtr_;
	edgePoint.setSize(edges.size());
	edgePoint = edge(-1, -1);
	//Info << "Number of edges is " << edges.size() << endl;

    index = 0;

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        if(
            newPointLabels[e.start()] != -1 &&
            newPointLabels[e.end()] == -1
        )
        {
            pointIndexHit ss =
                so.intersection(points[e.start()], points[e.end()]);

            if( ss.hit() )
            {
                edgePoint[eI].start() = nPoints_;
                edgePoint[eI].end() = -1;
				
				# ifdef DEBUGCutter
				Info << nl << "1. Points for edge " << eI
					<< " are " << edges[eI] << endl;
				Info << "Start vertex " << points[e.start()] << endl;
				Info << "End vertex " << points[e.end()] << endl;
				Info << "Intersection " << ss.rawPoint() << endl;
				# endif

                newPoints.append(ss.rawPoint());
                pointTriIndex_.newElmt(nPoints_).append(ss.index());
                nPoints_++;
            }
            else
            {
                Info << "Vertex " << e.start() << "  "
                    << points[e.start()]
                    << " is inside" << endl;
                Info << "Vertex " << e.end() << "  " << points[e.end()]
                    << " is outside" << endl;
                FatalErrorIn
                (
                    "void createInternalMeshPointsAndFaces()"
                ) << "Intersection with the surface is not found"
                    << abort(FatalError);
            }
        }
        else if(
            newPointLabels[e.start()] == -1 &&
            newPointLabels[e.end()] != -1
        )
        {
            pointIndexHit ss =
                so.intersection(points[e.end()], points[e.start()]);

            if( ss.hit() )
            {
                edgePoint[eI].start() = -1;
                edgePoint[eI].end() = nPoints_;
				
				# ifdef DEBUGCutter
				Info << nl << "2. Points for edge " << eI
					<< " are " << edges[eI] << endl;
				Info << "Start vertex " << points[e.start()] << endl;
				Info << "End vertex " << points[e.end()] << endl;
				Info << "Intersection " << ss.rawPoint() << endl;
				# endif

                newPoints.append(ss.rawPoint());
                pointTriIndex_.newElmt(nPoints_).append(ss.index());
                nPoints_++;
            }
            else
            {
                Info << "Vertex " << e.end() << "  "
                    << points[e.end()]
                    << " is inside" << endl;
                Info << "Vertex " << e.start() << "  "
                    << points[e.start()] << " is outside" << endl;
                FatalErrorIn
                (
                    "void createInternalMeshPointsAndFaces()"
                ) << "Intersection with the surface is not found"
                    << abort(FatalError);
            }
        }
        else
        {
            pointIndexHit ss, ss1;

            so.intersection
            (
                points[e.start()],
                points[e.end()],
                ss,
                ss1
            );

            if( ss.hit() && ss1.hit() && (ss.index() != ss1.index()) )
            {
                edgePoint[eI].start() = nPoints_;
                newPoints.append(ss.rawPoint());

				# ifdef DEBUGCutter
				Info << nl << "3. Points for edge " << eI
					<< " are " << edges[eI] << endl;
				Info << "Start vertex " << points[e.start()] << endl;
				Info << "End vertex " << points[e.end()] << endl;
				Info << "Intersection1 " << ss.rawPoint() << endl;
				Info << "Intersection2 " << ss1.rawPoint() << endl;
				# endif
                pointTriIndex_.newElmt(nPoints_).append(ss.index());
                nPoints_++;

                edgePoint[eI].end() = nPoints_;
                newPoints.append(ss1.rawPoint());

                pointTriIndex_.newElmt(nPoints_).append(ss1.index());
                nPoints_++;
            }
            else if( ss.hit() || ss1.hit() )
            {
				Info << "ss " << ss << endl;
				Info << "ss1 " << ss1 << endl;
                FatalErrorIn
                (
                    "void meshSurfaceCutter::createInternalMeshPointsAndFaces()"
                ) << "I expect 2 or zero intersections! Aborting"
                    << exit(FatalError);
            }
            else
            {
                edgePoint[eI].start() = -1;
                edgePoint[eI].end() = -1;
            }
        }

        if( index++ >= 0.1 * edges.size() )
        {
            Serr << ".";
            index = 0;
        }
    }

    Info << endl;

    # ifdef DEBUGCutter
    forAll(edgePoint, eI)
    {
        const edge& e = edgePoint[eI];
        if(
            e.start() != -1 &&
            (
                pointTriIndex_[e.start()].size() > 1 ||
                pointTriIndex_[e.start()][0] < 0
            )
        )
            FatalErrorIn("void meshSurfaceCutter::createInternalMeshPointsAndFaces()")
                << "Edge " << eI << " has got incorrect intersections"
                    << abort(FatalError);

        if(
            e.end() != -1 &&
            (
                pointTriIndex_[e.end()].size() > 1 ||
                pointTriIndex_[e.end()][0] < 0
            )
        )
            FatalErrorIn("void meshSurfaceCutter::createInternalMeshPointsAndFaces()")
                << "Edge " << eI << " has got incorrect intersections"
                    << abort(FatalError);
    }

    Info << "2. New points with boundary points " << newPoints << endl;
    Info << "Number of points before edge points are created are "
        << nPoints_ << endl;

    Info << "edge patch triangles " << endl;
    forAll(edgePoint, eI)
    {
        Info << "edgePoint["<<eI<<"] are "<< edgePoint[eI] << endl;
    }
    # endif
	
}

void meshSurfaceCutter::createInternalMeshPointsAndFaces()
{
	findInternalPointsAndIntersections();

    const labelList owner = mesh_.owner();
    const labelList neighbour = mesh_.neighbour();
	const faceListPMG& faces = mesh_.faces();

    //- Finally, create new internal faces
    Serr << "Creating internal faces:";
    if( !facesFromFacePtr_ )
    {
        facesFromFacePtr_ =
            new faceListList(neighbour.size());
    }
    else
    {
        FatalErrorIn
        (
            "void meshSurfaceCutter::createInternalMeshPointsAndFaces()"
        ) << "facesFromFacePtr_ ia already allocated!!" << abort(FatalError);
    }

    faceListList& facesFromFace_ = *facesFromFacePtr_;

    label index(0);

    forAll(facesFromFace_, faceI)
    {
        if( neighbour[faceI] != -1 )
        {
            # ifdef DEBUGCutter
            Info << "Cutting face " << faceI << " for cell "
                << owner[faceI] << " and cell " << neighbour[faceI] << endl;
			Info << "Original face " << faces[faceI] << endl;
            # endif

            //- create a mesh face from a face of a Voronoi polyhedron
            if( faceCutter(faceI) )
            {
                boundaryCell_[owner[faceI]] = true;
                boundaryCell_[neighbour[faceI]] = true;
            }
			
			# ifdef DEBUGCutter
			Info << "New faces for face are " << facesFromFace_[faceI] << endl;
			# endif
			
            if( index++ >= neighbour.size() / 10 )
            {
                Serr << ".";
                index = 0;
            }
        }
    }

    Info << endl;

    //- delete surface topology addressing
    const_cast<triSurface&>(surface_).clearTopology();

    //- store internal faces
    forAll(nFacesInCell_, cI)
    {
        nFacesInCell_[cI] = 0;
    }
	
    nIntFaces_ = 0;
	
	polyMeshGenModifier meshModifier(mesh_);
	cellListPMG& newCells = meshModifier.cellsAccess();
	faceListPMG& newFaces = meshModifier.facesAccess();

    forAll(facesFromFace_, faceI)
    {
        # ifdef DEBUGCutter
        Info << "Original face " << faces[faceI] << endl;
        forAll(faces[faceI], pI)
            Info << pI << " " << (*newPointLabelPtr_)[faces[faceI][pI]] << ",";
        Info << endl;
        Info << "Faces for face " << faceI << " are "
            << facesFromFace_[faceI] << endl;
        # endif

        //- if the face is cut into two mesh faces then it may happen
        //- that the owner and/or neighbour cells consist of two parts
        if( facesFromFace_[faceI].size() > 1 )
        {
			# ifdef DEBUCutter
            Info << "face " << faceI << " " << faces[faceI]
                << nl << "owner " << owner[faceI]
                << nl << "neighbour " << neighbour[faceI]
                << nl << " is decomposed into "
                << facesFromFace_[faceI] << endl;
			# endif
            problematicTopology_[owner[faceI]] = true;
            problematicTopology_[neighbour[faceI]] = true;
        }

        forAll(facesFromFace_[faceI], fJ)
        {
            const face& f = facesFromFace_[faceI][fJ];
            # ifdef DEBUGCutter
            Info << "Storing face " << nIntFaces_ << "  " << f 
                << " into cells " << owner[faceI]
                << " and " << neighbour[faceI] << endl;
            # endif
            if( f.size() > 2 )
			{
				//- store cut internal faces
				const label own = owner[faceI];
				const label nei = neighbour[faceI];
				newCells[own][nFacesInCell_[own]++] = nIntFaces_;
				newCells[nei][nFacesInCell_[nei]++] = nIntFaces_;
				newFaces.newElmt(nIntFaces_++) = f;
			}
                
        }
    }

    Info << "Number of internal faces is " << nIntFaces_ << endl;

    //- facesFromFace_ are not needed any more
    deleteDemandDrivenData(facesFromFacePtr_);

    # ifdef DEBUGCutter
    Info << "Number of points after internal faces are created is "
        << nPoints_ << endl;

    for(label i=0;i<nIntFaces_;i++)
    {
        Info << "New face " << i << " is " << faces[i] << endl;

        forAll(faces[i], pI)
            Info << "Face point " << pI <<  " has coordinates "
                << (*newPointsPtr_)[faces[i][pI]] << endl;
    }
    # endif

    //- store new points into the mesh
    for(index=0;index<nPoints_;++index)
        meshModifier.pointsAccess().newElmt(index) = (*newPointsPtr_)[index];

    deleteDemandDrivenData(newPointsPtr_);

    # ifdef DEBUGCutter
    //- check cell topology
    //- every internal edge should appear twice in one cell
    //- the boundary ones should appear once
    forAll(newCells, cellI)
        newCells[cellI].setSize(nFacesInCell_[cellI]);

    meshModifier.pointsAccess().setSize(nPoints_);

    newFaces.setSize(nIntFaces_);
    Info << "Number of polyhedral cells is " << newCells.size() << endl;

    forAll(newCells, cellI)
    {
        Info << "Checking closedness for cell " << cellI
            << " consisting of newFaces " << newCells[cellI] << endl;
        
        const edgeList edges = newCells[cellI].edges(newFaces);
        labelList nAppearances(edges.size(), 0);

        Info << "Cell edges are " << edges << endl;

        const cell& c = newCells[cellI];

        forAll(c, fI)
        {
            Info << "Face " << c[fI] << " is " << newFaces[c[fI]] << endl;
            const edgeList fEdges = newFaces[c[fI]].edges();

            forAll(fEdges, eI)
                for(short eJ=0;eJ<edges.size();eJ++)
                    if( edges[eJ] == fEdges[eI] )
                    {
                        nAppearances[eJ]++;
                        break;
                    }
        }

        forAll(edges, eI)
            if( nAppearances[eI] > 2 )
            {
                FatalErrorIn
                (
                    "void createInternalMeshPointsAndFaces()"
                ) << "Edge " << edges[eI]
                    << " appears more than two times in cell " << cellI
                    << abort(FatalError);
            }
            else if
            (
                (nAppearances[eI] == 1) &&
                (
                    !pointTriIndex_[edges[eI].start()].size() ||
                    !pointTriIndex_[edges[eI].end()].size()
                )
            )
            {
                Info << nl << "Point " << edges[eI].start() << " is "
                    << meshModifier.pointsAccess()[edges[eI].start()]
                    << " tri index " << pointTriIndex_[edges[eI].start()]
                    << endl;
                Info << "Point " << edges[eI].end() << " is "
                    << meshModifier.pointsAccess()[edges[eI].end()]
                    << " tri index " << pointTriIndex_[edges[eI].end()]
                    << endl;
                FatalErrorIn
                (
                    "void createInternalMeshPointsAndFaces()"
                ) << "Edge " << edges[eI]
                    << " appears once in cell " << cellI
                    << " and it is not on the boundary" << abort(FatalError);
            } 
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
