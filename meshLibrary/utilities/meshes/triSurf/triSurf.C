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

#include "triSurf.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"

#include "gzstream.h"
#include "STLtriangle.H"

#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfaceCleanupDuplicates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triangleIO
{
    // Private data
        //- coordinates
        triangle<Vector<float>, Vector<float> > triangle_;

        //- normal
        Vector<float> normal_;

public:

    // Constructore
        inline triangleIO()
        :
            triangle_
            (
                Vector<float>(),
                Vector<float>(),
                Vector<float>()
            ),
            normal_()
        {}

        inline triangleIO(const triangle<Vector<float>, Vector<float> >& t)
        :
            triangle_(t.a(), t.b(), t.c()),
            normal_(0.5*((b()-a()) ^ (c()-a())))
        {}

    // Member functions

        //- access to a
        inline Vector<float>& a()
        {
            return const_cast<Vector<float>&>(triangle_.a());
        }

        //- access to b
        inline Vector<float>& b()
        {
            return const_cast<Vector<float>&>(triangle_.b());
        }

        //- access to c
        inline Vector<float>& c()
        {
            return const_cast<Vector<float>&>(triangle_.c());
        }

    // Friend operators

        friend inline ostream& operator<<(ostream& os, const triangleIO& t)
        {
            const triangle<Vector<float>, Vector<float> >& tria = t.triangle_;
            os << "  facet normal ";
            os << t.normal_.x() << ' ' << t.normal_.y() << ' ' << t.normal_.z();
            os << nl << "    outer loop";
            os << nl << "      vertex "
               << tria.a().x() << ' ' << tria.a().y() << ' ' << tria.a().z();
            os << nl << "      vertex "
               << tria.b().x() << ' ' << tria.b().y() << ' ' << tria.b().z();
            os << nl << "      vertex "
               << tria.c().x() << ' ' << tria.c().y() << ' ' << tria.c().z();
            os << nl << "    endloop";
            os << "  endfacet";

            return os;
        }

        inline bool read(istream& is)
        {
            char c;

            //- read facet
            if( !((is >> c) && (c == 'f' || c == 'F')) )
            {
                is.putback(c);
                return false;
            }
            is.putback(c);

            //- read facet
            word tok;
            if( !((is >> tok) && (tok == "facet" || tok == "FACET")) )
                return false;

            //- read normal
            if( !((is >> tok) && (tok == "normal" || tok == "NORMAL")) )
                return false;

            //- read coordinates of the normal vector
            is >> normal_.x() >> normal_.y() >> normal_.z();
            if( !is.good() )
                return false;

            //- read outer
            if( !((is >> tok) && (tok == "outer" || tok == "OUTER")) )
                return false;

            //- read loop
            if( !((is >> tok) && (tok == "loop" || tok == "LOOP")) )
                return false;

            //- read vertex
            if( !((is >> tok) && (tok == "vertex" || tok == "VERTEX")) )
                return false;
            is >> a().x() >> a().y() >> a().z();
            if( !is.good() )
                return false;

            //- read vertex
            if( !((is >> tok) && (tok == "vertex" || tok == "VERTEX")) )
                return false;
            is >> b().x() >> b().y() >> b().z();
            if( !is.good() )
                return false;

            //- read vertex
            if( !((is >> tok) && (tok == "vertex" || tok == "VERTEX")) )
                return false;
            is >> this->c().x() >> this->c().y() >> this->c().z();
            if( !is.good() )
                return false;

            //- read endloop
            if( !((is >> tok) && (tok == "endloop" || tok == "ENDLOOP")) )
                return false;

            //- read endfacet
            if( !(is >> tok) )
                return false;

            if( tok != "endfacet" || tok == "ENDFACET" )
                FatalErrorIn
                (
                    "friend inline istream& operator>>(istream&, triangleIO&)"
                ) << "Cannot read triangle" << exit(FatalError);

            return true;
        }
};

void triSurf::readFromSTL(const fileName& fName)
{
    bool compressed = false;

    autoPtr<istream> filePtr(new ifstream(fName.c_str()));

    //- check if the file is compressed or not
    if( !filePtr->good() && isFile(fName + ".gz", false))
    {
        compressed = true;
        filePtr.reset(new igzstream((fName + ".gz").c_str()));
    }
    istream& file = filePtr();

    if( !file.good() )
    {
        FatalErrorIn("void triSurf::readFromSTL(const fileName&)")
            << "Cannot read file " << fName
            << " or file " << fName + ".gz"
            << exit(FatalError);
    }

    //- set the reference to the points
    LongList<point> points;
    LongList<word> patchNames;

    do
    {
        //- read starting token
        word solid;

        if( (file >> solid) && (solid == "solid" || solid == "SOLID") )
        {
            //- read the name of the region
            word solidName;
            if( !(file >> solidName) )
                FatalErrorIn
                (
                    "void triSurf::readFromSTL(const fileName&)"
                ) << "Cannot read solit name" << exit(FatalError);

            triangleIO tri;
            while( tri.read(file) )
            {
                labelledTri ltri;
                ltri[0] = points.size();
                points.append(point(tri.a().x(), tri.a().y(), tri.a().z()));

                ltri[1] = points.size();
                points.append(point(tri.b().x(), tri.b().y(), tri.b().z()));

                ltri[2] = points.size();
                points.append(point(tri.c().x(), tri.c().y(), tri.c().z()));

                ltri.region() = patchNames.size();

                triangles_.append(ltri);
            }

            if( (file >> solid) && (solid == "endsolid" || solid == "ENDSOLID") )
            {
                word endSolid;
                if( (file >> endSolid) && (solidName == endSolid) )
                {
                    patchNames.append(solidName);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void triSurf::readFromSTL(const fileName&)"
                    ) << "Solid name does not end properly"
                      << solidName << exit(FatalError);
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void triSurf::readFromSTL(const fileName&)"
                ) << "Cannot read solid " << solidName << exit(FatalError);
            }
        }
        else
        {
            if( file.eof() )
                break;

            FatalErrorIn
            (
                "void triSurf::readFromSTL(const fileName&)"
            ) << "Cannot read solid " << solid << exit(FatalError);
        }
    } while( file.good() );

    triSurfPoints::points_.setSize(points.size());
    forAll(points, pI)
        triSurfPoints::points_[pI] = points[pI];

    patches_.setSize(patchNames.size());
    forAll(patchNames, patchI)
        patches_[patchI].name() = patchNames[patchI];

    //- merge identical points
    meshOctree octree(*this);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(15, 30);
    triSurfaceCleanupDuplicates(octree).mergeIdentities();

    //- delete addressing which is invalid
    this->clearAddressing();
    this->clearGeometry();
}

void triSurf::writeToSTL(const fileName& fName) const
{
    FatalErrorIn
    (
        "void triSurf::writeToSTL(const fileName& fName) const"
    ) << "Not implemented" << exit(FatalError);
}

void triSurf::readFromSTLB(const fileName& fName)
{
    bool compressed = false;

    autoPtr<istream> filePtr(new ifstream(fName.c_str(), std::ios::binary));

    // If the file is compressed, decompress it before reading.
    if( !filePtr->good() && isFile(fName + ".gz", false))
    {
        compressed = true;
        filePtr.reset(new igzstream((fName + ".gz").c_str()));
    }
    istream& file = filePtr();

    if( !file.good() )
    {
        FatalErrorIn("void triSurf::readFromSTLB(const fileName&)")
            << "Cannot read file " << fName
            << " or file " << fName + ".gz"
            << exit(FatalError);
    }

    //- check the header
    const label STLheaderSize(80);
    char header[STLheaderSize];
    file.read(header, STLheaderSize);

    // Check that stream is OK, if not this maybe an ASCII file
    if( !file )
        return;

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    file.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    if( !file || nTris < 0 )
        return;

    // Compare the size of the file with that expected from the number of tris
    // If the comparison is not sensible then it maybe an ASCII file
    if( !compressed )
    {
        const label dataFileSize = Foam::fileSize(fName) - 80;

        if( nTris < dataFileSize/50 || nTris > dataFileSize/25 )
            return;
    }

    //- set size of points
    pointField& points = triSurfPoints::points_;
    points.setSize(3*nTris);

    // Allocate storage for triangles
    triangles_.setSize(nTris);

    label pointI(0);

    forAll(triangles_, triI)
    {
        labelledTri& tri = triangles_[triI];

        //- read the triangles from the stream
        STLtriangle stlTri(file);

        // Set the rawPoints to the vertices of the STL triangle
        // and set the point labels of the labelledTri
        points[pointI] = stlTri.a();
        tri[0] = pointI++;

        points[pointI] = stlTri.b();
        tri[1] = pointI++;

        points[pointI] = stlTri.c();
        tri[2] = pointI++;

        tri.region() = stlTri.region();
    }

    //- merge identical points
    meshOctree octree(*this);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(15, 30);
    triSurfaceCleanupDuplicates(octree).mergeIdentities();

    //- delete addressing which is invalid
    this->clearAddressing();
    this->clearGeometry();
}

void triSurf::writeToSTLB(const fileName& fName) const
{
    FatalErrorIn
    (
        "void triSurf::writeToSTLB(const fileName& fName) const"
    ) << "Not implemented" << exit(FatalError);
}

void triSurf::readFromFTR(const fileName& fName)
{
    IFstream fStream(fName);

    fStream >> triSurfFacets::patches_;

    fStream >> triSurfPoints::points_;

    fStream >> triSurfFacets::triangles_;

    token c;
    fStream >> c;
    if( fStream.eof() )
    {
        return;
    }
    else
    {
        fStream.putBack(c);
    }

    List<meshSubset> subsets;

    //- read point subsets
    fStream >> subsets;
    forAll(subsets, subsetI)
        triSurfPoints::pointSubsets_.insert(subsetI, subsets[subsetI]);

    subsets.clear();

    //- read facet subsets
    fStream >> subsets;
    forAll(subsets, subsetI)
        triSurfFacets::facetSubsets_.insert(subsetI, subsets[subsetI]);
}

void triSurf::writeToFTR(const fileName& fName) const
{
    OFstream fStream(fName);

    fStream << triSurfFacets::patches_;

    fStream << nl;

    fStream << triSurfPoints::points_;

    fStream << nl;

    fStream << triSurfFacets::triangles_;

    fStream << nl;

    fStream << triSurfPoints::pointSubsets_;

    fStream << nl;

    fStream << triSurfFacets::facetSubsets_;
}

/*
void triSurf::calculateFaceGroups() const
{
    nFaceGroups_ = 0;

    faceGroupPtr_ = new labelListPMG(triSurface::size(), -1);
    labelListPMG& faceGroup = *faceGroupPtr_;

    const labelListList& faceEdges = this->faceEdges();
    const labelListList& edgeFaces = this->edgeFaces();

    labelListPMG front;

    forAll(faceGroup, fI)
    {
        if( faceGroup[fI] != -1 )
            continue;

        front.clear();
        front.append(fI);
        faceGroup[fI] = nFaceGroups_;

        while( front.size() != 0 )
        {
            const label fLabel = front.removeLastElement();

            const labelList& fEdges = faceEdges[fLabel];
            forAll(fEdges, feI)
            {
                const label eI = fEdges[feI];

                if( edgeFaces[eI].size() != 2 )
                    continue;

                label nei = edgeFaces[eI][0];
                if( nei == fLabel )
                    nei = edgeFaces[eI][1];

                if( faceGroup[nei] == -1 )
                {
                    faceGroup[nei] = nFaceGroups_;
                    front.append(nei);
                }
            }
        }

        ++nFaceGroups_;
    }
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurf::triSurf()
:
    triSurfPoints(),
    triSurfFacets(),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{}

//- Construct from parts
triSurf::triSurf
(
    const LongList<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& points
)
:
    triSurfPoints(points),
    triSurfFacets(triangles, patches),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{}

//- Read from file
triSurf::triSurf(const fileName& fName)
:
    triSurfPoints(),
    triSurfFacets(),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{
    readSurface(fName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurf::~triSurf()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurf::readSurface(const fileName& fName)
{
    if( fName.ext() == "ftr" || fName.ext() == "FTR")
    {
        readFromFTR(fName);
    }
    else if( fName.ext() == "stl" || fName.ext() == "STL" )
    {
        readFromSTL(fName);
    }
    else if( fName.ext() == "stlb" || fName.ext() == "STLB" )
    {
        readFromSTLB(fName);
    }
    else
    {
        WarningIn("void triSurf::readSurface(const fileName& fName)")
            << "Unknown surface format for file " << fName << endl;
    }
}

void triSurf::writeSurface(const fileName& fName) const
{
    if( fName.ext() == "ftr" || fName.ext() == "FTR" )
    {
        writeToFTR(fName);
    }
    else if( fName.ext() == "stl" || fName.ext() == "STL" )
    {
        writeToSTL(fName);
    }
    else if( fName.ext() == "stlb" || fName.ext() == "STLB" )
    {
        writeToSTLB(fName);
    }
    else
    {
        WarningIn("void triSurf::readSurface(const fileName& fName)")
            << "Unknown surface format for file " << fName << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
