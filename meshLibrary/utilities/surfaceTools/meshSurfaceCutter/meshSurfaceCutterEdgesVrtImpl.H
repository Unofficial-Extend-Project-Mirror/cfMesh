// Info << "Face intersects a vertex. This is tricky" << endl;
    const point& p = points[faces[searchTri][eI]];
    const constRow pf = surface_.pointFacets()[faces[searchTri][eI]];

    DynList<label> patches;
    patches.append(searchTri);

    //- check if the vertex is at a boundary edge

    forAll(pf, pfI)
        if( pf[pfI] != searchTri )
        {
            bool store(true);

            forAll(patches, patchI)
                if(
                    faces[pf[pfI]].region() ==
                    faces[patches[patchI]].region()
                )
                    store = false;

            if( store )
                patches.append(pf[pfI]);
        }

    bool foundTri(false);
    //- examine if any of the neighbouring triangles is a lastTri inward
    forAll(pf, pfI)
        if( pf[pfI] == lastTri )
        {
            intersectedNeighbours[counter] = lastTri;
            intersectionPoints[counter] = p;
            foundTri = true;
        }
        else if( pf[pfI] == inwardNext )
        {
            intersectedNeighbours[counter] = inwardNext;
            intersectionPoints[counter] = p;
            foundTri = true;
        }
        else if( pf[pfI] == inward )
        {
            intersectedNeighbours[counter] = inward;
            intersectionPoints[counter] = p;
            foundTri = true;
        }

    if( !foundTri )
        forAll(pf, pfI)
            if( pf[pfI] != searchTri )
            {
                trianglePlaneIntersections tpj
                (
                    normal,
                    meshPoints[f[0]],
                    surface_,
                    pf[pfI]
                );

                if(
                    tpj.triInfluence() ==
                    trianglePlaneIntersections::VERTEX_AND_EDGE
                )
                {
                    intersectedNeighbours[counter] = pf[pfI];
                    intersectionPoints[counter] = p;
                    break;
                }
                else if(
                    tpj.triInfluence() ==
                    trianglePlaneIntersections::TWO_VERTICES
                )
                {
                    intersectedNeighbours[counter] = pf[pfI];
                    intersectionPoints[counter] = p;

                    //- check that the intersected vertices are not the same
                    //- vertices as in the searchTri
                    bool neighbour(true);
                    if(
                        tpi.triInfluence() ==
                        trianglePlaneIntersections::TWO_VERTICES
                    )
                    {

                        forAll(intersectedPoints, ipI)
                        {
                            forAll(tpj.intersectedPoints(), ipJ)
                                if(
                                    intersectedPoints[ipI] &&
                                    tpj.intersectedPoints()[ipJ] &&
                                    (
                                        faces[searchTri][ipI] !=
                                        faces[pf[pfI]][ipJ]
                                    )
                                )
                                {
                                    neighbour = false;
                                    break;
                                }

                            if( !neighbour ) break;
                        }
                    }

                    //- finish searching if the triangle is not a neighbour
                    if( !neighbour ) break;
                }
            }

    //- store patches
    iPointPatches[counter++] = patches;
