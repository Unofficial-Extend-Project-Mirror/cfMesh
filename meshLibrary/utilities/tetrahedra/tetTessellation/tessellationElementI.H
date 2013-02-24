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

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void tessellationElement::setNode(const direction i, const label nodeI)
{
    data_[i] = nodeI;
}

inline label tessellationElement::neighbour(const direction i) const
{
    return neighbours_[i];
}

inline void tessellationElement::setNeighbour
(
    const direction i,
    const label neiI
)
{
    neighbours_[i] = neiI;
}

inline triFace tessellationElement::face(const direction i) const
{
    switch( i )
    {
        case 0:
        {
            return triFace(data_[1], data_[3], data_[2]);
        } break;
        case 1:
        {
            return triFace(data_[0], data_[2], data_[3]);
        } break;
        case 2:
        {
            return triFace(data_[0], data_[3], data_[1]);
        } break;
        case 3:
        {
            return triFace(data_[0], data_[1], data_[2]);
        } break;
    };
    
    return triFace();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void tessellationElement::operator=(const tessellationElement& elmt)
{
    for(label i=0;i<4;++i)
    {
        data_[i] = elmt[i];
        neighbours_[i] = elmt.neighbour(i);
    }
    
    influence_ = elmt.influence_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
