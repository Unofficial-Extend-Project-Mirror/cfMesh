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
