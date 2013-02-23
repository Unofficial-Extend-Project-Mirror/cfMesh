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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void triSurf::existingFaceSubsets(DynList<word>& names) const
{
    names.clear();

    std::map<word, labelListPMG>::const_iterator iter;
    for(iter=triSubsets_.begin();iter!=triSubsets_.end();++iter)
        names.append(iter->first);
}

inline bool triSurf::doesFaceSubsetExist(const word& name) const
{
    if( triSubsets_.find(name) != triSubsets_.end() )
        return true;

    return false;
}

inline const labelListPMG& triSurf::facesInSubset(const word& name) const
{
    if( !doesFaceSubsetExist(name) )
        FatalErrorIn
        (
            "inline const labelListPMG& triSurf::"
            "facesInSubset(const word& name) const"
        ) << "Subset " << name << " does not exist!" << exit(FatalError);

    std::map<word, labelListPMG>::const_iterator iter = triSubsets_.find(name);

    return iter->second;
}

inline void triSurf::addFacetsToSubset
(
    const word& name,
    const labelListPMG& newElements
)
{
    if( doesFaceSubsetExist(name) )
    {
        std::map<word, labelListPMG>::iterator iter = triSubsets_.find(name);

        labelListPMG& elmts = iter->second;
        boolList alreadyInside(this->size(), false);
        forAll(elmts, elI)
            alreadyInside[elmts[elI]] = true;

        forAll(newElements, elI)
            if(
                (newElements[elI] >= 0) &&
                (newElements[elI] < this->size()) &&
                !alreadyInside[newElements[elI]]
            )
            {
                elmts.append(newElements[elI]);
                alreadyInside[newElements[elI]] = true;
            }
    }
    else
    {
        triSubsets_.insert(std::pair<word, labelListPMG>(name, newElements));
    }
}

inline void triSurf::removeFaceSubset(const word& name)
{
    if( doesFaceSubsetExist(name) )
        triSubsets_.erase(name);
}

inline label triSurf::numberOfFaceGroups() const
{
    if( !faceGroupPtr_ )
        calculateFaceGroups();

    return nFaceGroups_;
}

inline const labelListPMG& triSurf::faceGroups() const
{
    if( !faceGroupPtr_ )
        calculateFaceGroups();

    return *faceGroupPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
