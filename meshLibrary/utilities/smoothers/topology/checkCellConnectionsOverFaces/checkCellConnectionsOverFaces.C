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

#include "checkCellConnectionsOverFaces.H"
#include "labelListPMG.H"
#include "labelledPair.H"
#include "helperFunctionsPar.H"

#include <omp.h>
#include <set>
#include <map>

#include "helperFunctions.H"
#include "writeMeshFPMA.H"

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

typedef std::pair<label, label> lPair;
typedef std::pair<lPair, lPair> lPairPair;
typedef std::pair<lPair, label> lPairLabel;

inline Ostream& operator<<(Ostream& os, const lPair& lp)
{
    os << token::BEGIN_LIST;
    os << lp.first;
    os << token::SPACE;
    os << lp.second;
    os << token::END_LIST;

    return os;
}

inline Istream& operator>>(Istream& is, lPair& lp)
{
    is.readBegin("lPair");

    is >> lp.first;
    is >> lp.second;

    is.readEnd("lPair");
    is.check("operator>>(Istream&, lPair&");

    return is;
}

inline Ostream& operator<<(Ostream& os, const lPairPair& lp)
{
    os << token::BEGIN_LIST;
    os << lp.first;
    os << token::SPACE;
    os << lp.second;
    os << token::END_LIST;

    return os;
}

inline Istream& operator>>(Istream& is, lPairPair& lp)
{
    is.readBegin("lPairPair");

    is >> lp.first;
    is >> lp.second;

    is.readEnd("lPairPair");
    is.check("operator>>(Istream&, lPairPair&");

    return is;
}

inline Ostream& operator<<(Ostream& os, const lPairLabel& lp)
{
    os << token::BEGIN_LIST;
    os << lp.first;
    os << token::SPACE;
    os << lp.second;
    os << token::END_LIST;

    return os;
}

inline Istream& operator>>(Istream& is, lPairLabel& lp)
{
    is.readBegin("lPairLabel");

    is >> lp.first;
    is >> lp.second;

    is.readEnd("lPairLabel");
    is.check("operator>>(Istream&, lPairLabel&");

    return is;
}

template<>
inline bool contiguous<lPairPair>() {return true;}

template<>
inline bool contiguous<lPairLabel>() {return true;}

class groupingOp
{
    // Private data
        //- map containing group connections
        std::map<lPair, std::set<lPair> >& groupMap_;

    // Private member functions
        //- sort neighbouring groups
        void sortNeighbouringGroups()
        {
            typedef std::set<lPair> pairSet;
            typedef std::map<lPair, pairSet> gMap;

            DynList<lPair> eraseRow;

            for
            (
                gMap::reverse_iterator it=groupMap_.rbegin();
                it!=groupMap_.rend();
                ++it
            )
            {
                const pairSet& neiGroups = it->second;

                gMap::reverse_iterator rIt = it;
                for(++rIt;rIt!=groupMap_.rend();++rIt)
                {
                    pairSet& otherGroups = rIt->second;

                    if( otherGroups.find(it->first) != otherGroups.end() )
                    {
                        //- add all neighbours of current group to
                        //- otherGroups
                        for
                        (
                            pairSet::const_iterator sIter=neiGroups.begin();
                            sIter!=neiGroups.end();
                            ++sIter
                        )
                            otherGroups.insert(*sIter);

                        eraseRow.append(it->first);
                    }
                    else
                    {
                        for
                        (
                            pairSet::const_iterator sIter=neiGroups.begin();
                            sIter!=neiGroups.end();
                            ++sIter
                        )
                        {
                            const lPair& lp = *sIter;

                            if( otherGroups.find(lp) != otherGroups.end() )
                            {
                                eraseRow.append(it->first);
                                otherGroups.insert(it->first);

                                pairSet::const_iterator sIt;
                                for
                                (
                                    sIt=neiGroups.begin();
                                    sIt!=neiGroups.end();
                                    ++sIt
                                )
                                    rIt->second.insert(*sIt);

                                break;
                            }
                        }
                    }
                }
            }

            forAll(eraseRow, rowI)
                groupMap_.erase(eraseRow[rowI]);
        }

    public:

    // Constructor
        //- construct from map
        groupingOp(std::map<lPair, std::set<lPair> >& groupMap)
        :
            groupMap_(groupMap)
        {
            sortNeighbouringGroups();
        }

    // Public member functions
        //- prepare information for sending to other processors
        void dataForNeighbourProc
        (
            DynList<lPairLabel>& groupLabels,
            const label procNo,
            const labelList& globalGroupLabel
        ) const
        {
            groupLabels.clear();

            typedef std::map<lPair, std::set<lPair> > gMap;
            typedef std::set<lPair> pairSet;

            pairSet alreadyAdded;

            forAllConstIter(gMap, groupMap_, it)
            {
                if( it->first.first != Pstream::myProcNo() )
                    continue;

                const pairSet& neiGroups = it->second;

                forAllConstIter(pairSet, neiGroups, sIt)
                {
                    if( alreadyAdded.find(*sIt) != alreadyAdded.end() )
                        continue;

                    groupLabels.append
                    (
                        std::make_pair
                        (
                            *sIt,
                            globalGroupLabel[it->first.second]
                        )
                    );

                    alreadyAdded.insert(*sIt);
                }
            }
        }

        template<class ListType>
        bool updateGlobalGroupLabels
        (
            labelList& globalGroupLabel,
            const ListType& receivedGroups
        ) const
        {
            bool changed(false);

            forAll(receivedGroups, i)
            {
                const lPairLabel& lpl = receivedGroups[i];

                if( lpl.first.first != Pstream::myProcNo() )
                    continue;

                if( lpl.second < globalGroupLabel[lpl.first.second] )
                {
                    globalGroupLabel[lpl.first.second] = lpl.second;
                    changed = true;
                }
            }

            if( updateGlobalGroupLabels(globalGroupLabel) )
                changed = true;

            return changed;
        }

        bool updateGlobalGroupLabels(labelList& globalGroupLabel) const
        {
            typedef std::map<lPair, std::set<lPair> > gMap;
            typedef std::set<lPair> pairSet;

            bool changed(false);

            forAllConstIter(gMap, groupMap_, it)
            {
                if( it->first.first != Pstream::myProcNo() )
                    continue;

                const label groupI = globalGroupLabel[it->first.second];

                const pairSet& neiGroups = it->second;
                forAllConstIter(pairSet, neiGroups, sIt)
                {
                    if( sIt->first != Pstream::myProcNo() )
                        continue;

                    if( globalGroupLabel[sIt->second] < groupI )
                    {
                        forAllConstIter(pairSet, neiGroups, sIt2)
                        {
                            if( sIt2->first == Pstream::myProcNo() )
                                globalGroupLabel[sIt2->second] =
                                    globalGroupLabel[sIt->second];
                        }

                        changed = true;
                    }
                    else if( globalGroupLabel[sIt->second] > groupI )
                    {
                        forAllConstIter(pairSet, neiGroups, sIt2)
                        {
                            if( sIt2->first == Pstream::myProcNo() )
                                globalGroupLabel[sIt2->second] = groupI;
                        }

                        changed = true;
                    }
                }
            }

            return changed;
        }

        bool exchangeDataTopToBottom(labelList& globalGroupLabel) const
        {
            bool changed(false);

            const Pstream::commsStruct& myComm =
                Pstream::treeCommunication()[Pstream::myProcNo()];

            //- propagate data to the processors below
            LongList<lPairLabel> propagateData;

            if( myComm.above() != -1 )
            {
                List<lPairLabel> receivedGroupLabels;
                IPstream fromOtherProc
                (
                    Pstream::scheduled,
                    myComm.above()
                );

                fromOtherProc >> receivedGroupLabels;

                forAll(receivedGroupLabels, i)
                {
                    if(
                        receivedGroupLabels[i].first.first <=
                        Pstream::myProcNo()
                    )
                        continue;

                    propagateData.append(receivedGroupLabels[i]);
                }

                bool check =
                    updateGlobalGroupLabels
                    (
                        globalGroupLabel,
                        receivedGroupLabels
                    );

                if( check )
                    changed = true;
            }

            forAll(myComm.below(), belowI)
            {
                //- send the data to the processors below
                DynList<lPairLabel> sendGroupLabels;
                dataForNeighbourProc
                (
                    sendGroupLabels,
                    myComm.below()[belowI],
                    globalGroupLabel
                );

                forAll(propagateData, i)
                    sendGroupLabels.append(propagateData[i]);

                OPstream toOtherProc
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    sendGroupLabels.byteSize()
                );

                toOtherProc << sendGroupLabels;
            }

            reduce(changed, maxOp<bool>());

            return changed;
        }

        bool exchangeDataBottomToTop(labelList& globalGroupLabel) const
        {
            bool changed(false);

            const Pstream::commsStruct& myComm =
                Pstream::treeCommunication()[Pstream::myProcNo()];

            //- propagate data to the processors above
            LongList<lPairLabel> propagateData;
            forAll(myComm.below(), belowI)
            {
                IPstream fromOtherProc
                (
                    Pstream::scheduled,
                    myComm.below()[belowI]
                );

                List<lPairLabel> receivedGroups;
                fromOtherProc >> receivedGroups;

                forAll(receivedGroups, i)
                {
                    if( receivedGroups[i].first.first >= Pstream::myProcNo() )
                        continue;

                    propagateData.append(receivedGroups[i]);
                }

                bool check =
                    updateGlobalGroupLabels
                    (
                        globalGroupLabel,
                        receivedGroups
                    );
                if( check )
                    changed = true;
            }

            if( myComm.above() != -1 )
            {
                DynList<lPairLabel> sendData;
                dataForNeighbourProc
                (
                    sendData,
                    myComm.above(),
                    globalGroupLabel
                );

                forAll(propagateData, i)
                    sendData.append(propagateData[i]);

                OPstream toOtherProc
                (
                    Pstream::scheduled,
                    myComm.above(),
                    sendData.byteSize()
                );

                toOtherProc << sendData;
            }

            reduce(changed, maxOp<bool>());

            return changed;
        }

        bool exchangeDataNeighbours(labelList& globalGroupLabel)
        {
            bool changed(false);

            typedef std::map<label, LongList<lPairLabel> > dMap;
            dMap exchangeData;
            typedef std::map<lPair, std::set<lPair> > gMap;
            typedef std::set<lPair> pairSet;
            forAllConstIter(gMap, groupMap_, it)
            {
                if( it->first.first != Pstream::myProcNo() )
                    continue;

                const label groupLabel = globalGroupLabel[it->first.second];
                const pairSet& neiGroups = it->second;
                forAllConstIter(pairSet, neiGroups, sIt)
                {
                    const label neiProcNo = sIt->first;

                    if( neiProcNo == Pstream::myProcNo() )
                        continue;

                    dMap::iterator mapIt = exchangeData.find(neiProcNo);
                    if( mapIt == exchangeData.end() )
                    {
                        exchangeData.insert
                        (
                            std::make_pair(neiProcNo, LongList<lPairLabel>())
                        );
                        mapIt = exchangeData.find(neiProcNo);
                    }

                    mapIt->second.append(lPairLabel(*sIt, groupLabel));
                }
            }

            LongList<lPairLabel> receivedGroups;
            help::exchangeMap(exchangeData, receivedGroups);

            changed = updateGlobalGroupLabels(globalGroupLabel, receivedGroups);
            reduce(changed, maxOp<bool>());

            return changed;
        }
};

void checkCellConnectionsOverFaces::findCellGroups()
{
    Info << "Checking cell connections" << endl;

    //mesh_.write();
    //returnReduce(1, sumOp<label>());
    //::exit(1);

    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    labelListPMG front, communicationFaces;
    label chunkI(0), nChunks, chunkSize;
    VRWGraph neighbouringGroups;

    # ifdef USE_OMP
    # pragma omp parallel if( cells.size() > 1000 ) \
    private(front, communicationFaces)
    # endif
    {
        //- set the number of chunks and their size
        //- the number of chunks is greater than the number of threads
        //- in order to enhance load balancing
        # ifdef USE_OMP
        # pragma omp master
        {
            nChunks = 3 * omp_get_num_threads();
            chunkSize = cells.size() / nChunks + 1;
        }

        # pragma omp flush(nChunks, chunkSize)

        # pragma omp barrier
        # else
        nChunks = 1;
        chunkSize = cells.size();
        # endif

        while( chunkI < nChunks )
        {
            label minCell, maxCell;
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            minCell = chunkI++ * chunkSize;

            if( minCell >= cells.size() )
                break;

            maxCell = Foam::min(cells.size(), minCell + chunkSize);

            for(label cellI=minCell;cellI<maxCell;++cellI)
            {
                if( cellGroup_[cellI] != -1 )
                    continue;

                label groupI;
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                groupI = nGroups_++;

                front.clear();
                front.append(cellI);
                cellGroup_[cellI] = groupI;

                while( front.size() )
                {
                    const label fLabel = front.removeLastElement();

                    const cell& c = cells[fLabel];

                    forAll(c, fI)
                    {
                        label nei = owner[c[fI]];
                        if( nei == fLabel )
                            nei = neighbour[c[fI]];

                        if( nei < 0 )
                            continue;

                        if( (nei < minCell) || (nei >= maxCell) )
                        {
                            communicationFaces.append(c[fI]);
                        }
                        else if( cellGroup_[nei] == -1 )
                        {
                            cellGroup_[nei] = groupI;
                            front.append(nei);
                        }
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        {
            neighbouringGroups.setSize(nGroups_);
            forAll(neighbouringGroups, groupI)
                neighbouringGroups.append(groupI, groupI);
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- find group to neighbouring groups addressing
        forAll(communicationFaces, cfI)
        {
            const label faceI = communicationFaces[cfI];
            const label groupI = cellGroup_[owner[faceI]];
            const label neiGroup = cellGroup_[neighbour[faceI]];

            if( (neiGroup >= nGroups_) || (groupI >= nGroups_) )
                FatalError << "neiGroup " << neiGroup
                    << " groupI " << groupI << " are >= than "
                    << "nGroups " << nGroups_ << abort(FatalError);

            if( neiGroup != -1 )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                {
                    if( !neighbouringGroups.contains(groupI, neiGroup) )
                        neighbouringGroups.append(groupI, neiGroup);
                }
            }
        }
    }

    //- check global label for each group
    chunkI = 0;
    if( Pstream::parRun() )
    {
        labelList nGroupsAtProc(Pstream::nProcs());
        nGroupsAtProc[Pstream::myProcNo()] = nGroups_;

        Pstream::gatherList(nGroupsAtProc);
        Pstream::scatterList(nGroupsAtProc);

        for(label procI=0;procI<Pstream::myProcNo();++procI)
            chunkI += nGroupsAtProc[procI];
    }

    globalGroupLabel_.setSize(nGroups_);
    forAll(globalGroupLabel_, i)
        globalGroupLabel_[i] = chunkI++;

    //- check global group labels
    //- resolve groups for SMP parallelisation
    bool changed;
    do
    {
        changed = false;

        forAllReverse(neighbouringGroups, groupI)
        {
            label minGroup = globalGroupLabel_[groupI];

            forAllRow(neighbouringGroups, groupI, i)
            {
                const label currGroup = neighbouringGroups(groupI, i);
                minGroup = Foam::min(minGroup, globalGroupLabel_[currGroup]);
            }

            forAllRow(neighbouringGroups, groupI, i)
            {
                const label currGroup = neighbouringGroups(groupI, i);

                if( globalGroupLabel_[currGroup] != minGroup )
                {
                    globalGroupLabel_[currGroup] = minGroup;
                    changed = true;
                }
            }
        }
    } while( changed );

    if( Pstream::parRun() )
    {
        //- send global group labels
        const PtrList<writeProcessorPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- create ranges on processor boundaries
        //- each pair represent the groups which are neighbours over
        //- processor boundaries
        std::map<lPair, std::set<lPair> > groupMap;

        //- insert local groups
        forAll(neighbouringGroups, groupI)
        {
            std::set<lPair>& neiGroups =
                groupMap[std::make_pair(Pstream::myProcNo(), groupI)];
            forAllRow(neighbouringGroups, groupI, i)
            {
                neiGroups.insert
                (
                    std::make_pair
                    (
                        Pstream::myProcNo(),
                        neighbouringGroups(groupI, i)
                    )
                );
            }
        }

        //- create ranges for sending
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            LongList<labelledPair> dataToSend;

            label currGroup = cellGroup_[owner[start]];
            label groupStart(0);
            for(label faceI=start;faceI<end;++faceI)
            {
                if( currGroup != cellGroup_[owner[faceI]] )
                {
                    labelPair lp(groupStart, faceI-start);
                    dataToSend.append(labelledPair(currGroup, lp));

                    currGroup = cellGroup_[owner[faceI]];
                    groupStart = faceI - start;
                }
            }

            dataToSend.append
            (
                labelledPair(currGroup, labelPair(groupStart, end-start))
            );

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
        }

        //- receive ranges and create neighbour pairs
        forAll(procBoundaries, patchI)
        {
            List<labelledPair> receivedGroups;

            const writeProcessorPatch& procPatch = procBoundaries[patchI];

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procPatch.neiProcNo()
            );
            fromOtherProc >> receivedGroups;

            const label start = procPatch.patchStart();
            forAll(receivedGroups, i)
            {
                const label neiGroup = receivedGroups[i].pairLabel();
                const labelPair& pair = receivedGroups[i].pair();

                label currGroup = cellGroup_[owner[start+pair.first()]];
                for(label fI=pair.first();fI<pair.second();++fI)
                {
                    const label groupI = cellGroup_[owner[start+fI]];
                    std::set<lPair>& createdGroupPairs =
                        groupMap[lPair(Pstream::myProcNo(), currGroup)];

                    if( groupI != currGroup )
                    {
                        createdGroupPairs.insert
                        (
                            lPair(procPatch.neiProcNo(), neiGroup)
                        );
                        currGroup = groupI;
                    }
                }

                groupMap[lPair(Pstream::myProcNo(), currGroup)].insert
                (
                    lPair(procPatch.neiProcNo(), neiGroup)
                );
            }
        }

        //- resolve local groups based on inter-procesor connections
        //- if two group of a processor share the same group at some other
        //- processor then these two groups are the same group
        groupingOp gop(groupMap);
        gop.updateGlobalGroupLabels(globalGroupLabel_);

        //- the process of determining global group labels is iterative
        //- and is performed in multigrid-like manner
        //- a single iteration of the process consists of the following:
        //- 1. send data to the processors below in the tree structure
        //- 2. Each processor update global group labels based on neighbour data
        //- and the group labels received from other processors
        //- 3. send data to the processors above in the tree structure
        do
        {
            changed = false;

            if( gop.exchangeDataTopToBottom(globalGroupLabel_) )
                changed = true;

            if( gop.exchangeDataNeighbours(globalGroupLabel_) )
                changed = true;

            if( gop.exchangeDataBottomToTop(globalGroupLabel_) )
                changed = true;

            reduce(changed, maxOp<bool>());
        } while( changed );
    }

    nGroups_ = Foam::max(globalGroupLabel_) + 1;
    reduce(nGroups_, maxOp<label>());

    Info << "Finished checking cell connections" << endl;
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkCellConnectionsOverFaces::checkCellConnectionsOverFaces(polyMeshGen& mesh)
:
    mesh_(mesh),
    cellGroup_(mesh.cells().size(), -1),
    globalGroupLabel_(),
    nGroups_(0)
{
    findCellGroups();
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkCellConnectionsOverFaces::~checkCellConnectionsOverFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool checkCellConnectionsOverFaces::checkCellGroups()
{
    if( nGroups_ == 1 )
        return false;

    Warning << "Mesh has " << nGroups_ << " unconnected regions" << endl;

    labelList nCellsInGroup(nGroups_, 0);

    forAll(cellGroup_, cI)
        ++nCellsInGroup[globalGroupLabel_[cellGroup_[cI]]];

    if( Pstream::parRun() )
    {
        forAll(nCellsInGroup, groupI)
            reduce(nCellsInGroup[groupI], sumOp<label>());
    }

    //- find groups which has most cells this group will be kept
    label maxGroup(-1);
    forAll(nCellsInGroup, groupI)
        if( nCellsInGroup[groupI] > maxGroup )
        {
            maxGroup = nCellsInGroup[groupI];
            nGroups_ = groupI;
        }

    //- remove cells which are not in the group which has max num of cells
    boolList removeCell(mesh_.cells().size(), false);
    forAll(cellGroup_, cellI)
        if( globalGroupLabel_[cellGroup_[cellI]] != nGroups_ )
            removeCell[cellI] = true;

    polyMeshGenModifier(mesh_).removeCells(removeCell);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
