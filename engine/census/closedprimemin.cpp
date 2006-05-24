
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2006, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

#include <sstream>
#include "census/ncensus.h"
#include "census/ngluingpermsearcher.h"
#include "triangulation/nedge.h"
#include "triangulation/nfacepair.h"
#include "triangulation/ntriangulation.h"
#include "utilities/boostutils.h"
#include "utilities/memutils.h"

namespace regina {

const unsigned NClosedPrimeMinSearcher::EDGE_CHAIN_END = 1;
const unsigned NClosedPrimeMinSearcher::EDGE_CHAIN_INTERNAL_FIRST = 2;
const unsigned NClosedPrimeMinSearcher::EDGE_CHAIN_INTERNAL_SECOND = 3;
const unsigned NClosedPrimeMinSearcher::EDGE_DOUBLE_FIRST = 4;
const unsigned NClosedPrimeMinSearcher::EDGE_DOUBLE_SECOND = 5;
const unsigned NClosedPrimeMinSearcher::EDGE_MISC = 6;

const char NClosedPrimeMinSearcher::VLINK_CLOSED = 1;
const char NClosedPrimeMinSearcher::VLINK_NON_SPHERE = 2;

const char NClosedPrimeMinSearcher::ECLASS_TWISTED = 1;
const char NClosedPrimeMinSearcher::ECLASS_LOWDEG = 2;
const char NClosedPrimeMinSearcher::ECLASS_HIGHDEG = 4;
const char NClosedPrimeMinSearcher::ECLASS_CONE = 8;
const char NClosedPrimeMinSearcher::ECLASS_L31 = 16;

const int NClosedPrimeMinSearcher::vertexLinkNextFace[4][4] = {
    { -1, 2, 3, 1},
    { 3, -1, 0, 2},
    { 1, 3, -1, 0},
    { 1, 2, 0, -1}
};

const unsigned NClosedPrimeMinSearcher::coneEdge[12][2] = {
    { 0, 1 }, { 0, 2 }, { 1, 2 }, { 0, 3 }, { 0, 4 }, { 3, 4 },
    { 1, 3 }, { 1, 5 }, { 3, 5 }, { 2, 4 }, { 2, 5 }, { 4, 5 },
};

const char NClosedPrimeMinSearcher::coneNoTwist[12] = {
    1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1
};

const char NClosedPrimeMinSearcher::dataTag_ = 'c';

void NClosedPrimeMinSearcher::TetVertexState::dumpData(std::ostream& out)
        const {
    // Be careful with twistUp, which is a char but which should be
    // written as an int.
    out << parent << ' ' << rank << ' ' << bdry << ' '
        << (twistUp ? 1 : 0) << ' ' << (hadEqualRank ? 1 : 0) << ' '
        << static_cast<int>(bdryEdges) << ' '
        << bdryNext[0] << ' ' << bdryNext[1] << ' '
        << static_cast<int>(bdryTwist[0]) << ' '
        << static_cast<int>(bdryTwist[1]) << ' '
        << bdryNextOld[0] << ' ' << bdryNextOld[1] << ' '
        << static_cast<int>(bdryTwistOld[0]) << ' '
        << static_cast<int>(bdryTwistOld[1]);
}

bool NClosedPrimeMinSearcher::TetVertexState::readData(std::istream& in,
        unsigned long nStates) {
    in >> parent >> rank >> bdry;

    // twistUp is a char, but we need to read it as an int.
    int twist;
    in >> twist;
    twistUp = twist;

    // hadEqualRank is a bool, but we need to read it as an int.
    int bRank;
    in >> bRank;
    hadEqualRank = bRank;

    // More chars to ints coming.
    int bVal;

    in >> bVal; bdryEdges = bVal;
    in >> bdryNext[0] >> bdryNext[1];
    in >> bVal; bdryTwist[0] = bVal;
    in >> bVal; bdryTwist[1] = bVal;
    in >> bdryNextOld[0] >> bdryNextOld[1];
    in >> bVal; bdryTwistOld[0] = bVal;
    in >> bVal; bdryTwistOld[1] = bVal;

    if (parent < -1 || parent >= static_cast<long>(nStates))
        return false;
    if (rank >= nStates)
        return false;
    if (bdry > 3 * nStates)
        return false;
    if (twist != 1 && twist != 0)
        return false;
    if (bRank != 1 && bRank != 0)
        return false;
    if (bdryEdges > 3) /* Never < 0 since this is unsigned. */
        return false;
    if (bdryNext[0] < 0 || bdryNext[0] >= static_cast<long>(nStates))
        return false;
    if (bdryNext[1] < 0 || bdryNext[1] >= static_cast<long>(nStates))
        return false;
    if (bdryNextOld[0] < -1 || bdryNext[0] >= static_cast<long>(nStates))
        return false;
    if (bdryNextOld[1] < -1 || bdryNextOld[1] >= static_cast<long>(nStates))
        return false;
    if (bdryTwist[0] < 0 || bdryTwist[0] > 1)
        return false;
    if (bdryTwist[1] < 0 || bdryTwist[1] > 1)
        return false;
    if (bdryTwistOld[0] < 0 || bdryTwistOld[0] > 1)
        return false;
    if (bdryTwistOld[1] < 0 || bdryTwistOld[1] > 1)
        return false;

    return true;
}

void NClosedPrimeMinSearcher::TetEdgeState::dumpData(std::ostream& out) const {
    // Be careful with twistUp, which is a char but which should be
    // written as an int.
    out << parent << ' ' << rank << ' ' << size << ' '
        << (bounded ? 1 : 0) << ' ' << (twistUp ? 1 : 0) << ' '
        << (hadEqualRank ? 1 : 0);
}

bool NClosedPrimeMinSearcher::TetEdgeState::readData(std::istream& in,
        unsigned long nStates) {
    in >> parent >> rank >> size;

    // bounded is a bool, but we need to read it as an int.
    int bBounded;
    in >> bBounded;
    bounded = bBounded;

    // twistUp is a char, but we need to read it as an int.
    int twist;
    in >> twist;
    twistUp = twist;

    // hadEqualRank is a bool, but we need to read it as an int.
    int bRank;
    in >> bRank;
    hadEqualRank = bRank;

    if (parent < -1 || parent >= static_cast<long>(nStates))
        return false;
    if (rank >= nStates)
        return false;
    if (size >= nStates)
        return false;
    if (bBounded != 1 && bBounded != 0)
        return false;
    if (twist != 1 && twist != 0)
        return false;
    if (bRank != 1 && bRank != 0)
        return false;

    return true;
}

NClosedPrimeMinSearcher::NClosedPrimeMinSearcher(const NFacePairing* pairing,
        const NFacePairingIsoList* autos, bool orientableOnly,
        UseGluingPerms use, void* useArgs) :
        NGluingPermSearcher(pairing, autos, orientableOnly,
            true /* finiteOnly */,
            NCensus::PURGE_NON_MINIMAL_PRIME | NCensus::PURGE_P2_REDUCIBLE,
            use, useArgs) {
    initOrder();
}

void NClosedPrimeMinSearcher::initOrder() {
    // Preconditions:
    //     Only closed prime minimal P2-irreducible triangulations are needed.
    //     The given face pairing is closed with order >= 3.

    // ---------- Selecting an ordering of faces ----------

    // We fill permutations in the order:
    //     1. One-ended chains (== layered solid tori) from loop to
    //        boundary, though chains may be interlaced in the
    //        processing order;
    //     2. Everything else ordered by tetrahedron faces.
    //
    // Both permutations for each double edge will be processed
    // consecutively, the permutation for the smallest face involved
    // in the double edge being processed first.
    //
    // Note from the tests above that there are no triple edges.

    unsigned nTets = getNumberOfTetrahedra();

    order = new NTetFace[nTets * 2];
    orderType = new unsigned[nTets * 2];

    bool* orderAssigned = new bool[nTets * 4];
        /**< Have we placed a tetrahedron face or its partner in the
             order[] array yet? */

    // Hunt for structures within the face pairing graph.

    NTetFace face, adj;
    unsigned orderDone = 0;
    std::fill(orderAssigned, orderAssigned + 4 * nTets, false);

    // Begin by searching for tetrahedra that are joined to themselves.
    // Note that each tetrahedra can be joined to itself at most once,
    // since we are guaranteed that the face pairing is connected with
    // order >= 3.

    for (face.setFirst(); ! face.isPastEnd(nTets, true); face++) {
        if (orderAssigned[face.tet * 4 + face.face])
            continue;

        adj = (*pairing)[face];
        if (adj.tet != face.tet)
            continue;

        order[orderDone] = face;
        orderType[orderDone] = EDGE_CHAIN_END;
        orderAssigned[face.tet * 4 + face.face] = true;
        orderAssigned[adj.tet * 4 + adj.face] = true;
        orderDone++;
    }

    // Record the number of one-ended chains.

    unsigned nChains = orderDone;

    // Continue by following each one-ended chain whose base was
    // identified in the previous loop.

    unsigned i;
    int tet;
    NTetFace dest1, dest2;
    NFacePair faces;
    for (i = 0; i < nChains; i++) {
        tet = order[i].tet;
        faces = NFacePair(order[i].face,
            (*pairing)[order[i]].face).complement();
        dest1 = pairing->dest(tet, faces.lower());
        dest2 = pairing->dest(tet, faces.upper());

        // Currently tet and faces refer to the two faces of the base
        // tetrahedron that are pointing outwards.
        while (dest1.tet == dest2.tet && dest1.tet != tet &&
                (! orderAssigned[tet * 4 + faces.lower()]) &&
                (! orderAssigned[tet * 4 + faces.upper()])) {
            // Insert this pair of edges into the ordering and follow
            // the chain.
            orderType[orderDone] = EDGE_CHAIN_INTERNAL_FIRST;
            orderType[orderDone + 1] = EDGE_CHAIN_INTERNAL_SECOND;

            if (tet < dest1.tet) {
                order[orderDone] = NTetFace(tet, faces.lower());
                order[orderDone + 1] = NTetFace(tet, faces.upper());
            }

            orderAssigned[tet * 4 + faces.lower()] = true;
            orderAssigned[tet * 4 + faces.upper()] = true;
            orderAssigned[dest1.tet * 4 + dest1.face] = true;
            orderAssigned[dest2.tet * 4 + dest2.face] = true;

            faces = NFacePair(dest1.face, dest2.face);

            if (dest1.tet < tet) {
                order[orderDone] = NTetFace(dest1.tet, faces.lower());
                order[orderDone + 1] = NTetFace(dest1.tet, faces.upper());
            }

            faces = faces.complement();
            tet = dest1.tet;

            dest1 = pairing->dest(tet, faces.lower());
            dest2 = pairing->dest(tet, faces.upper());

            orderDone += 2;
        }
    }

    // Record the number of edges in the face pairing graph
    // belonging to one-ended chains.
    nChainEdges = orderDone;

    // Run through the remaining faces.
    for (face.setFirst(); ! face.isPastEnd(nTets, true); face++)
        if (! orderAssigned[face.tet * 4 + face.face]) {
            order[orderDone] = face;
            if (face.face < 3 && pairing->dest(boost::next(face)).tet ==
                    pairing->dest(face).tet)
                orderType[orderDone] = EDGE_DOUBLE_FIRST;
            else if (face.face > 0 && pairing->dest(boost::prior(face)).tet ==
                    pairing->dest(face).tet)
                orderType[orderDone] = EDGE_DOUBLE_SECOND;
            else
                orderType[orderDone] = EDGE_MISC;
            orderDone++;

            adj = (*pairing)[face];
            orderAssigned[face.tet * 4 + face.face] = true;
            orderAssigned[adj.tet * 4 + adj.face] = true;
        }

    // All done for the order[] array.  Tidy up.
    delete[] orderAssigned;

    // ---------- Calculating the possible gluing permutations ----------

    // For each face in the order[] array of type EDGE_CHAIN_END or
    // EDGE_CHAIN_INTERNAL_FIRST, we calculate the two gluing permutations
    // that must be tried.
    //
    // For the remaining faces we try all possible permutations.

    chainPermIndices = (nChainEdges == 0 ? 0 : new int[nChainEdges * 2]);

    NFacePair facesAdj, comp, compAdj;
    NPerm trial1, trial2;
    for (i = 0; i < nChainEdges; i++) {
        if (orderType[i] == EDGE_CHAIN_END) {
            faces = NFacePair(order[i].face, pairing->dest(order[i]).face);
            comp = faces.complement();

            // order[i].face == faces.lower(),
            // pairing->dest(order[i]).face == faces.upper().
            chainPermIndices[2 * i] = gluingToIndex(order[i],
                NPerm(faces.lower(), faces.upper(),
                      faces.upper(), comp.lower(),
                      comp.lower(), comp.upper(),
                      comp.upper(), faces.lower()));
            chainPermIndices[2 * i + 1] = gluingToIndex(order[i],
                NPerm(faces.lower(), faces.upper(),
                      faces.upper(), comp.upper(),
                      comp.upper(), comp.lower(),
                      comp.lower(), faces.lower()));
        } else if (orderType[i] == EDGE_CHAIN_INTERNAL_FIRST) {
            faces = NFacePair(order[i].face, order[i + 1].face);
            comp = faces.complement();
            facesAdj = NFacePair(pairing->dest(order[i]).face,
                pairing->dest(order[i + 1]).face);
            compAdj = facesAdj.complement();

            // order[i].face == faces.lower(),
            // order[i + 1].face == faces.upper(),
            // pairing->dest(order[i]).face == facesAdj.lower().
            // pairing->dest(order[i + 1]).face == facesAdj.upper().
            trial1 = NPerm(faces.lower(), facesAdj.lower(),
                           faces.upper(), compAdj.lower(),
                           comp.lower(), compAdj.upper(),
                           comp.upper(), facesAdj.upper());
            trial2 = NPerm(faces.lower(), facesAdj.lower(),
                           faces.upper(), compAdj.upper(),
                           comp.lower(), compAdj.lower(),
                           comp.upper(), facesAdj.upper());
            if (trial1.compareWith(trial2) < 0) {
                chainPermIndices[2 * i] = gluingToIndex(order[i], trial1);
                chainPermIndices[2 * i + 2] = gluingToIndex(order[i + 1],
                    NPerm(faces.lower(), compAdj.upper(),
                          faces.upper(), facesAdj.upper(),
                          comp.lower(), facesAdj.lower(),
                          comp.upper(), compAdj.lower()));
            } else {
                chainPermIndices[2 * i] = gluingToIndex(order[i], trial2);
                chainPermIndices[2 * i + 2] = gluingToIndex(order[i + 1],
                    NPerm(faces.lower(), compAdj.lower(),
                          faces.upper(), facesAdj.upper(),
                          comp.lower(), facesAdj.lower(),
                          comp.upper(), compAdj.upper()));
            }

            trial1 = NPerm(faces.lower(), facesAdj.lower(),
                           faces.upper(), compAdj.lower(),
                           comp.lower(), facesAdj.upper(),
                           comp.upper(), compAdj.upper());
            trial2 = NPerm(faces.lower(), facesAdj.lower(),
                           faces.upper(), compAdj.upper(),
                           comp.lower(), facesAdj.upper(),
                           comp.upper(), compAdj.lower());
            if (trial1.compareWith(trial2) < 0) {
                chainPermIndices[2 * i + 1] = gluingToIndex(order[i], trial1);
                chainPermIndices[2 * i + 3] = gluingToIndex(order[i + 1],
                    NPerm(faces.lower(), compAdj.upper(),
                          faces.upper(), facesAdj.upper(),
                          comp.lower(), compAdj.lower(),
                          comp.upper(), facesAdj.lower()));
            } else {
                chainPermIndices[2 * i + 1] = gluingToIndex(order[i], trial2);
                chainPermIndices[2 * i + 3] = gluingToIndex(order[i + 1],
                    NPerm(faces.lower(), compAdj.lower(),
                          faces.upper(), facesAdj.upper(),
                          comp.lower(), compAdj.upper(),
                          comp.upper(), facesAdj.lower()));
            }
        }
    }

    // ---------- Tracking of vertex / edge equivalence classes ----------

    nVertexClasses = nTets * 4;
    vertexState = new TetVertexState[nTets * 4];
    vertexStateChanged = new int[nTets * 8];
    std::fill(vertexStateChanged, vertexStateChanged + nTets * 8, -1);
    for (i = 0; i < nTets * 4; i++) {
        vertexState[i].bdryEdges = 3;
        vertexState[i].bdryNext[0] = vertexState[i].bdryNext[1] = i;
        vertexState[i].bdryTwist[0] = vertexState[i].bdryTwist[1] = 0;
        // Initialise the backup members also so we're not writing
        // uninitialised data via dumpData().
        vertexState[i].bdryNextOld[0] = vertexState[i].bdryNextOld[1] = -1;
        vertexState[i].bdryTwistOld[0] = vertexState[i].bdryTwistOld[1] = 0;
    }

    nEdgeClasses = nTets * 6;
    edgeState = new TetEdgeState[nTets * 6];
    edgeStateChanged = new int[nTets * 8];
    std::fill(edgeStateChanged, edgeStateChanged + nTets * 8, -1);

#if PRUNE_HIGH_DEG_EDGE_SET
    highDegSum = 0;
    highDegBound = 3 * nTets - 3;
#endif
}

// TODO (net): See what was removed when we brought in vertex link checking.
void NClosedPrimeMinSearcher::runSearch(long maxDepth) {
    // Preconditions:
    //     Only closed prime minimal P2-irreducible triangulations are needed.
    //     The given face pairing is closed with order >= 3.

    unsigned nTets = getNumberOfTetrahedra();
    if (maxDepth < 0) {
        // Larger than we will ever see (and in fact grossly so).
        maxDepth = nTets * 4 + 1;
    }

    if (! started) {
        // Search initialisation.
        started = true;

        // Begin by testing for face pairings that can never lead to such a
        // triangulation.
        if (pairing->hasTripleEdge() ||
                pairing->hasBrokenDoubleEndedChain() ||
                pairing->hasOneEndedChainWithDoubleHandle() ||
                pairing->hasOneEndedChainWithStrayBigon() ||
                pairing->hasWedgedDoubleEndedChain() ||
                pairing->hasTripleOneEndedChain()) {
            use_(0, useArgs_);
            return;
        }

        orderElt = 0;
        if (nChainEdges < nTets * 2)
            orientation[order[nChainEdges].tet] = 1;
    }

    // Is it a partial search that has already finished?
    if (orderElt == static_cast<int>(nTets) * 2) {
        if (isCanonical())
            use_(this, useArgs_);
        use_(0, useArgs_);
        return;
    }

    // ---------- Selecting the individual gluing permutations ----------

    // Observe that in a canonical face pairing, one-ended chains always
    // follow an increasing sequence of tetrahedra from boundary to end,
    // or follow the sequence of tetrahedra 0, 1, ..., k from end to
    // boundary.
    //
    // In particular, this means that for any tetrahedron not internal
    // to a one-ended chain (with the possible exception of tetrahedron
    // order[nChainEdges].tet), face 0 of this tetrahedron is not
    // involved in a one-ended chain.

    // In this generation algorithm, each orientation is simply +/-1.
    // We won't bother assigning orientations to the tetrahedra internal
    // to the one-ended chains.

    int minOrder = orderElt;
    int maxOrder = orderElt + maxDepth;

    NTetFace face, adj;
    bool generic;
    int mergeResult;
    while (orderElt >= minOrder) {
        face = order[orderElt];
        adj = (*pairing)[face];

        // TODO (long-term): Check for cancellation.

        // Move to the next permutation.
        if (orderType[orderElt] == EDGE_CHAIN_END ||
                orderType[orderElt] == EDGE_CHAIN_INTERNAL_FIRST) {
            // Choose from one of the two permutations stored in array
            // chainPermIndices[].
            generic = false;
            if (permIndex(face) < 0)
                permIndex(face) = chainPermIndices[2 * orderElt];
            else if (permIndex(face) == chainPermIndices[2 * orderElt])
                permIndex(face) = chainPermIndices[2 * orderElt + 1];
            else
                permIndex(face) = 6;
        } else if (orderType[orderElt] == EDGE_CHAIN_INTERNAL_SECOND) {
            // The permutation is predetermined.
            generic = false;
            if (permIndex(face) < 0) {
                if (permIndex(order[orderElt - 1]) ==
                        chainPermIndices[2 * orderElt - 2])
                    permIndex(face) = chainPermIndices[2 * orderElt];
                else
                    permIndex(face) = chainPermIndices[2 * orderElt + 1];
            } else
                permIndex(face) = 6;
        } else {
            // Generic case.
            generic = true;

            // Be sure to preserve the orientation of the permutation if
            // necessary.
            if ((! orientableOnly_) || pairing->dest(face).face == 0)
                permIndex(face)++;
            else
                permIndex(face) += 2;
        }

        // Are we out of ideas for this face?
        if (permIndex(face) >= 6) {
            // Head back down to the previous face.
            permIndex(face) = -1;
            permIndex(adj) = -1;
            orderElt--;

            // Pull apart vertex and edge links at the previous level.
            if (orderElt >= minOrder) {
                splitVertexClasses();
                splitEdgeClasses();
            }

            continue;
        }

        // We are sitting on a new permutation to try.
        permIndex(adj) = allPermsS3Inv[permIndex(face)];

        // In the following code we use several results from
        // "Face pairing graphs and 3-manifold enumeration", B. A. Burton,
        // J. Knot Theory Ramifications 13 (2004).
        //
        // These include:
        //
        // - We cannot have an edge of degree <= 2, or an edge of degree 3
        //   meeting three distinct tetrahedra (section 2.1);
        // - We must have exactly one vertex (lemma 2.6);
        // - We cannot have a face with two edges identified to form a
        //   cone (lemma 2.8);
        // - We cannot have a face with all three edges identified to
        //   form an L(3,1) spine (lemma 2.5).

        // Merge edge links and run corresponding tests.
        if (mergeEdgeClasses()) {
            // We created a structure that should not appear in a final
            // census triangulation (e.g., a low-degree or invalid edge,
            // or a face whose edges are identified in certain ways).
            splitEdgeClasses();
            continue;
        }
        // The final triangulation should have precisely (nTets + 1) edges
        // (since it must have precisely one vertex).
        if (nEdgeClasses < nTets + 1) {
            // We already have too few edge classes, and the count can
            // only get smaller.
            // Note that the triangulations we are pruning include ideal
            // triangulations (with vertex links of Euler characteristic < 2).
            splitEdgeClasses();
            continue;
        }
        // In general, one can prove that (assuming no invalid edges or
        // boundary faces) we will end up with (<= nTets + nVertices) edges
        // (with strictly fewer edges if some vertex links are non-spherical).
        // If we must end up with (> nTets + 1) edges we can therefore
        // prune since we won't have a one-vertex triangulation.
        if (nEdgeClasses > nTets + 1 + 3 * (nTets * 2 - orderElt - 1)) {
            // We have (2n - orderElt - 1) more gluings to choose.
            // Since each merge can reduce the number of edge classes
            // by at most 3, there is no way we can end up with just
            // (nTets + 1) edges at the end.
            splitEdgeClasses();
            continue;
        }

        // Merge vertex links and run corresponding tests.
        mergeResult = mergeVertexClasses();
        if (mergeResult & VLINK_CLOSED) {
            // We closed off a vertex link, which means we will end up
            // with more than one vertex (unless this was our very last
            // gluing).
            if (orderElt + 1 < static_cast<int>(nTets) * 2) {
                splitVertexClasses();
                splitEdgeClasses();
                continue;
            }
        }
        if (mergeResult & VLINK_NON_SPHERE) {
            // Our vertex link will never be a 2-sphere.  Stop now.
            splitVertexClasses();
            splitEdgeClasses();
            continue;
        }
        if (nVertexClasses > 1 + 3 * (nTets * 2 - orderElt - 1)) {
            // We have (2n - orderElt - 1) more gluings to choose.
            // Since each merge can reduce the number of vertex classes
            // by at most 3, there is no way we can end up with just one
            // vertex at the end.
            splitVertexClasses();
            splitEdgeClasses();
            continue;
        }

        // Fix the orientation if appropriate.
        if (generic && adj.face == 0 && orientableOnly_) {
            // It's the first time we've hit this tetrahedron.
            if ((permIndex(face) + (face.face == 3 ? 0 : 1) +
                    (adj.face == 3 ? 0 : 1)) % 2 == 0)
                orientation[adj.tet] = -orientation[face.tet];
            else
                orientation[adj.tet] = orientation[face.tet];
        }

        // Move on to the next face.
        orderElt++;

        // If we're at the end, try the solution and step back.
        if (orderElt == static_cast<int>(nTets) * 2) {
            // We in fact have an entire triangulation.
            // Run through the automorphisms and check whether our
            // permutations are in canonical form.
            if (isCanonical())
                use_(this, useArgs_);

            // Back to the previous face.
            orderElt--;

            // Pull apart vertex and edge links at the previous level.
            if (orderElt >= minOrder) {
                splitVertexClasses();
                splitEdgeClasses();
            }
        } else {
            // Not a full triangulation; just one level deeper.

            // We've moved onto a new face.
            // Be sure to get the orientation right.
            face = order[orderElt];
            if (orientableOnly_ && pairing->dest(face).face > 0) {
                // permIndex(face) will be set to -1 or -2 as appropriate.
                adj = (*pairing)[face];
                if (orientation[face.tet] == orientation[adj.tet])
                    permIndex(face) = 1;
                else
                    permIndex(face) = 0;

                if ((face.face == 3 ? 0 : 1) + (adj.face == 3 ? 0 : 1) == 1)
                    permIndex(face) = (permIndex(face) + 1) % 2;

                permIndex(face) -= 2;
            }

            if (orderElt == maxOrder) {
                // We haven't found an entire triangulation, but we've
                // gone as far as we need to.
                // Process it, then step back.
                use_(this, useArgs_);

                // Back to the previous face.
                permIndex(face) = -1;
                orderElt--;

                // Pull apart vertex links at the previous level.
                if (orderElt >= minOrder) {
                    splitVertexClasses();
                    splitEdgeClasses();
                }
            }
        }
    }

    // And the search is over.

    // Some extra sanity checking.
    if (minOrder == 0) {
        // Our vertex classes had better be 4n standalone vertices.
        if (nVertexClasses != 4 * nTets)
            std::cerr << "ERROR: nVertexClasses == "
                << nVertexClasses << " at end of search!" << std::endl;
        for (int i = 0; i < static_cast<int>(nTets) * 4; i++) {
            if (vertexState[i].parent != -1)
                std::cerr << "ERROR: vertexState[" << i << "].parent == "
                    << vertexState[i].parent << " at end of search!"
                    << std::endl;
            if (vertexState[i].rank != 0)
                std::cerr << "ERROR: vertexState[" << i << "].rank == "
                    << vertexState[i].rank << " at end of search!" << std::endl;
            if (vertexState[i].bdry != 3)
                std::cerr << "ERROR: vertexState[" << i << "].bdry == "
                    << vertexState[i].bdry << " at end of search!" << std::endl;
            if (vertexState[i].hadEqualRank)
                std::cerr << "ERROR: vertexState[" << i << "].hadEqualRank == "
                    "true at end of search!" << std::endl;
            if (vertexState[i].bdryEdges != 3)
                std::cerr << "ERROR: vertexState[" << i << "].bdryEdges == "
                    << static_cast<int>(vertexState[i].bdryEdges)
                    << " at end of search!" << std::endl;
            if (vertexState[i].bdryNext[0] != i)
                std::cerr << "ERROR: vertexState[" << i << "].bdryNext[0] == "
                    << vertexState[i].bdryNext[0] << " at end of search!"
                    << std::endl;
            if (vertexState[i].bdryNext[1] != i)
                std::cerr << "ERROR: vertexState[" << i << "].bdryNext[1] == "
                    << vertexState[i].bdryNext[1] << " at end of search!"
                    << std::endl;
            if (vertexState[i].bdryTwist[0])
                std::cerr << "ERROR: vertexState[" << i << "].bdryTwist == "
                    << static_cast<int>(vertexState[i].bdryTwist[0])
                    << " at end of search!" << std::endl;
            if (vertexState[i].bdryTwist[1])
                std::cerr << "ERROR: vertexState[" << i << "].bdryTwist == "
                    << static_cast<int>(vertexState[i].bdryTwist[1])
                    << " at end of search!" << std::endl;
        }
        for (unsigned i = 0; i < nTets * 8; i++)
            if (vertexStateChanged[i] != -1)
                std::cerr << "ERROR: vertexStateChanged[" << i << "] == "
                    << vertexStateChanged[i] << " at end of search!"
                    << std::endl;

        // And our edge classes had better be 6n standalone edges.
        if (nEdgeClasses != 6 * nTets)
            std::cerr << "ERROR: nEdgeClasses == "
                << nEdgeClasses << " at end of search!" << std::endl;
        for (unsigned i = 0; i < nTets * 6; i++) {
            if (edgeState[i].parent != -1)
                std::cerr << "ERROR: edgeState[" << i << "].parent == "
                    << edgeState[i].parent << " at end of search!"
                    << std::endl;
            if (edgeState[i].rank != 0)
                std::cerr << "ERROR: edgeState[" << i << "].rank == "
                    << edgeState[i].rank << " at end of search!" << std::endl;
            if (edgeState[i].size != 1)
                std::cerr << "ERROR: edgeState[" << i << "].size == "
                    << edgeState[i].size << " at end of search!" << std::endl;
            if (! edgeState[i].bounded)
                std::cerr << "ERROR: edgeState[" << i << "].bounded == "
                    << edgeState[i].bounded << " at end of search!"
                    << std::endl;
            if (edgeState[i].hadEqualRank)
                std::cerr << "ERROR: edgeState[" << i << "].hadEqualRank == "
                    "true at end of search!" << std::endl;
        }
        for (unsigned i = 0; i < nTets * 8; i++)
            if (edgeStateChanged[i] != -1)
                std::cerr << "ERROR: edgeStateChanged[" << i << "] == "
                    << edgeStateChanged[i] << " at end of search!"
                    << std::endl;

#if PRUNE_HIGH_DEG_EDGE_SET
        if (highDegSum != 0)
            std::cerr << "ERROR: highDegSum == " << highDegSum
                << " at end of search!" << std::endl;
#endif
    }

    use_(0, useArgs_);
}

void NClosedPrimeMinSearcher::dumpData(std::ostream& out) const {
    NGluingPermSearcher::dumpData(out);

    unsigned nTets = getNumberOfTetrahedra();
    unsigned i;

    for (i = 0; i < 2 * nTets; i++) {
        if (i)
            out << ' ';
        out << order[i].tet << ' ' << order[i].face << ' ' << orderType[i];
    }
    out << std::endl;

    out << nChainEdges << std::endl;
    if (nChainEdges) {
        for (i = 0; i < 2 * nChainEdges; i++) {
            if (i)
                out << ' ';
            out << chainPermIndices[i];
        }
        out << std::endl;
    }

    out << orderElt << std::endl;

    out << nVertexClasses << std::endl;
    for (i = 0; i < 4 * nTets; i++) {
        vertexState[i].dumpData(out);
        out << std::endl;
    }
    for (i = 0; i < 8 * nTets; i++) {
        if (i)
            out << ' ';
        out << vertexStateChanged[i];
    }
    out << std::endl;

    out << nEdgeClasses << std::endl;
    for (i = 0; i < 6 * nTets; i++) {
        edgeState[i].dumpData(out);
        out << std::endl;
    }
    for (i = 0; i < 8 * nTets; i++) {
        if (i)
            out << ' ';
        out << edgeStateChanged[i];
    }
    out << std::endl;

#if PRUNE_HIGH_DEG_EDGE_SET
    out << highDegSum << ' ' << highDegBound << std::endl;
#endif
}

NClosedPrimeMinSearcher::NClosedPrimeMinSearcher(std::istream& in,
        UseGluingPerms use, void* useArgs) :
        NGluingPermSearcher(in, use, useArgs),
        order(0), orderType(0), nChainEdges(0), chainPermIndices(0),
        nVertexClasses(0), vertexState(0), vertexStateChanged(0),
        nEdgeClasses(0), edgeState(0), edgeStateChanged(0),
        orderElt(0) {
    if (inputError_)
        return;

    unsigned nTets = getNumberOfTetrahedra();
    unsigned i;

    order = new NTetFace[2 * nTets];
    orderType = new unsigned[nTets * 2];
    for (i = 0; i < 2 * nTets; i++) {
        in >> order[i].tet >> order[i].face >> orderType[i];
        if (order[i].tet >= static_cast<int>(nTets) || order[i].tet < 0 ||
                order[i].face >= 4 || order[i].face < 0) {
            inputError_ = true; return;
        }
    }

    in >> nChainEdges;
    /* Unnecessary since nChainEdges is unsigned.
    if (nChainEdges < 0) {
        inputError_ = true; return;
    } */
    if (nChainEdges) {
        chainPermIndices = new int[nChainEdges * 2];
        for (i = 0; i < 2 * nChainEdges; i++) {
            in >> chainPermIndices[i];
            if (chainPermIndices[i] < 0 || chainPermIndices[i] >= 6) {
                inputError_ = true; return;
            }
        }
    }

    in >> orderElt;

    in >> nVertexClasses;
    if (nVertexClasses > 4 * nTets) {
        inputError_ = true; return;
    }

    vertexState = new TetVertexState[4 * nTets];
    for (i = 0; i < 4 * nTets; i++)
        if (! vertexState[i].readData(in, 4 * nTets)) {
            inputError_ = true; return;
        }

    vertexStateChanged = new int[8 * nTets];
    for (i = 0; i < 8 * nTets; i++) {
        in >> vertexStateChanged[i];
        if (vertexStateChanged[i] < -1 ||
                 vertexStateChanged[i] >= 4 * static_cast<int>(nTets)) {
            inputError_ = true; return;
        }
    }

    in >> nEdgeClasses;
    if (nEdgeClasses > 6 * nTets) {
        inputError_ = true; return;
    }

    edgeState = new TetEdgeState[6 * nTets];
    for (i = 0; i < 6 * nTets; i++)
        if (! edgeState[i].readData(in, 6 * nTets)) {
            inputError_ = true; return;
        }

    edgeStateChanged = new int[8 * nTets];
    for (i = 0; i < 8 * nTets; i++) {
        in >> edgeStateChanged[i];
        if (edgeStateChanged[i] < -1 ||
                 edgeStateChanged[i] >= 6 * static_cast<int>(nTets)) {
            inputError_ = true; return;
        }
    }

#if PRUNE_HIGH_DEG_EDGE_SET
    in >> highDegSum >> highDegBound;
    if (highDegSum < 0 || highDegSum > 6 * static_cast<int>(nTets) ||
            highDegBound != 3 * static_cast<int>(nTets) - 3) {
        inputError_ = true; return;
    }
#endif

    // Did we hit an unexpected EOF?
    if (in.eof())
        inputError_ = true;
}

int NClosedPrimeMinSearcher::mergeVertexClasses() {
    // Merge all three vertex pairs for the current face.
    NTetFace face = order[orderElt];
    NTetFace adj = (*pairing)[face];

    int retVal = 0;

    int v, w;
    int vIdx, wIdx, tmpIdx, nextIdx;
    unsigned orderIdx;
    int vRep, wRep;
    int vNext[2], wNext[2];
    char vTwist[2], wTwist[2];
    NPerm p = gluingPerm(face);
    char parentTwists, hasTwist, tmpTwist;
    for (v = 0; v < 4; v++) {
        if (v == face.face)
            continue;

        w = p[v];
        vIdx = v + 4 * face.tet;
        wIdx = w + 4 * adj.tet;
        orderIdx = v + 4 * orderElt;

        // Are the natural 012 representations of the two faces joined
        // with reversed orientations?
        // Here we combine the sign of permutation p with the mappings
        // from 012 to the native tetrahedron vertices, i.e., v <-> 3 and
        // w <-> 3.
        hasTwist = (p.sign() < 0 ? 0 : 1);
        if ((v == 3 && w != 3) || (v != 3 && w == 3))
            hasTwist ^= 1;

        parentTwists = 0;
        for (vRep = vIdx; vertexState[vRep].parent >= 0;
                vRep = vertexState[vRep].parent)
            parentTwists ^= vertexState[vRep].twistUp;
        for (wRep = wIdx; vertexState[wRep].parent >= 0;
                wRep = vertexState[wRep].parent)
            parentTwists ^= vertexState[wRep].twistUp;

        if (vRep == wRep) {
            vertexState[vRep].bdry -= 2;
            if (vertexState[vRep].bdry == 0)
                retVal |= VLINK_CLOSED;

            // Have we made the vertex link non-orientable?
            if (hasTwist ^ parentTwists)
                retVal |= VLINK_NON_SPHERE;

            vertexStateChanged[orderIdx] = -1;

            // Examine the cycles of boundary components.
            if (vIdx == wIdx) {
                // Ignore this case; it implies either a one-face cone
                // or a low degree edge, both of which should have
                // already been picked up in the edge link tests.
                std::cerr << "ERROR: vIdx == wIdx" << std::endl;
            } else {
                // We are joining two distinct tetrahedron vertices that
                // already contribute to the same vertex link.
                if (vertexState[vIdx].bdryEdges == 2)
                    vtxBdryBackup(vIdx);
                if (vertexState[wIdx].bdryEdges == 2)
                    vtxBdryBackup(wIdx);

                if (vtxBdryLength1(vIdx) && vtxBdryLength1(wIdx)) {
                    // We are joining together two boundaries of length one.
                    // Do nothing and mark the non-trivial genus.
                    // std::cerr << "NON-SPHERE: 1 >-< 1" << std::endl;
                    retVal |= VLINK_NON_SPHERE;
                } else if (vtxBdryLength2(vIdx, wIdx)) {
                    // We are closing off a single boundary of length two.
                    // All good.
                } else {
                    vtxBdryNext(vIdx, face.tet, v, face.face, vNext, vTwist);
                    vtxBdryNext(wIdx, adj.tet, w, adj.face, wNext, wTwist);

                    if (vNext[0] == wIdx && wNext[1 ^ vTwist[0]] == vIdx) {
                        // We are joining two adjacent edges of the vertex link.
                        // Simply eliminate them.
                        vtxBdryJoin(vNext[1], 0 ^ vTwist[1],
                            wNext[0 ^ vTwist[0]],
                            (vTwist[0] ^ wTwist[0 ^ vTwist[0]]) ^ vTwist[1]);
                    } else if (vNext[1] == wIdx &&
                            wNext[0 ^ vTwist[1]] == vIdx) {
                        // Again, joining two adjacent edges of the vertex link.
                        vtxBdryJoin(vNext[0], 1 ^ vTwist[0],
                            wNext[1 ^ vTwist[1]],
                            (vTwist[1] ^ wTwist[1 ^ vTwist[1]]) ^ vTwist[0]);
                    } else {
                        // See if we are joining two different boundary cycles
                        // together; if so, we have created non-trivial genus in
                        // the vertex link.
                        tmpIdx = vertexState[vIdx].bdryNext[0];
                        tmpTwist = vertexState[vIdx].bdryTwist[0];
                        while (tmpIdx != vIdx && tmpIdx != wIdx) {
                            nextIdx = vertexState[tmpIdx].
                                bdryNext[0 ^ tmpTwist];
                            tmpTwist ^= vertexState[tmpIdx].
                                bdryTwist[0 ^ tmpTwist];
                            tmpIdx = nextIdx;
                        }

                        if (tmpIdx == vIdx) {
                            // Different boundary cycles.
                            // Don't touch anything; just flag a
                            // high genus error.
                            // std::cerr << "NON-SPHERE: (X)" << std::endl;
                            retVal |= VLINK_NON_SPHERE;
                        } else {
                            // Same boundary cycle.
                            vtxBdryJoin(vNext[0], 1 ^ vTwist[0],
                                wNext[1 ^ hasTwist],
                                vTwist[0] ^ (hasTwist ^ wTwist[1 ^ hasTwist]));
                            vtxBdryJoin(vNext[1], 0 ^ vTwist[1],
                                wNext[0 ^ hasTwist],
                                vTwist[1] ^ (hasTwist ^ wTwist[0 ^ hasTwist]));
                        }
                    }
                }

                vertexState[vIdx].bdryEdges--;
                vertexState[wIdx].bdryEdges--;
            }
        } else {
            // We are joining two distinct vertices together and merging
            // their vertex links.
            if (vertexState[vRep].rank < vertexState[wRep].rank) {
                // Join vRep beneath wRep.
                vertexState[vRep].parent = wRep;
                vertexState[vRep].twistUp = hasTwist ^ parentTwists;

                vertexState[wRep].bdry = vertexState[wRep].bdry +
                    vertexState[vRep].bdry - 2;
                if (vertexState[wRep].bdry == 0)
                    retVal |= VLINK_CLOSED;

                vertexStateChanged[orderIdx] = vRep;
            } else {
                // Join wRep beneath vRep.
                vertexState[wRep].parent = vRep;
                vertexState[wRep].twistUp = hasTwist ^ parentTwists;
                if (vertexState[vRep].rank == vertexState[wRep].rank) {
                    vertexState[vRep].rank++;
                    vertexState[wRep].hadEqualRank = true;
                }

                vertexState[vRep].bdry = vertexState[vRep].bdry +
                    vertexState[wRep].bdry - 2;
                if (vertexState[vRep].bdry == 0)
                    retVal |= VLINK_CLOSED;

                vertexStateChanged[orderIdx] = wRep;
            }

            nVertexClasses--;

            // Adjust the cycles of boundary components.
            if (vertexState[vIdx].bdryEdges == 2)
                vtxBdryBackup(vIdx);
            if (vertexState[wIdx].bdryEdges == 2)
                vtxBdryBackup(wIdx);

            if (vtxBdryLength1(vIdx)) {
                if (vtxBdryLength1(wIdx)) {
                    // Both vIdx and wIdx form entire boundary components of
                    // length one; these are joined together and the vertex
                    // link is closed off.
                    // No changes to make for the boundary cycles.
                } else {
                    // Here vIdx forms a boundary component of length one,
                    // and wIdx does not.  Ignore vIdx, and simply excise the
                    // relevant edge from wIdx.
                    // There is nothing to do here unless wIdx only has one
                    // boundary edge remaining (in which case we know it
                    // joins to some different tetrahedron vertex).
                    if (vertexState[wIdx].bdryEdges == 1) {
                        wNext[0] = vertexState[wIdx].bdryNext[0];
                        wNext[1] = vertexState[wIdx].bdryNext[1];
                        wTwist[0] = vertexState[wIdx].bdryTwist[0];
                        wTwist[1] = vertexState[wIdx].bdryTwist[1];

                        vtxBdryJoin(wNext[0], 1 ^ wTwist[0], wNext[1],
                            wTwist[0] ^ wTwist[1]);
                    }
                }
            } else if (vtxBdryLength1(wIdx)) {
                // As above, but with the two vertices the other way around.
                if (vertexState[vIdx].bdryEdges == 1) {
                    vNext[0] = vertexState[vIdx].bdryNext[0];
                    vNext[1] = vertexState[vIdx].bdryNext[1];
                    vTwist[0] = vertexState[vIdx].bdryTwist[0];
                    vTwist[1] = vertexState[vIdx].bdryTwist[1];

                    vtxBdryJoin(vNext[0], 1 ^ vTwist[0], vNext[1],
                        vTwist[0] ^ vTwist[1]);
                }
            } else {
                // Each vertex belongs to a boundary component of length
                // at least two.  Merge the components together.
                vtxBdryNext(vIdx, face.tet, v, face.face, vNext, vTwist);
                vtxBdryNext(wIdx, adj.tet, w, adj.face, wNext, wTwist);

                vtxBdryJoin(vNext[0], 1 ^ vTwist[0], wNext[1 ^ hasTwist],
                    vTwist[0] ^ (hasTwist ^ wTwist[1 ^ hasTwist]));
                vtxBdryJoin(vNext[1], 0 ^ vTwist[1], wNext[0 ^ hasTwist],
                    vTwist[1] ^ (hasTwist ^ wTwist[0 ^ hasTwist]));
            }

            vertexState[vIdx].bdryEdges--;
            vertexState[wIdx].bdryEdges--;
        }
    }

    return retVal;
}

void NClosedPrimeMinSearcher::splitVertexClasses() {
    // Split all three vertex pairs for the current face.
    NTetFace face = order[orderElt];
    NTetFace adj = (*pairing)[face];

    int v, w;
    int vIdx, wIdx;
    unsigned orderIdx;
    int rep, subRep;
    NPerm p = gluingPerm(face);
    // Do everything in reverse.  This includes the loop over vertices.
    for (v = 3; v >= 0; v--) {
        if (v == face.face)
            continue;

        w = p[v];
        vIdx = v + 4 * face.tet;
        wIdx = w + 4 * adj.tet;
        orderIdx = v + 4 * orderElt;

        if (vertexStateChanged[orderIdx] < 0) {
            for (rep = vIdx; vertexState[rep].parent >= 0;
                    rep = vertexState[rep].parent)
                ;
            vertexState[rep].bdry += 2;
        } else {
            subRep = vertexStateChanged[orderIdx];
            rep = vertexState[subRep].parent;

            vertexState[subRep].parent = -1;
            if (vertexState[subRep].hadEqualRank) {
                vertexState[subRep].hadEqualRank = false;
                vertexState[rep].rank--;
            }

            vertexState[rep].bdry = vertexState[rep].bdry + 2 -
                vertexState[subRep].bdry;

            vertexStateChanged[orderIdx] = -1;
            nVertexClasses++;
        }

        // Restore cycles of boundary components.
        if (vIdx == wIdx) {
            // We did nothing during the merge; do nothing during the split.
        } else {
            vertexState[wIdx].bdryEdges++;
            vertexState[vIdx].bdryEdges++;

            switch (vertexState[wIdx].bdryEdges) {
                case 3: vertexState[wIdx].bdryNext[0] =
                            vertexState[wIdx].bdryNext[1] = wIdx;
                        vertexState[wIdx].bdryTwist[0] =
                            vertexState[wIdx].bdryTwist[1] = 0;
                        break;

                case 2: vtxBdryRestore(wIdx);
                        // Fall through to the next case, so we can
                        // adjust the neighbours.

                case 1: // Nothing was changed for wIdx during the merge,
                        // so there is nothing there to restore.

                        // Adjust neighbours to point back to wIdx.
                        vtxBdryFixAdj(wIdx);
            }

            switch (vertexState[vIdx].bdryEdges) {
                case 3: vertexState[vIdx].bdryNext[0] =
                            vertexState[vIdx].bdryNext[1] = vIdx;
                        vertexState[vIdx].bdryTwist[0] =
                            vertexState[vIdx].bdryTwist[1] = 0;
                        break;

                case 2: vtxBdryRestore(vIdx);
                        // Fall through to the next case, so we can
                        // adjust the neighbours.

                case 1: // Nothing was changed for vIdx during the merge,
                        // so there is nothing there to restore.

                        // Adjust neighbours to point back to vIdx.
                        vtxBdryFixAdj(vIdx);
            }
        }
    }
}

int NClosedPrimeMinSearcher::mergeEdgeClasses() {
    NTetFace face = order[orderElt];
    NTetFace adj = (*pairing)[face];

    int retVal = 0;

    NPerm p = gluingPerm(face);
    int v1, w1, v2, w2;
    int e, f;
    int orderIdx;
    int eRep, fRep;
    int middleTet;

    v1 = face.face;
    w1 = p[v1];

    char parentTwists, hasTwist;
    for (v2 = 0; v2 < 4; v2++) {
        if (v2 == v1)
            continue;

        w2 = p[v2];

        // Look at the edge opposite v1-v2.
        e = 5 - edgeNumber[v1][v2];
        f = 5 - edgeNumber[w1][w2];

        orderIdx = v2 + 4 * orderElt;

        // We declare the natural orientation of an edge to be smaller
        // vertex to larger vertex.
        hasTwist = (p[edgeStart[e]] > p[edgeEnd[e]] ? 1 : 0);

        parentTwists = 0;
        eRep = findEdgeClass(e + 6 * face.tet, parentTwists);
        fRep = findEdgeClass(f + 6 * adj.tet, parentTwists);

        if (eRep == fRep) {
            edgeState[eRep].bounded = false;

            if (edgeState[eRep].size <= 2)
                retVal |= ECLASS_LOWDEG;
            else if (edgeState[eRep].size == 3) {
                // Flag as LOWDEG only if three distinct tetrahedra are used.
                middleTet = pairing->dest(face.tet, v2).tet;
                if (face.tet != adj.tet && adj.tet != middleTet &&
                        middleTet != face.tet)
                    retVal |= ECLASS_LOWDEG;
            }
            if (hasTwist ^ parentTwists)
                retVal |= ECLASS_TWISTED;

            edgeStateChanged[orderIdx] = -1;
        } else {
#if PRUNE_HIGH_DEG_EDGE_SET
            if (edgeState[eRep].size >= 3) {
                if (edgeState[fRep].size >= 3)
                    highDegSum += 3;
                else
                    highDegSum += edgeState[fRep].size;
            } else if (edgeState[fRep].size >= 3)
                highDegSum += edgeState[eRep].size;
            else if (edgeState[eRep].size == 2 && edgeState[fRep].size == 2)
                ++highDegSum;
#endif

            if (edgeState[eRep].rank < edgeState[fRep].rank) {
                // Join eRep beneath fRep.
                edgeState[eRep].parent = fRep;
                edgeState[eRep].twistUp = hasTwist ^ parentTwists;

                edgeState[fRep].size += edgeState[eRep].size;
#if PRUNE_HIGH_DEG_EDGE_SET
#else
                if (edgeState[fRep].size > 3 * getNumberOfTetrahedra())
                    retVal |= ECLASS_HIGHDEG;
#endif

                edgeStateChanged[orderIdx] = eRep;
            } else {
                // Join fRep beneath eRep.
                edgeState[fRep].parent = eRep;
                edgeState[fRep].twistUp = hasTwist ^ parentTwists;
                if (edgeState[eRep].rank == edgeState[fRep].rank) {
                    edgeState[eRep].rank++;
                    edgeState[fRep].hadEqualRank = true;
                }

                edgeState[eRep].size += edgeState[fRep].size;
#if PRUNE_HIGH_DEG_EDGE_SET
#else
                if (edgeState[eRep].size > 3 * getNumberOfTetrahedra())
                    retVal |= ECLASS_HIGHDEG;
#endif

                edgeStateChanged[orderIdx] = fRep;
            }

#if PRUNE_HIGH_DEG_EDGE_SET
            if (highDegSum > highDegBound)
                retVal |= ECLASS_HIGHDEG;
#endif

            nEdgeClasses--;
        }
    }

    // If we've already found something bad, exit now.  No sense in
    // looking for even more bad structures, since we're only going to
    // discard the triangulation anyway.
    if (retVal)
        return retVal;

    // Find representatives of the equivalence classes for all six edges
    // of the current tetrahedron (instead of calculating them each time
    // we want them).
    int tRep[6];
    char tTwist[6];
    for (e = 0; e < 6; e++)
        tRep[e] = findEdgeClass(e + 6 * face.tet, tTwist[e] = 0);

    // Test for cones in all possible positions on all possible faces.
    // Apologies for the tightness of the code; this part is being
    // micro-optimised since it is run so very frequently.  The old,
    // more readable version of this code is in the commented block below.
    for (e = 0; e < 12; e++)
        if (tRep[coneEdge[e][0]] == tRep[coneEdge[e][1]] && (coneNoTwist[e] ^
            (tTwist[coneEdge[e][0]] ^ tTwist[coneEdge[e][1]])))
                return ECLASS_CONE;

    /*
    // Test for cones on edges v1->w1->v2.
    for (w1 = 0; w1 < 4; w1++)
        for (v1 = 0; v1 < 3; v1++) {
            if (v1 == w1)
                continue;
            for (v2 = v1 + 1; v2 < 4; v2++) {
                if (v2 == w1)
                    continue;

                parentTwists = tTwist[edgeNumber[v1][w1]] ^
                    tTwist[edgeNumber[v2][w1]];

                if (tRep[edgeNumber[v1][w1]] == tRep[edgeNumber[v2][w1]]) {
                    hasTwist = (v1 < w1 && w1 < v2 ? 0 : 1);
                    if (hasTwist ^ parentTwists) {
                        return ECLASS_CONE;
                    }
                }
            }
        }
    */

    // Test for L(3,1) spines.
    // Don't bother checking the directions of the edges -- if it's not an
    // L(3,1) spine then it includes a cone, which we've already tested for.

    // L(3,1) on face 012:
    if (tRep[0] == tRep[1] && tRep[1] == tRep[3])
        return ECLASS_L31;
    // L(3,1) on face 013:
    if (tRep[0] == tRep[2] && tRep[2] == tRep[4])
        return ECLASS_L31;
    // L(3,1) on face 023:
    if (tRep[1] == tRep[2] && tRep[2] == tRep[5])
        return ECLASS_L31;
    // L(3,1) on face 123:
    if (tRep[3] == tRep[4] && tRep[4] == tRep[5])
        return ECLASS_L31;

    // Nothing bad was found.
    return 0;
}

void NClosedPrimeMinSearcher::splitEdgeClasses() {
    NTetFace face = order[orderElt];

    int v1, v2;
    int e;
    int eIdx, orderIdx;
    int rep, subRep;

    v1 = face.face;

    for (v2 = 3; v2 >= 0; v2--) {
        if (v2 == v1)
            continue;

        // Look at the edge opposite v1-v2.
        e = 5 - edgeNumber[v1][v2];

        eIdx = e + 6 * face.tet;
        orderIdx = v2 + 4 * orderElt;

        if (edgeStateChanged[orderIdx] < 0)
            edgeState[findEdgeClass(eIdx)].bounded = true;
        else {
            subRep = edgeStateChanged[orderIdx];
            rep = edgeState[subRep].parent;

            edgeState[subRep].parent = -1;
            if (edgeState[subRep].hadEqualRank) {
                edgeState[subRep].hadEqualRank = false;
                edgeState[rep].rank--;
            }

            edgeState[rep].size -= edgeState[subRep].size;
#if PRUNE_HIGH_DEG_EDGE_SET
            if (edgeState[rep].size >= 3) {
                if (edgeState[subRep].size >= 3)
                    highDegSum -= 3;
                else
                    highDegSum -= edgeState[subRep].size;
            } else if (edgeState[subRep].size >= 3)
                highDegSum -= edgeState[rep].size;
            else if (edgeState[rep].size == 2 && edgeState[subRep].size == 2)
                --highDegSum;
#endif

            edgeStateChanged[orderIdx] = -1;
            nEdgeClasses++;
        }
    }
}

void NClosedPrimeMinSearcher::vtxBdryConsistencyCheck() {
    int adj, id, end;
    for (id = 0; id < static_cast<int>(getNumberOfTetrahedra()) * 4; id++)
        if (vertexState[id].bdryEdges > 0)
            for (end = 0; end < 2; end++) {
                adj = vertexState[id].bdryNext[end];
                if (vertexState[adj].bdryEdges == 0)
                    std::cerr << "CONSISTENCY ERROR: Vertex link boundary "
                        << id << '/' << end
                        << " runs into an internal vertex." << std::endl;
                if (vertexState[adj].bdryNext[(1 ^ end) ^
                        vertexState[id].bdryTwist[end]] != id)
                    std::cerr << "CONSISTENCY ERROR: Vertex link boundary "
                        << id << '/' << end
                        << " has a mismatched adjacency." << std::endl;
                if (vertexState[adj].bdryTwist[(1 ^ end) ^
                        vertexState[id].bdryTwist[end]] !=
                        vertexState[id].bdryTwist[end])
                    std::cerr << "CONSISTENCY ERROR: Vertex link boundary "
                        << id << '/' << end
                        << " has a mismatched twist." << std::endl;
            }
}

void NClosedPrimeMinSearcher::vtxBdryDump(std::ostream& out) {
    for (unsigned id = 0; id < getNumberOfTetrahedra() * 4; id++) {
        if (id > 0)
            out << ' ';
        out << vertexState[id].bdryNext[0]
            << (vertexState[id].bdryTwist[0] ? '~' : '-')
            << id
            << (vertexState[id].bdryTwist[1] ? '~' : '-')
            << vertexState[id].bdryNext[1];
    }
    out << std::endl;
}

} // namespace regina
