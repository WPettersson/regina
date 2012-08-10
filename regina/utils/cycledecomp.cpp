
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2011, Ben Burton                                   *
 *                2012-2012, William Pettersson                           *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                  or William Pettersson (william.pettersson@gmail.com). *
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

using namespace regina;

namespace CycleDecompSearcher {

const unsigned faceEdges[4][3] = {
    {3,4,5}, {1,2,5}, {0,2,4}, {0,1,3}
};

const unsigned otherFace[4][3] = {
    {3,2,1}, {3,2,0}, {3,1,0}, {2,1,0}
};

const signed edgeParity[6][6] = {
    {-1,2,1,4,3,-1}, {2,-1,0,5,-1,3}, {1,0,-1,-1,5,4},
    {4,5,-1,-1,0,1}, {3,-1,5,0,-1,2}, {-1,3,4,1,2,-1}
};

const char dataTag_ = 'd';



CycleDecompSearcher(const NFacePairing* pairing,
        const NFacePairing::IsoList* autos, bool orientableOnly,
        UseGluingPerms use, void* useArgs) :
        NCompactSearcher(pairing, autos, orientableOnly,
            NCensus::PURGE_NON_MINIMAL_PRIME | NCensus::PURGE_P2_REDUCIBLE,
            use, useArgs) {

    nTets = size();
    nEdges = 2*nTets;
    nEnds = 2*nEdges;
    
    tets = new Tetrahedron[nTets];
    edges = new Edge[nEdges];
    ends = new EdgeEnds[nEnds];
    nextColour = 0;


    bool* orderAssigned = new bool[nTets * 4];
        /**< Have we placed a tetrahedron face or its partner in the
             order[] array yet? */

    // Hunt for structures within the face pairing graph.

    NTetFace face, adj;
    unsigned edgesDone = 0;
    std::fill(orderAssigned, orderAssigned + 4 * nTets, false);

    // Begin by searching for tetrahedra that are joined to themselves.
    // Note that each tetrahedra can be joined to itself at most once,
    // since we are guaranteed that the face pairing is connected with
    // order >= 3.

    for (face.setFirst(); ! face.isPastEnd(nTets, true); face++) {
        if (orderAssigned[face.simp * 4 + face.facet])
            continue;

        adj = (*pairing)[face];
        if (adj.simp != face.simp)
            continue;

        ends[2*edgesDone].tet = &test[face.simp];
        ends[2*edgesDone+1].tet = &test[adj.simp];
        ends[2*edgesDone].face = face.facet;
        ends[2*edgesDone+1].face = adj.facet;
        ends[2*edgesDone].edge = &edges[edgesDone];
        ends[2*edgesDone+1].edge = &edges[edgesDone];


        edges[edgesDone].ends[0]=&ends[2*edgesDone];
        edges[edgesDone].ends[1]=&ends[2*edgesDone+1];
        edges[edgesDone].index = edgesDone;
                
        tets[face.simp].externalEdgeEnd[face.facet] = 
            edges[edgesDone].ends[0];

        tets[adj.simp].externalEdgeEnd[adj.facet] = 
            edges[edgesDone].ends[1];

        orderAssigned[face.simp * 4 + face.facet] = true;
        orderAssigned[adj.simp * 4 + adj.facet] = true;
        edgesDone++;
    }

    // All done for the orderAssigned[] array.  Tidy up.
    delete[] orderAssigned;

}

void colourOnTetrahedra(tet) {
    unsigned edge;
    // Note empty for loop to find first unused internal edge.
    for (edge=0; tets[tet].internalEdges[edge] != 0; edge++) ;

    if (edge >= 6)
        return;
    nextColour++;
    tets[tet].internalEdges[edge] = nextColour;
    tets[tet].used++;
    
    unsigned startFace = NEdge::edgeVertex[5-edge][0];
    EdgeEnd *start = tets[tet].externalEdgeEnd[startFace];
    unsigned outFace = NEdge::edgeVertex[5-edge][1];
    unsigned nextEdgeEnd = tets[tet].externalEdgeEnd[outFace];
    Edge *nextEdge = ends[nextEdgeEnd].edge;
    ends[nextEdgeEnd].map[edge] = nextEdge->used;
    nextEdge->colour(nextColour);
    edgesLeft--;
    EdgeEnd *nextEdgeEnd = nextEdge->otherEnd(ends[nextEdgeEnd]);

    findPath(start, edge, nextEdgeEnd);
    
    nextEdge->unColour(); 
    ends[nextEdgeEnd].map[edge] = 0;
    edgesLeft++;
    tets[tet].internalEdges[edge] = 0;
    tets[tet].used--;
    nextColour--;
}

bool checkColourOk() {
    // No short cycles
    if ( cycleLength[nextColour] < 3) 
        return false;
    // A length 3 cycle on 3 distinct tetrahedra implies a 3-2 move.
    // Therefore one of the edges must be a loop.
    if ( cycleLength[nextColour] == 3 ) {
        bool foundLoop = false;
        for(unsigned i=0;i<3;i++) { 
            Edge *e = edges[cycles[nextColour][i]];
            if (e->ends[0]->tet == e->ends[1]->tet) {
                foundLoop = true;
                break;
            }
        }
        if (foundLoop == false) 
            return false;
    }

    return true;
}

void nextPath(EdgeEnd *start, unsigned firstEdge,  EdgeEnd *now) {
    unsigned nextInternal;
    EdgeEnd *nextEnd;
    Edge *nextEdge;
    Tetrahedron *nextTet = now->tet;
    // Count the number of edges left. We need nTet+1 cycles,
    // and are working on cycle number nextColour. We need 3 edges
    // for each of the (nTet+1 - nextColour) cycles after this one.
    if (edgesLeft < 3*(nTet+1 - nextColour)) {
        return;
    }
    for (unsigned i=0; i<3 ; i++) {
        nextInternal = faceEdges[now->face][i];
        if (nextTet->internalEdges[nextInternal] != 0) 
            continue;
        nextEnd = nextTet->externalEdgeEnd[otherface[now->face][nextInternal]];
        nextEdge = nextEnd->edge;
        
        // An edge has two ends, ends[0] and ends[1]
        // A "forward" direction along an edge is from ends[0] to ends[1]
        // and "backwards" is the reverse. Forwards is denoted by a positive
        // edge number and backwards by a negative edge number.
        signed dir = nextEdge->index;
        if (nextEnd == (nextEnd->edge->ends[1]) ) 
            dir = -1* nextEdge->index;
        // In the following loop, dir=0 means the edge has already been used in
        // the opposite direction.
        for(unsigned j=0; j< cycleLength[nextColour];j++) {
            if (cycles[nextColour][j] == -1*dir) {
                dir=0;
                break;  
            }
        }
        if (dir == 0) 
            continue;
        cycles[nextColour][cycleLength[nextColour]]=dir;
        cycleLength[nextColour]++;
        
        nextTet->internalEdges[nextInternal] = nextColour;
        nextEdge->colour(nextColour);
        now->map[nextInternal] = nextEdge->used;
        edgesLeft--;
        
        if (nextEnd == start) {
            start->map[firstEdge] = start->edge->used;
            if (checkColourOk()) {
                colourOnTetrahedra(findTetWithMostInternalEdgesUsed());
            }
            start->map[firstEdge] = 0;
        }

        nextPath(start, firstEdge, nextEnd);


        edgesLeft++;
        now->map[nextInternal]= 0;
        nextEdge->unColour(nextColour);
        nextTet->internalEdges[nextInternal] = 0;
}
 

bool checkComplete() {
    for(unsigned i=0; i< nEdges; i++) {
        if (edges[i].used != 3) 
            return false;
    }
}

unsigned findTetWithMostInternalEdgesUsed() {
    unsigned mostUsed=0;
    // If nothing is used, start on tetrahedra 0.
    unsigned tet=0;
    for(unsigned i=0; i< nTet; i++) {
        if (tets[i].used > mostUsed) {
            mostUsed = tets[i].used;
            tet = i;
        }
    }
    return tet;
}

void runSearch(long maxDepth) {
    // Preconditions:
    //     Only closed prime minimal P2-irreducible triangulations are needed.
    //     The given face pairing is closed with order >= 3.

    if (! started) {
        // Search initialisation.
        started = true;

        // Begin by testing for face pairings that can never lead to such a
        // triangulation.
        //if (pairing_->hasTripleEdge() ||
        //        pairing_->hasBrokenDoubleEndedChain() ||
        //        pairing_->hasOneEndedChainWithDoubleHandle() ||
        //        pairing_->hasOneEndedChainWithStrayBigon() ||
        //        pairing_->hasWedgedDoubleEndedChain() ||
        //        pairing_->hasTripleOneEndedChain()) {
        //    use_(0, useArgs_);
        //    return;
        //}

        orderElt = 0;
    }

    // Is it a partial search that has already finished?
    if (orderElt == static_cast<int>(nTets) * 2) {
        if (isCanonical())
            use_(this, useArgs_);
        use_(0, useArgs_);
        return;
    }

    unsigned tet = findTetWithMostInternalEdgesUsed();
    colourOnTetrahedra(tet);

    // ---------- Selecting the individual gluing permutations ----------

    // Observe that in a canonical face pairing, one-ended chains always
    // follow an increasing sequence of tetrahedra from boundary to end,
    // or follow the sequence of tetrahedra 0, 1, ..., k from end to
    // boundary.
    //
    // In particular, this means that for any tetrahedron not internal
    // to a one-ended chain (with the possible exception of tetrahedron
    // order[nChainEdges].simp), face 0 of this tetrahedron is not
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
        adj = (*pairing_)[face];

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
            if ((! orientableOnly_) || adj.facet == 0)
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
        permIndex(adj) = NPerm4::invS3[permIndex(face)];

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
        if (generic && adj.facet == 0 && orientableOnly_) {
            // It's the first time we've hit this tetrahedron.
            if ((permIndex(face) + (face.facet == 3 ? 0 : 1) +
                    (adj.facet == 3 ? 0 : 1)) % 2 == 0)
                orientation[adj.simp] = -orientation[face.simp];
            else
                orientation[adj.simp] = orientation[face.simp];
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
            if (orientableOnly_ && pairing_->dest(face).facet > 0) {
                // permIndex(face) will be set to -1 or -2 as appropriate.
                adj = (*pairing_)[face];
                if (orientation[face.simp] == orientation[adj.simp])
                    permIndex(face) = 1;
                else
                    permIndex(face) = 0;

                if ((face.facet == 3 ? 0 : 1) + (adj.facet == 3 ? 0 : 1) == 1)
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
                    "true at end of search!" << std::endl;
            if (vertexState[i].bdryTwist[1])
                std::cerr << "ERROR: vertexState[" << i << "].bdryTwist == "
                    "true at end of search!" << std::endl;
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
                    "false at end of search!" << std::endl;
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
    NCompactSearcher::dumpData(out);

    int i;
    for (i = 0; i < orderSize; i++) {
        if (i)
            out << ' ';
        out << orderType[i];
    }
    out << std::endl;

    out << nChainEdges << std::endl;
    if (nChainEdges) {
        for (i = 0; i < 2 * static_cast<int>(nChainEdges); i++) {
            if (i)
                out << ' ';
            out << chainPermIndices[i];
        }
        out << std::endl;
    }
}

NClosedPrimeMinSearcher::NClosedPrimeMinSearcher(std::istream& in,
        UseGluingPerms use, void* useArgs) :
        NCompactSearcher(in, use, useArgs),
        orderType(0), nChainEdges(0), chainPermIndices(0) {
    if (inputError_)
        return;

    unsigned nTets = getNumberOfTetrahedra();
    int i;

    orderType = new unsigned[2 * nTets];
    for (i = 0; i < orderSize; i++)
        in >> orderType[i];

    in >> nChainEdges;
    /* Unnecessary since nChainEdges is unsigned.
    if (nChainEdges < 0) {
        inputError_ = true; return;
    } */
    if (nChainEdges) {
        chainPermIndices = new int[nChainEdges * 2];
        for (i = 0; i < 2 * static_cast<int>(nChainEdges); i++) {
            in >> chainPermIndices[i];
            if (chainPermIndices[i] < 0 || chainPermIndices[i] >= 6) {
                inputError_ = true; return;
            }
        }
    }

#if PRUNE_HIGH_DEG_EDGE_SET
    in >> highDegLimit >> highDegSum >> highDegBound;
    if (highDegLimit < 3 || highDegLimit > 4 || highDegSum < 0 ||
            highDegSum > 6 * static_cast<int>(nTets) || highDegBound !=
                (6 - highDegLimit) * static_cast<int>(nTets) - highDegLimit) {
        inputError_ = true; return;
    }
#endif

    // Did we hit an unexpected EOF?
    if (in.eof())
        inputError_ = true;
}

int NClosedPrimeMinSearcher::mergeEdgeClasses() {
    NTetFace face = order[orderElt];
    NTetFace adj = (*pairing_)[face];

    int retVal = 0;

    NPerm4 p = gluingPerm(face);
    int v1, w1, v2, w2;
    int e, f;
    int orderIdx;
    int eRep, fRep;
    int middleTet;

    v1 = face.facet;
    w1 = p[v1];

    char parentTwists, hasTwist;
    for (v2 = 0; v2 < 4; v2++) {
        if (v2 == v1)
            continue;

        w2 = p[v2];

        // Look at the edge opposite v1-v2.
        e = 5 - NEdge::edgeNumber[v1][v2];
        f = 5 - NEdge::edgeNumber[w1][w2];

        orderIdx = v2 + 4 * orderElt;

        // We declare the natural orientation of an edge to be smaller
        // vertex to larger vertex.
        hasTwist = (p[NEdge::edgeVertex[e][0]] > p[NEdge::edgeVertex[e][1]] ?
            1 : 0);

        parentTwists = 0;
        eRep = findEdgeClass(e + 6 * face.simp, parentTwists);
        fRep = findEdgeClass(f + 6 * adj.simp, parentTwists);

        if (eRep == fRep) {
            edgeState[eRep].bounded = false;

            if (edgeState[eRep].size <= 2)
                retVal |= ECLASS_LOWDEG;
            else if (edgeState[eRep].size == 3) {
                // Flag as LOWDEG only if three distinct tetrahedra are used.
                middleTet = pairing_->dest(face.simp, v2).simp;
                if (face.simp != adj.simp && adj.simp != middleTet &&
                        middleTet != face.simp)
                    retVal |= ECLASS_LOWDEG;
            }
            if (hasTwist ^ parentTwists)
                retVal |= ECLASS_TWISTED;

            edgeStateChanged[orderIdx] = -1;
        } else {
#if PRUNE_HIGH_DEG_EDGE_SET
            if (edgeState[eRep].size >= highDegLimit) {
                if (edgeState[fRep].size >= highDegLimit)
                    highDegSum += highDegLimit;
                else
                    highDegSum += edgeState[fRep].size;
            } else if (edgeState[fRep].size >= highDegLimit)
                highDegSum += edgeState[eRep].size;
            else if (edgeState[eRep].size + edgeState[fRep].size >
                    highDegLimit)
                highDegSum += (edgeState[eRep].size + edgeState[fRep].size -
                    highDegLimit);
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

                if (edgeState[eRep].twistUp) {
                    edgeState[fRep].facesPos += edgeState[eRep].facesNeg;
                    edgeState[fRep].facesNeg += edgeState[eRep].facesPos;
                } else {
                    edgeState[fRep].facesPos += edgeState[eRep].facesPos;
                    edgeState[fRep].facesNeg += edgeState[eRep].facesNeg;
                }
                if (edgeState[fRep].facesPos.hasNonZeroMatch(
                        edgeState[fRep].facesNeg))
                    retVal |= ECLASS_CONE;
                if (edgeState[fRep].facesPos.has3() ||
                        edgeState[fRep].facesNeg.has3())
                    retVal |= ECLASS_L31;

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

                if (edgeState[fRep].twistUp) {
                    edgeState[eRep].facesPos += edgeState[fRep].facesNeg;
                    edgeState[eRep].facesNeg += edgeState[fRep].facesPos;
                } else {
                    edgeState[eRep].facesPos += edgeState[fRep].facesPos;
                    edgeState[eRep].facesNeg += edgeState[fRep].facesNeg;
                }
                if (edgeState[eRep].facesPos.hasNonZeroMatch(
                        edgeState[eRep].facesNeg))
                    retVal |= ECLASS_CONE;
                if (edgeState[eRep].facesPos.has3() ||
                        edgeState[eRep].facesNeg.has3())
                    retVal |= ECLASS_L31;

                edgeStateChanged[orderIdx] = fRep;
            }

#if PRUNE_HIGH_DEG_EDGE_SET
            if (highDegSum > highDegBound)
                retVal |= ECLASS_HIGHDEG;
#endif

            nEdgeClasses--;
        }
    }

    return retVal;
}

void NClosedPrimeMinSearcher::splitEdgeClasses() {
    NTetFace face = order[orderElt];

    int v1, v2;
    int e;
    int eIdx, orderIdx;
    int rep, subRep;

    v1 = face.facet;

    for (v2 = 3; v2 >= 0; v2--) {
        if (v2 == v1)
            continue;

        // Look at the edge opposite v1-v2.
        e = 5 - NEdge::edgeNumber[v1][v2];

        eIdx = e + 6 * face.simp;
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
            if (edgeState[rep].size >= highDegLimit) {
                if (edgeState[subRep].size >= highDegLimit)
                    highDegSum -= highDegLimit;
                else
                    highDegSum -= edgeState[subRep].size;
            } else if (edgeState[subRep].size >= highDegLimit)
                highDegSum -= edgeState[rep].size;
            else if (edgeState[rep].size + edgeState[subRep].size >
                    highDegLimit)
                highDegSum -= (edgeState[rep].size + edgeState[subRep].size
                    - highDegLimit);
#endif

            if (edgeState[subRep].twistUp) {
                edgeState[rep].facesPos -= edgeState[subRep].facesNeg;
                edgeState[rep].facesNeg -= edgeState[subRep].facesPos;
            } else {
                edgeState[rep].facesPos -= edgeState[subRep].facesPos;
                edgeState[rep].facesNeg -= edgeState[subRep].facesNeg;
            }

            edgeStateChanged[orderIdx] = -1;
            nEdgeClasses++;
        }
    }
}

} // namespace regina

