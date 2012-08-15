
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

#include <assert.h>
#include <sstream>
#include "census/ncensus.h"
#include "census/ngluingpermsearcher.h"
#include "triangulation/nedge.h"
#include "triangulation/nfacepair.h"
#include "triangulation/ntriangulation.h"
#include "utilities/boostutils.h"
#include "utilities/memutils.h"

#include "cycledecomp.h"

using namespace regina;

const unsigned CycleDecompSearcher::faceEdges[4][3] = {
    {3,4,5}, {1,2,5}, {0,2,4}, {0,1,3}
};

const unsigned CycleDecompSearcher::otherFace[4][3] = {
    {3,2,1}, {3,2,0}, {3,1,0}, {2,1,0}
};

const signed CycleDecompSearcher::edgeParity[6][6] = {
    {-1,2,1,4,3,-1}, {2,-1,0,5,-1,3}, {1,0,-1,-1,5,4},
    {4,5,-1,-1,0,1}, {3,-1,5,0,-1,2}, {-1,3,4,1,2,-1}
};

const char CycleDecompSearcher::dataTag_ = 'd';



CycleDecompSearcher::CycleDecompSearcher(const NFacePairing* pairing,
        const NFacePairing::IsoList* autos, bool orientableOnly,
        UseCycles use, void* useArgs) {

    use_ = use;
    useArgs_ = useArgs;

    nTets = pairing->size();
    nEdges = 2*nTets;
    nEnds = 2*nEdges;

    edgesLeft = nEdges;
    
    tets = new Tetrahedron[nTets];
    edges = new Edge[nEdges];
    ends = new EdgeEnd[nEnds];
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

        ends[2*edgesDone].tet = &tets[face.simp];
        ends[2*edgesDone+1].tet = &tets[adj.simp];
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

CycleDecompSearcher::~CycleDecompSearcher() {
    delete[] tets;
    delete[] edges;
    delete[] ends;
}

void CycleDecompSearcher::colourOnTetrahedra(unsigned tet) {
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
    EdgeEnd *outEdgeEnd = tets[tet].externalEdgeEnd[outFace];
    Edge *nextEdge = outEdgeEnd->edge;
    outEdgeEnd->map[edge] = nextEdge->used;
    nextEdge->colour(nextColour);
    edgesLeft--;
    EdgeEnd *nextEdgeEnd = nextEdge->otherEnd(outEdgeEnd);

    nextPath(start, edge, nextEdgeEnd);
    
    nextEdge->unColour(); 
    outEdgeEnd->map[edge] = 0;
    edgesLeft++;
    tets[tet].internalEdges[edge] = 0;
    tets[tet].used--;
    nextColour--;
}

bool CycleDecompSearcher::checkColourOk() {
    // No short cycles
    if ( cycleLengths[nextColour] < 3) 
        return false;
    // A length 3 cycle on 3 distinct tetrahedra implies a 3-2 move.
    // Therefore one of the edges must be a loop.
    if ( cycleLengths[nextColour] == 3 ) {
        bool foundLoop = false;
        for(unsigned i=0;i<3;i++) { 
            Edge *e = &(edges[cycles[nextColour][i]]);
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

void CycleDecompSearcher::nextPath(EdgeEnd *start, unsigned firstEdge,  EdgeEnd *now) {
    unsigned nextInternal;
    EdgeEnd *nextEnd;
    Edge *nextEdge;
    Tetrahedron *nextTet = now->tet;
    // Count the number of edges left. We need nTet+1 cycles,
    // and are working on cycle number nextColour. We need 3 edges
    // for each of the (nTet+1 - nextColour) cycles after this one.
    if (edgesLeft < 3*(nTets+1 - nextColour)) {
        return;
    }
    for (unsigned i=0; i<3 ; i++) {
        nextInternal = faceEdges[now->face][i];
        if (nextTet->internalEdges[nextInternal] != 0) 
            continue;
        nextEnd = nextTet->externalEdgeEnd[otherFace[now->face][nextInternal]];
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
        for(unsigned j=0; j< cycleLengths[nextColour];j++) {
            if (cycles[nextColour][j] == -1*dir) {
                dir=0;
                break;  
            }
        }
        if (dir == 0) 
            continue;
        cycles[nextColour][cycleLengths[nextColour]]=dir;
        cycleLengths[nextColour]++;
        
        nextTet->internalEdges[nextInternal] = nextColour;
        nextEdge->colour(nextColour);
        now->map[nextInternal] = nextEdge->used;
        edgesLeft--;
       
        // Try to complete the cycle
        if (nextEnd == start) {
            start->map[firstEdge] = start->edge->used;
            if (checkColourOk()) {
                if (checkComplete()) {
                    use_(this, useArgs_);
                } else {
                    colourOnTetrahedra(findTetWithMostInternalEdgesUsed());
                }
            }
            start->map[firstEdge] = 0;
        }

        nextPath(start, firstEdge, nextEnd);


        edgesLeft++;
        now->map[nextInternal]= 0;
        nextEdge->unColour();
        nextTet->internalEdges[nextInternal] = 0;
    }
}
 

bool CycleDecompSearcher::checkComplete() {
    for(unsigned i=0; i< nEdges; i++) {
        if (edges[i].used != 3) 
            return false;
    }
}

unsigned CycleDecompSearcher::findTetWithMostInternalEdgesUsed() {
    unsigned mostUsed=0;
    // If nothing is used, start on tetrahedra 0.
    unsigned tet=0;
    for(unsigned i=0; i< nTets; i++) {
        if (tets[i].used > mostUsed) {
            mostUsed = tets[i].used;
            tet = i;
        }
    }
    return tet;
}

void CycleDecompSearcher::runSearch(long maxDepth) {
    // Preconditions:
    //     Only closed prime minimal P2-irreducible triangulations are needed.
    //     The given face pairing is closed with order >= 3.

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


    unsigned tet = findTetWithMostInternalEdgesUsed();
    colourOnTetrahedra(tet);


    use_(0, useArgs_);
}

NTriangulation* CycleDecompSearcher::triangulate() {
    NTriangulation* ans = new NTriangulation;
    NTetrahedron** simp = new NTetrahedron*[nTets];
    unsigned t;
    for (t = 0; t < nTets; ++t)
        simp[t] = ans->newSimplex();



}

void CycleDecompSearcher::dumpData(std::ostream& out) const {
    //NCompactSearcher::dumpData(out);

    unsigned i;
    unsigned j;
    for (i = 0; i < nTets+1; i++) {
        out << i << ": ";
        for(j = 0; j < cycleLengths[i]; j++) {
            out << cycles[i][j];
            if (j != cycleLengths[i]-1) 
                out << ", ";
        out << std::endl;
        }
    }
}

inline void CycleDecompSearcher::Edge::colour(unsigned newColour) {
    assert(used < 3);
    colours[used++] = newColour;
}

inline void CycleDecompSearcher::Edge::unColour() {
    assert(used > 0);
    colours[used--] = 0;
}

inline CycleDecompSearcher::EdgeEnd* CycleDecompSearcher::Edge::otherEnd(EdgeEnd *one) {
    if ( ends[0] == one )
        return ends[1];
    else {
        assert(one == ends[1]);
        return ends[0];
    }
}


//ClosedPrimeMinSearcher(std::istream& in,
//        UseGluingPerms use, void* useArgs) :
//        NCompactSearcher(in, use, useArgs),
//        orderType(0), nChainEdges(0), chainPermIndices(0) {
//    if (inputError_)
//        return;
//
//    unsigned nTets = getNumberOfTetrahedra();
//    int i;
//
//    orderType = new unsigned[2 * nTets];
//    for (i = 0; i < orderSize; i++)
//        in >> orderType[i];
//
//    in >> nChainEdges;
//    /* Unnecessary since nChainEdges is unsigned.
//    if (nChainEdges < 0) {
//        inputError_ = true; return;
//    } */
//    if (nChainEdges) {
//        chainPermIndices = new int[nChainEdges * 2];
//        for (i = 0; i < 2 * static_cast<int>(nChainEdges); i++) {
//            in >> chainPermIndices[i];
//            if (chainPermIndices[i] < 0 || chainPermIndices[i] >= 6) {
//                inputError_ = true; return;
//            }
//        }
//    }
//
//#if PRUNE_HIGH_DEG_EDGE_SET
//    in >> highDegLimit >> highDegSum >> highDegBound;
//    if (highDegLimit < 3 || highDegLimit > 4 || highDegSum < 0 ||
//            highDegSum > 6 * static_cast<int>(nTets) || highDegBound !=
//                (6 - highDegLimit) * static_cast<int>(nTets) - highDegLimit) {
//        inputError_ = true; return;
//    }
//#endif
//
//    // Did we hit an unexpected EOF?
//    if (in.eof())
//        inputError_ = true;
//}


