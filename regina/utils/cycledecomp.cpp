
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
#include <stdlib.h>

#include "census/ncensus.h"
#include "census/ngluingpermsearcher.h"
#include "triangulation/nedge.h"
#include "triangulation/nfacepair.h"
#include "triangulation/ntriangulation.h"
#include "utilities/boostutils.h"
#include "utilities/memutils.h"

#include "cycledecomp.h"

using namespace regina;

const unsigned int CycleDecompSearcher::faceEdges[4][3] = {
    {3,4,5}, {1,2,5}, {0,2,4}, {0,1,3}
};

const unsigned int CycleDecompSearcher::otherFace[4][3] = {
    {3,2,1}, {3,2,0}, {3,1,0}, {2,1,0}
};

const signed int CycleDecompSearcher::edgeParity[6][6] = {
    {-1,2,1,4,3,-1}, {2,-1,0,5,-1,3}, {1,0,-1,-1,5,4},
    {4,5,-1,-1,0,1}, {3,-1,5,0,-1,2}, {-1,3,4,1,2,-1}
};

const signed int CycleDecompSearcher::otherVert[4][6] = {
    {-1,-1,-1,3,2,1},
    {-1,3,2,-1,-1,0},
    {3,-1,1,-1,0,-1},
    {2,1,-1,0,-1,-1}
};

const char CycleDecompSearcher::dataTag_ = 'd';



CycleDecompSearcher::CycleDecompSearcher(const NFacePairing* pairing,
        const NFacePairing::IsoList* autos, bool orientableOnly,
        UseCycles use, void* useArgs) : use_(use), useArgs_(useArgs), 
        pairing_(pairing), orientable(orientableOnly) {

    nTets = pairing->size();
    nEdges = 2*nTets;
    nEnds = 2*nEdges;

    edgesLeft = 3*nEdges;
    
    tets = new Tetrahedron[nTets];
    edges = new Edge[nEdges];
    ends = new EdgeEnd[nEnds];
    
    for( unsigned int i=0; i < nTets; i++ ) {
        tets[i].index = i;
    }

    nextColour = 0;

    // Note that we need (nTets+1) cycles,
    // but we also don't use index 0 to denote
    // a cycle so we need space for (nTets+2) 
    // cycle descriptors.
    cycleLengths = new unsigned int[nTets+2];
    
    cycles = new signed int*[nTets+2];
    for( unsigned int i=0; i < nTets+2; i++ ) {
        cycleLengths[i] = 0;
        cycles[i] = new int[3*nEdges];
    }
    
    if ( autos ) {
        nAutos = autos->size();
        automorphisms = new signed int*[nAutos];
        for( unsigned int i=0; i < nAutos; i++ ) {
            automorphisms[i] = new signed int[nEdges];
        }
    } else {
        automorphisms = 0;
    }

    bool* orderAssigned = new bool[nTets * 4];
        /**< Have we placed a tetrahedron face or its partner in the
             order[] array yet? */

    NTetFace face, adj;
    unsigned int edgesDone = 0;
    std::fill(orderAssigned, orderAssigned + 4 * nTets, false);

    // Begin by searching for tetrahedra that are joined to themselves.
    // Note that each tetrahedra can be joined to itself at most once,
    // since we are guaranteed that the face pairing is connected with
    // order >= 3.

    for (face.setFirst(); ! face.isPastEnd(nTets, true); face++) {
        if (orderAssigned[face.simp * 4 + face.facet])
            continue;
        
        adj = (*pairing)[face];

        // Why did I ever put this here?
        //if (adj.simp != face.simp)
        //    continue;

        ends[2*edgesDone].tet = &tets[face.simp];
        ends[2*edgesDone+1].tet = &tets[adj.simp];
        ends[2*edgesDone].face = face.facet;
        ends[2*edgesDone+1].face = adj.facet;
        ends[2*edgesDone].edge = &edges[edgesDone];
        ends[2*edgesDone+1].edge = &edges[edgesDone];


        edges[edgesDone].ends[0]=&ends[2*edgesDone];
        edges[edgesDone].ends[1]=&ends[2*edgesDone+1];
        edges[edgesDone].index = edgesDone+1;
                
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

    if ( automorphisms ) {
        unsigned int autoCounter=0;
        for (NFacePairing::IsoList::const_iterator it = autos->begin();
                    it != autos->end(); it++) {
            for (unsigned int i=0; i< edgesDone;i++) {
                unsigned int startFace = edges[i].ends[0]->face;
                unsigned int startTet = edges[i].ends[0]->tet->index;
                unsigned int endFace = edges[i].ends[1]->face;
                unsigned int endTet = edges[i].ends[1]->tet->index;

                unsigned int newStartFace = (*it)->facePerm(startTet)[startFace];
                unsigned int newStartTet = (*it)->tetImage(startTet);
                unsigned int newEndFace = (*it)->facePerm(endTet)[endFace];
                unsigned int newEndTet = (*it)->tetImage(endTet);
                
                for(unsigned int j=0; j< edgesDone;j++) {
                    if (( newStartTet == edges[j].ends[0]->tet->index) &&
                        ( newStartFace == edges[j].ends[0]->face) &&
                        ( newEndTet == edges[j].ends[1]->tet->index) &&
                        ( newEndFace == edges[j].ends[1]->face)) {
                        assert((0 <= j) && (j < edgesDone));
                        automorphisms[autoCounter][i] = j;
                        break;
                    }
                    if (( newEndTet == edges[j].ends[0]->tet->index) &&
                        ( newEndFace == edges[j].ends[0]->face) &&
                        ( newStartTet == edges[j].ends[1]->tet->index) &&
                        ( newStartFace == edges[j].ends[1]->face)) {
                        assert((0 <= j) && (j < edgesDone));
                        automorphisms[autoCounter][i] = -j;
                        break;
                    }
                }
            }
            autoCounter+=1;
        }
        std::cout << "Original " << std::endl;

        cycleLengths[1] = 4;
        cycles[1][0] = 3;
        cycles[1][1] = 5;
        cycles[1][2] = 7;
        cycles[1][3] = 9;
        std::cout << "1: ";
        for(unsigned int i=0; i < cycleLengths[1]; i++) {
            std::cout << cycles[1][i];
            if (i < (cycleLengths[1]-1))
                std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "Automorphisms " << std::endl;
        
        for (autoCounter=0; autoCounter< nAutos; autoCounter++) {
            signed int newCycle[nEdges];
            std::cout << "1: ";
            for(unsigned int i=0; i < cycleLengths[1]; i++) {
                unsigned int index = abs(cycles[1][i]);
                if ( cycles[1][i] < 0 ) {
                    newCycle[i] = -1* automorphisms[autoCounter][index];
                } else {
                    newCycle[i] = automorphisms[autoCounter][index];
                }
                std::cout << newCycle[i];
                if (i < (cycleLengths[1]-1))
                    std::cout << ", ";
            }
            std::cout << std::endl;
            signed int internalEdges[6*nTet];
            for(unsigned int i=1; i < cycleLengths[1]; i++) {
                if edge[i-1] used negative:
                    edge[i-1].end[0].tet is this tet
                    edge[i-1].end[0].face is one face
                else:
                    edge[i-1].end[1].tet is this tet
                    edge[i-1].end[1].face is one face
                if edge[i] used negative
                    edge[i].end[1].face is other face
                else
                    edge[i].end[0].face is other face
                find common edge to one face and other face
                set internalEdges[this tet][common edge] = 1
            }
            print all the things


                edge[i-1] to edge[1];
                
        }
        cycleLengths[1] = 0;
        cycles[1][0] = 0;
        cycles[1][1] = 0;
        cycles[1][2] = 0;
        cycles[1][3] = 0;
    }
}

CycleDecompSearcher::~CycleDecompSearcher() {
    for( unsigned int i=0; i < nTets+1; i++ ) {
        delete[] cycles[i];
    }
    delete[] cycles;
    delete[] cycleLengths;    
    delete[] tets;
    delete[] edges;
    delete[] ends;
}

void CycleDecompSearcher::colourOnTetrahedra(unsigned int tet) {
    unsigned int edge;
    // Note empty for loop to find first unused internal edge.
    for (edge=0; tets[tet].internalEdges[edge] != 0; edge++) 
        assert(edge<=6);

    nextColour++;

    assert(cycleLengths[nextColour] == 0);

    
    tets[tet].internalEdges[edge] = nextColour;
    tets[tet].used++;
    unsigned int startFace = NEdge::edgeVertex[5-edge][0];
    EdgeEnd *start = tets[tet].externalEdgeEnd[startFace];
    unsigned int outFace = NEdge::edgeVertex[5-edge][1];
    EdgeEnd *outEdgeEnd = tets[tet].externalEdgeEnd[outFace];
    Edge *nextEdge = outEdgeEnd->edge;
    nextEdge->colour(nextColour);
    outEdgeEnd->map[edge] = nextEdge->used;
    edgesLeft--;
    EdgeEnd *nextEdgeEnd = nextEdge->otherEnd(outEdgeEnd);
    
    // An edge has two ends, ends[0] and ends[1]
    // A "forward" direction along an edge is from ends[0] to ends[1]
    // and "backwards" is the reverse. Forwards is denoted by a positive
    // edge number and backwards by a negative edge number.
    //
    // The below seems to be reversed, but the output makes sense this way.
    signed int dir = nextEdge->index;
    if (outEdgeEnd == (nextEdge->ends[1]) ) 
        dir = -1* nextEdge->index;
    
    cycles[nextColour][cycleLengths[nextColour]]=dir;
    cycleLengths[nextColour]++;

    nextPath(start, edge, nextEdgeEnd);
    
    cycles[nextColour][cycleLengths[nextColour]]=0;
    cycleLengths[nextColour]--;
    
    edgesLeft++;
    outEdgeEnd->map[edge] = 0;
    nextEdge->unColour(); 
    tets[tet].used--;
    tets[tet].internalEdges[edge] = 0;
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
        for(unsigned int i=0;i<3;i++) { 
            unsigned int ind = abs(cycles[nextColour][i]);
            // edges[0] has index 1
            Edge *e = &(edges[ind-1]);
            if (e->ends[0]->tet == e->ends[1]->tet) {
                foundLoop = true;
                break;
            }
        }
        if (foundLoop == false) 
            return false;
    }

    // Count the number of edges left. We need nTet+1 cycles,
    // and have just finished cycle number nextColour. We need 3 edges
    // for each of the (nTets+1 - nextColour) cycles after this one.
    
    if (edgesLeft < 3*(nTets+1 - nextColour)) {
        return false;
    }
    return true;
}

void CycleDecompSearcher::nextPath(EdgeEnd *start, unsigned int firstEdge,  EdgeEnd *now) {
    unsigned int nextInternal;
    EdgeEnd *nextEnd, *outEnd;
    Edge *nextEdge;
    Tetrahedron *nextTet = now->tet;

    unsigned int edgesLeftTest = edgesLeft;
    // Count the number of edges left. We need nTet+1 cycles,
    // and are working on cycle number nextColour. We need 3 edges
    // for each of the (nTets+1 - nextColour) cycles after this one.
    
    if (edgesLeft < 3*(nTets+1 - nextColour)) {
        return;
    }
    for (unsigned int i=0; i<3 ; i++) {
        nextInternal = faceEdges[now->face][i];
        if (nextTet->internalEdges[nextInternal] != 0) 
            continue;
        outEnd = nextTet->externalEdgeEnd[otherFace[now->face][i]];
        nextEdge = outEnd->edge;
        nextEnd = nextEdge->otherEnd(outEnd);
        
        // An edge has two ends, ends[0] and ends[1]
        // A "forward" direction along an edge is from ends[0] to ends[1]
        // and "backwards" is the reverse. Forwards is denoted by a positive
        // edge number and backwards by a negative edge number.
        signed dir = nextEdge->index;
        if (nextEnd == (nextEnd->edge->ends[0]) ) 
            dir = -1* nextEdge->index;
        // In the following loop, dir=0 means the edge has already been used in
        // the opposite direction.
        if (orientable) {
            for(unsigned j=0; j< cycleLengths[nextColour];j++) {
                if (cycles[nextColour][j] == -1*dir) {
                    dir=0;
                    break;  
                }
            }
            if (dir == 0) 
                continue;
        }
        
        assert(cycleLengths[nextColour] < 3*nEdges);
        cycles[nextColour][cycleLengths[nextColour]]=dir;
        cycleLengths[nextColour]++;
        
        nextTet->internalEdges[nextInternal] = nextColour;
        nextTet->used++;
        now->map[nextInternal] = now->edge->used;
        nextEdge->colour(nextColour);
        outEnd->map[nextInternal] = nextEdge->used;
        edgesLeft--;
       
        // Try to complete the cycle
        if (nextEnd == start) {
            start->map[firstEdge] = start->edge->used;
            if (checkColourOk()) {
                if ( (edgesLeft == 0) &&  checkComplete()) {
                    use_(this, useArgs_);
                } else {
                    // At most (nTets+1) cycles
                    if (nextColour <= nTets)
                        colourOnTetrahedra(findTetWithMostInternalEdgesUsed());
                }
            }
            start->map[firstEdge] = 0;
        }

        nextPath(start, firstEdge, nextEnd);
        

        cycleLengths[nextColour]--;
        cycles[nextColour][cycleLengths[nextColour]]=0;

        edgesLeft++;
        outEnd->map[nextInternal] = 0;
        now->map[nextInternal] = 0;
        nextEdge->unColour();
        nextTet->internalEdges[nextInternal] = 0;
        nextTet->used--;
    }
    assert(edgesLeft == edgesLeftTest);
}
 

bool CycleDecompSearcher::checkComplete() {
    for(unsigned int i=0; i< nEdges; i++) {
        if (edges[i].used != 3) {
            return false;
        }
    }
    return true;
}

unsigned int CycleDecompSearcher::findTetWithMostInternalEdgesUsed() {
    unsigned int mostUsed=0;
    // If nothing is used, start on tetrahedra 0.
    unsigned int tet=0;
    for(unsigned int i=0; i< nTets; i++) {
        if ((tets[i].used < 6) && (tets[i].used > mostUsed)) {
            mostUsed = tets[i].used;
            tet = i;
        }
    }
    assert(tets[tet].used < 6);
    return tet;
}

void CycleDecompSearcher::runSearch(long maxDepth) {
    // Preconditions:
    //     Only closed prime minimal P2-irreducible triangulations are needed.
    //     The given face pairing is closed with order >= 3.

    // Begin by testing for face pairings that can never lead to such a
    // triangulation.
    if (pairing_->hasTripleEdge() ||
            pairing_->hasBrokenDoubleEndedChain() ||
            pairing_->hasOneEndedChainWithDoubleHandle() ||
            pairing_->hasOneEndedChainWithStrayBigon() ||
            pairing_->hasWedgedDoubleEndedChain() ||
            pairing_->hasTripleOneEndedChain()) {
        use_(0, useArgs_);
        return;
    }


    unsigned tet = findTetWithMostInternalEdgesUsed();
    colourOnTetrahedra(tet);
    use_(0, useArgs_);
}

NTriangulation* CycleDecompSearcher::triangulate() const {
    NTriangulation* ans = new NTriangulation;
    NTetrahedron** simp = new NTetrahedron*[nTets];
    unsigned int t,a,b;
    signed int k;
    for (t = 0; t < nTets; ++t)
        simp[t] = ans->newSimplex();
    int perms[4];
    Edge *e;
    dumpData(std::cout);
    for (t = 0; t < nEdges; ++t) {
        e = &(edges[t]);
        // The edge e has 3 colours in it. These are denoted by k,
        // and can be 1, 2 or 3.
        for (k = 1; k <= 3; k++) {
            // The edge e represents 3 distinct edges.
            // The following loops work out, for each distinct edge,
            // which internal edge on either end is connected to
            // this edge.
            for (a = 0; (e->ends[0]->map[a] != k); a++) 
                assert(a<6);
            for (b = 0; (e->ends[1]->map[b] != k); b++) 
                assert(b<6);

            // Use otherVert to determine the "other" vertex on this face
            // of the tetrahedron, the vertex that is not on this edge.
            // Note that we don't know how edges a and b are glued together 
            // (in which direction that is) but we do know that the third
            // vertex on each face must be identified.
            int vertA = otherVert[e->ends[0]->face][a];
            assert(vertA >= 0);
            int vertB = otherVert[e->ends[1]->face][b];
            assert(vertB >= 0);

            perms[vertA] = vertB;
        }
        perms[e->ends[0]->face] = e->ends[1]->face;
        
        int thisTet = e->ends[0]->tet->index;
        int otherTet = e->ends[1]->tet->index;
        int thisFace = e->ends[0]->face;
        NPerm4 gluing(perms[0],perms[1],perms[2],perms[3]);
        simp[thisTet]->joinTo(thisFace, simp[otherTet], gluing);
    }
    return ans;
}

void CycleDecompSearcher::dumpData(std::ostream& out) const {
    //NCompactSearcher::dumpData(out);

    unsigned i;
    unsigned j;
    for (i = 1; i < nTets+2; i++) {
        out << i << ": ";
        for(j = 0; j < cycleLengths[i]; j++) {
            out << cycles[i][j];
            if (j != cycleLengths[i]-1) 
                out << ", ";
        }
        out << std::endl;
    }
    for (i = 0; i < nTets; i++) {
        for(j = 0; j < 6; j++) {
            out << tets[i].internalEdges[j] << " ";
        }
        out << " : ";
    }
    //for (i = 0; i < nEdges; i++) {
    //    out << "Edge " << i << ": [";
    //    for(j = 0; j < 6; j++) {
    //        unsigned int c = edges[i].ends[0]->map[j];
    //        assert(c < 3);
    //        out << std::endl << i << " " << j << " " << c << std::endl;
    //        if (c > 0) {
    //            out << c << "(" << edges[i].colours[c] << ")";
    //        } else {
    //            out << "0";
    //        }
    //        if (j != 5)
    //            out << ", ";
    //    }
    //    out << "] -> [";
    //    for(j = 0; j < 6; j++) {
    //        int c = edges[i].ends[1]->map[j];
    //        if (c > 0) {
    //            out << c << "(" << edges[i].colours[c] << ")";
    //        } else {
    //            out << "0";
    //        }
    //        if (j != 5)
    //            out << ", ";
    //    }
    //    out << "]" << std::endl;
    //}
}

inline void CycleDecompSearcher::Edge::colour(unsigned newColour) {
    assert(0 <= used && used < 3);
    colours[used] = newColour;
    used+=1;
}

inline void CycleDecompSearcher::Edge::unColour() {
    assert(0 < used && used <= 3);
    used-=1;
    colours[used] = 0;
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


