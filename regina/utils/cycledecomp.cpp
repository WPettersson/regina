
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
        UseCycles use, void* useArgs, bool minimal) : 
        use_(use), useArgs_(useArgs), 
        pairing_(pairing), orientable(orientableOnly),minimal_(minimal) {

    nTets = pairing->size();
    nEdges = 2*nTets;
    nEnds = 2*nEdges;

    edgesLeft = 3*nEdges;
    
    tets = new Tetrahedron[nTets];
    edges = new Edge[nEdges];
    ends = new EdgeEnd[nEnds];
  
    if (minimal) 
        nCycles = nTets+1;
    else
        nCycles = 3*nEdges;

    for( unsigned int i=0; i < nTets; i++ ) {
        tets[i].index = i;
    }
    

    nextColour = 0;

    // Note that we don't use index 0 to denote
    // a cycle so we need space for (nCycles+1) 
    // cycle descriptors.
    cycleLengths = new unsigned int[nCycles+1];
    
    cycles = new signed int*[nCycles+1];
    for( unsigned int i=0; i < nCycles+1; i++ ) {
        cycleLengths[i] = 0;
        cycles[i] = new signed int[3*nEdges];
    }
    

    bool* orderAssigned = new bool[nTets * 4];
        /* Have we placed a tetrahedron face or its partner in the
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


    if ( autos ) {
        nAutos = 0;
        automorphisms = new Automorphism*[autos->size()];
        for( NFacePairing::IsoList::const_iterator it = autos->begin();
                it != autos->end(); it++) {
            automorphisms[nAutos] = new Automorphism( (*it),edges,nEdges);
            nAutos++;
        }
    } else {
        automorphisms = 0;
    }
}

CycleDecompSearcher::~CycleDecompSearcher() {
    for( unsigned int i=0; i < nTets+1; i++ ) {
        delete[] cycles[i];
    }
    for( unsigned int i=0; i < nAutos; i++ ) {
        delete automorphisms[i];
    }
    delete[] automorphisms;
    delete[] cycles;
    delete[] cycleLengths;    
    delete[] tets;
    delete[] edges;
    delete[] ends;
}

void CycleDecompSearcher::colourLowestEdge() {
    unsigned int edge;
    // Note empty for loop to find first edge with space.
    for (edge=0; edges[edge].used == 3; edge++) 
        assert(edge<=nEdges);
    Edge *nextEdge = &edges[edge];
    unsigned int tet = nextEdge->ends[0]->tet->index;

    nextColour++;
    // Try each unused internal edge to start this cycle.
    for(unsigned int i=0; i < 3; i++) {
        unsigned int iEdge = faceEdges[nextEdge->ends[0]->face][i];
        if (tets[tet].internalEdges[iEdge] != 0)
            continue;

        assert(cycleLengths[nextColour] == 0);

        tets[tet].internalEdges[iEdge] = nextColour;
        tets[tet].used++;
        unsigned int startFace = NEdge::edgeVertex[5-iEdge][0];
        EdgeEnd *start = tets[tet].externalEdgeEnd[startFace];
        if ( start == nextEdge->ends[0] ) {
            startFace = NEdge::edgeVertex[5-iEdge][1];
            start = tets[tet].externalEdgeEnd[startFace];
        }


        nextEdge->colour(nextColour);
        nextEdge->ends[0]->map[iEdge] = nextEdge->used;
        edgesLeft--;
        EdgeEnd *nextEdgeEnd = nextEdge->ends[1];
        

        signed int dir = nextEdge->index;
        
        
        cycles[nextColour][cycleLengths[nextColour]]=dir;
        cycleLengths[nextColour]++;

        // Minimal triangulations won't have degree one edges.
        if (! minimal_ ) {
            // Try to complete the cycle
            if (nextEdgeEnd == start) {
                nextEdgeEnd->map[iEdge] = nextEdgeEnd->edge->used;
                if (checkColourOk() && isCanonical()) {
                    if ( (edgesLeft == 0) &&  checkComplete()) {
                        use_(this, useArgs_);
                    } else {
                        // At most nCycles cycles.  Remember nextColour starts at
                        // 1, and colourOnLowestEdge() increments it by one.
                        if (nextColour < nCycles) {
                            colourLowestEdge();
                        }
                    }
                }
                nextEdgeEnd->map[iEdge] = 0;
            }
        }
        // Try to find more paths.
        nextPath(start, iEdge, nextEdgeEnd);
        
        cycles[nextColour][cycleLengths[nextColour]]=0;
        cycleLengths[nextColour]--;
        
        edgesLeft++;
        nextEdge->ends[0]->map[iEdge] = 0;
        nextEdge->unColour(); 
        tets[tet].used--;
        tets[tet].internalEdges[iEdge] = 0;
    }
    nextColour--;
}

bool CycleDecompSearcher::checkColourOk() {
    // The tests below only apply if we are looking for minimal triangulations.
    if (! minimal_) 
        return true;
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

    return true;
}

void CycleDecompSearcher::nextPath(EdgeEnd *start, unsigned int firstEdge, 
        EdgeEnd *now) {
    unsigned int nextInternal;
    EdgeEnd *nextEnd, *outEnd;
    Edge *nextEdge;
    Tetrahedron *nextTet = now->tet;
    // Count the number of edges left. We are working on cycle number 
    // nextColour. We need 3 edges for each of the (nCycles - nextColour) 
    // cycles after this one. Note that this only applies if we are looking for
    // a minimal triangulation (else we won't actually know how many cycles we
    // are looking for.
   
    if (minimal_) {
        if (edgesLeft < 3*(nCycles - nextColour)) 
            return;
    } else { 
        if (edgesLeft <= 0)
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
            if (nTets >=3) {
                for(unsigned j=0; j< cycleLengths[nextColour];j++) {
                    if (cycles[nextColour][j] == dir) {
                        dir=0;
                        break;  
                    }
                }
                if (dir == 0)  {
                    continue;
                }
            }
        }
        
        //if ((nextColour == 1) && (! isCanonical(nextTet->index, nextInternal) )) {
        //    continue;
        //}
        
        
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
            if (checkColourOk() && isCanonical()) {
                if ( (edgesLeft == 0) && checkComplete()) {
                    use_(this, useArgs_);
                } else {
                    // At most nCycles cycles.  Remember nextColour starts at
                    // 1, and colourOnLowestEdge() increments it by one.
                    if (nextColour < nCycles)
                        colourLowestEdge();
                }
            } 
            start->map[firstEdge] = 0;
        }
        // Try to find more paths.
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
    if (minimal_ && (pairing_->hasTripleEdge() ||
            pairing_->hasBrokenDoubleEndedChain() ||
            pairing_->hasOneEndedChainWithDoubleHandle() ||
            pairing_->hasOneEndedChainWithStrayBigon() ||
            pairing_->hasWedgedDoubleEndedChain() ||
            pairing_->hasTripleOneEndedChain())) {
        use_(0, useArgs_);
        return;
    }


    colourLowestEdge();
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
    //dumpData(std::cout);
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

    unsigned int i;
    unsigned int j;
    for (i = 1; i < nCycles && cycleLengths[i] > 0; i++) {
        out << i << ": ";
        for(j = 0; j < cycleLengths[i]; j++) {
            out << cycles[i][j];
            if (j != cycleLengths[i]-1) 
                out << ", ";
        }
        out << std::endl;
    }
    out << "------------------" << std::endl;
    //for (i = 0; i < nTets; i++) {
    //    for(j = 0; j < 6; j++) {
    //        out << tets[i].internalEdges[j] << " ";
    //    }
    //    out << ": ";
    //}
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

bool CycleDecompSearcher::isCanonical() {
    unsigned int autoNo;
    bool thisAuto;
    unsigned int *cycleList[nextColour+1];
    unsigned int cycleListLengths[nextColour+1];
    unsigned int offset[nCycles];
    for(unsigned int i=1; i<=nextColour; i++) {
        cycleList[i] = new unsigned int[3*nEdges];//cycleLengths[i]];
        //std::cout << "cycleList["<<i<<"] has size " << cycleLengths[i];
        //std::cout << " At " << cycleList[i] << std::endl;
    }
    //std::cout << "Checking" << std::endl;
    //for(unsigned int k=1; k<=nextColour; k++) {
    //    for(unsigned int l=0; l < cycleLengths[k]; l++) {
    //        std::cout << cycles[k][l] << " ";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "---" << std::endl;

    //bool debug=true;
    
    for(autoNo=0; autoNo < nAutos; autoNo++) {
        // As long as thisAuto is true, we keep looking at this automorphism.
        // If it becomes false, it means that this automorphism results in a
        // less canonical representation.
        thisAuto=true;
        // Generate new cycle lists
        for(unsigned int i=1; i<=nextColour; i++) {
            unsigned int min = nEdges;
            offset[i]=0;
            // offset[i] is going to hold the current "starting" place for the
            // cycle.

            unsigned int nextEdge;
            // nextEdge will store the "second" edge in the cycle, if offset[i]
            // is the first edge.  Note that the "second" edge may be the one
            // before the first edge, if the first edge is negative (since
            // swapping signs also reverses direction).
            signed int tempNewEdge = 2*(*automorphisms[autoNo])[cycles[i][0]];
            if (tempNewEdge < 0) {
                tempNewEdge = (-tempNewEdge) + 1;
            }
            assert(tempNewEdge>0);
            unsigned int newEdge = tempNewEdge;
            cycleList[i][0] = newEdge;

            bool checkNextPair=false;
            bool setNextEdge=false;
            for(unsigned int j=1; j < cycleLengths[i]; j++) {
                tempNewEdge = 2*(*automorphisms[autoNo])[cycles[i][j]];
                if (tempNewEdge < 0) {
                    tempNewEdge = (-tempNewEdge) + 1;
                }
                newEdge = tempNewEdge;
                // Initialise the "nextEdge" variable.
                if ( j==1 ) {
                    nextEdge = newEdge;
                }

                // If checkNextPair is true, that means the last edge we
                // checked was equal-smallest.  We check this new edge against
                // the "next" edge from the current smallest, and if we have a
                // smaller edge, update the offset.
                if (checkNextPair) {
                    if (newEdge < nextEdge) {
                        offset[i]=j-1;
                        nextEdge = newEdge;
                    }
                    checkNextPair=false;
                }
                // If setNextEdge is true, that means we reset the offset in
                // the last iteration through, so we need to update what the
                // "next" edge after the offset is.
                if (setNextEdge) {
                    nextEdge = newEdge;
                    setNextEdge = false;
                }
                // New lowest edge used in this cycle.
                if ( newEdge < min ) { 
                    min = newEdge;
                    offset[i] = j;
                    if (min%2 ==1) {
                        // Negative
                            
                        // If i%2 == 0, then 1 - 2*(i%2) == 1 
                        // If i%2 == 1, then 1 - 2*(i%2) == -1
                        //
                        // So the line below adds one (if i%2 ==0) or
                        // subtracts one (if i%2 == 1), which flips the
                        // sign.
                        nextEdge = cycleList[i][j-1] + (1- 2*(cycleList[i][j-1]%2));
                    } else {
                        // Positive
                        setNextEdge = true;
                    }
                } else {
                    // Check to see if we have a tie.
                    if ( newEdge == min ) { 
                        // Have to check whether the next edge is smaller or
                        // not.  Don't forget that "next" might be previous if
                        // the lowest edge is a -ve.
                       
                        if ( newEdge %2 == 1) { // lowest edge is a -ve
                            // If the current lowest is negative, we will be
                            // flipping the signs of all edges, so don't forget
                            // that in this comparison.

                            // If i%2 == 0, then 1 - 2*(i%2) == 1 
                            // If i%2 == 1, then 1 - 2*(i%2) == -1
                            //
                            // So the line below adds one (if i%2 ==0) or
                            // subtracts one (if i%2 == 1), which flips the
                            // sign.
                            unsigned int newNext = cycleList[i][j-1] + (1 - 2*(cycleList[i][j-1]%2));

                            if (newNext < nextEdge) {
                                    offset[i]=j;
                                    nextEdge = newNext;
                                    // min stays the same.
                            }
                        } else { 
                                // min and current edges are positive. 

                            if ( j == cycleLengths[i]-1) {
                                // "next" edge is first
                                unsigned int newNext = cycleList[i][0];
                                if  ( newNext < nextEdge) {
                                    offset[i]=j;
                                    // no need to change nextEdge
                                }

                            } else {
                                // We need to check the next edge against the
                                // "next edge" from the current minimum, so
                                // make this note.
                                checkNextPair = true;
                            }
                        }
                    } else {
                        // This check here is for when the current "lowest edge in
                        // the cycle" is positive, but the edge we're checking now
                        // is negative but has same absolute value.
                        if ( (newEdge == (min+1)) && (min%2 == 0)) {
                        
                            // If i%2 == 0, then 1 - 2*(i%2) == 1 
                            // If i%2 == 1, then 1 - 2*(i%2) == -1
                            //
                            // So the line below adds one (if i%2 ==0) or
                            // subtracts one (if i%2 == 1), which flips the
                            // sign.
                            unsigned int newNext = cycleList[i][j-1] + (1 - 2*(cycleList[i][j-1]%2));
                            if ( newNext < nextEdge) {
                                min = newEdge;
                                nextEdge = newNext;
                                offset[i]=j;
                            }
                        }
                    }
                }
                cycleList[i][j] = newEdge;
            }
            cycleListLengths[i] = cycleLengths[i];
        }
        //if (debug) {
        //    std::cout << "Before sorting" << std::endl;
        //    for(unsigned int k=1; k<=nextColour; k++) {
        //        for(unsigned int l=0; l < cycleListLengths[k]; l++) {
        //            signed int val = cycleList[k][l];
        //            if ( l == offset[k] ) {
        //                std::cout << "(";
        //            }
        //            if (val % 2 == 1) {
        //                val = - ( val-1)/2;
        //            } else {
        //                val = val/2;
        //            }
        //            std::cout << val;
        //            if ( l == offset[k] ) {
        //                std::cout << ")";
        //            }
        //            std::cout << " ";
        //        }
        //        std::cout << std::endl;
        //    }
        //    std::cout << "---" << std::endl;
        //}

        // Sort cycleList based on values of cycleList[i][offset[i]] 
        
        // Reverse bubble sort order so position[0] is definitely lowest
        // after one run through.  We can then check to see if this is more or
        // less canonical after 1 run through.
        for(unsigned int i=1; thisAuto && i <= nextColour ; i++ ) {
            // Do some sorting!
            for(unsigned int j=nextColour; j > i; j--) {
                if ( compareCycles(cycleList[j],        cycleList[j-1],
                                   cycleListLengths[j], cycleListLengths[j-1],
                                   offset[j],           offset[j-1]) == 0) {
                    unsigned int *temp;
                    temp = cycleList[j];
                    cycleList[j] = cycleList[j-1];
                    cycleList[j-1] = temp;
                    unsigned int temp2 = offset[j];
                    offset[j] = offset[j-1];
                    offset[j-1] = temp2;
                    temp2 = cycleListLengths[j];
                    cycleListLengths[j] = cycleListLengths[j-1];
                    cycleListLengths[j-1] = temp2;
                }

            }
            //if (debug) {
            //    std::cout << "Sorted " << i << " time(s)" << std::endl;
            //    for(unsigned int k=1; k<=nextColour; k++) {
            //        for(unsigned int l=0; l < cycleListLengths[k]; l++) {
            //            signed int val = cycleList[k][l];
            //            if ( l == offset[k] ) {
            //                std::cout << "(";
            //            }
            //            if (val % 2 == 1) {
            //                val = - ( val-1)/2;
            //            } else {
            //                val = val/2;
            //            }
            //            std::cout << val;
            //            if ( l == offset[k] ) {
            //                std::cout << ")";
            //            }
            //            std::cout << " ";
            //        }
            //        std::cout << std::endl;
            //    }
            //    std::cout << "---" << std::endl;
            //}
            unsigned int counterA=offset[i];
            unsigned int maxLength = cycleLengths[i] < cycleListLengths[i] ? cycleLengths[i] : cycleListLengths[i];
            for(unsigned int j=0; j < maxLength; j++) {
                signed int t = 2*cycles[i][j];
                if (t < 0) {
                    t = (-t)+1;
                }
                unsigned int e = t;
                unsigned newE = cycleList[i][counterA];
                if (cycleList[i][offset[i]] %2 == 1) {
                    // If i%2 == 0, then 1 - 2*(i%2) == 1 
                    // If i%2 == 1, then 1 - 2*(i%2) == -1
                    //
                    // So the line below adds one (if i%2 ==0) or
                    // subtracts one (if i%2 == 1), which flips the
                    // sign.
                    newE += (1 - 2*(newE%2));
                }

                // the automorphism gives a more canonical representation
                if (e > newE) {
                    //std::cout << " bigger on " << j << std::endl;
                    //for(unsigned int k=1; k<=nextColour; k++) {
                    //    for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                    //        signed int val = cycleList[k][l];
                    //        if ( l == offset[k] ) {
                    //            std::cout << "(";
                    //        }
                    //        if (val % 2 == 1) {
                    //            val = - ( val-1)/2;
                    //        } else {
                    //            val = val/2;
                    //        }
                    //        std::cout << val;
                    //        if ( l == offset[k] ) {
                    //            std::cout << ")";
                    //        }
                    //        std::cout << " ";
                    //    }
                    //    std::cout << std::endl;
                    //}
                    //std::cout << "---" << std::endl;
                    for(unsigned int k=1; k<=nextColour; k++)
                        delete[] cycleList[k];
                    return false;
                }
                // this automorphism is definitely worse, so go to next.
                if (e < newE) {
                    thisAuto=false;
                    break;
                }

                // Increment counter in appropriate direction.
                if (cycleList[i][offset[i]] %2 == 1) { // Negative
                    if ( counterA == 0) {
                        counterA = cycleListLengths[i]-1;
                    } else {
                        counterA--;
                    }
                } else {
                    counterA++;
                    if ( counterA == cycleListLengths[i] ) {
                        counterA = 0;
                    }
                }
            }
            // Check cycle lengths, shorter cycles are more canonical

            if ( cycleLengths[i] < cycleListLengths[i] ) {
                thisAuto=false;
            }
            if ( thisAuto && ( cycleListLengths[i] < cycleLengths[i] )) {
                //std::cout << "Beaten by length at cycle " << i << std::endl;
                //for(unsigned int k=1; k<=nextColour; k++) {
                //    for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                //        signed int val = cycleList[k][l];
                //        if ( l == offset[k] ) {
                //            std::cout << "(";
                //        }
                //        if (val % 2 == 1) {
                //            val = - (val-1)/2;
                //        } else {
                //            val = val/2;
                //        }
                //        std::cout << val;
                //        if ( l == offset[k] ) {
                //            std::cout << ")";
                //        }
                //        std::cout << " ";
                //    }
                //    std::cout << std::endl;
                //}
                //std::cout << "---" << std::endl;
                for(unsigned int k=1; k<=nextColour; k++)
                    delete[] cycleList[k];
                return false;
            }
        }
    }
    for(unsigned int i=1; i<=nextColour; i++)
        delete[] cycleList[i];
    return true;
}
// Check to see which of two cycles is more canonical.  Return values are
//  0 if A < B
//  1 if B > A
//  2 if A == B
unsigned int CycleDecompSearcher::compareCycles(unsigned int *cycleListA, 
        unsigned int *cycleListB, unsigned int lengthA, unsigned int lengthB,
        unsigned int offsetA, unsigned int offsetB) {
    unsigned int maxLength = lengthA < lengthB ? lengthA : lengthB;
    unsigned int counterA=offsetA;
    unsigned int counterB=offsetB;
    for(unsigned int i=0 ; i < maxLength; i++) {
        // Find next edges
        unsigned int edgeA = cycleListA[counterA];
        unsigned int edgeB = cycleListB[counterB];

        if (cycleListA[offsetA] %2 == 1) {
            // If i%2 == 0, then 1 - 2*(i%2) == 1 
            // If i%2 == 1, then 1 - 2*(i%2) == -1
            //
            // So the line below adds one (if i%2 ==0) or
            // subtracts one (if i%2 == 1), which flips the
            // sign.
            edgeA += (1 - 2*(edgeA%2));
        }

        if (cycleListB[offsetB] %2 == 1) {
            // If i%2 == 0, then 1 - 2*(i%2) == 1 
            // If i%2 == 1, then 1 - 2*(i%2) == -1
            //
            // So the line below adds one (if i%2 ==0) or
            // subtracts one (if i%2 == 1), which flips the
            // sign.
            edgeB += (1 - 2*(edgeB%2));
        }
        // Compare edges
        if (edgeA < edgeB)
            return 0;
        if (edgeB < edgeA)
            return 1;

        // Move counters along.  Note that we have to check whether the first
        // edge is positive or negative.
        if (cycleListA[offsetA] %2 == 1) { // Negative
            if ( counterA == 0) {
                counterA = lengthA-1;
            } else {
                counterA--;
            }
        } else {
            counterA++;
            if ( counterA == lengthA ) {
                counterA = 0;
            }
        }
        if (cycleListB[offsetB] %2 == 1) { // Negative
            if ( counterB == 0) {
                counterB = lengthB-1;
            } else {
                counterB--;
            }
        } else {
            counterB++;
            if ( counterB == lengthB ) {
                counterB = 0;
            }
        }
    }
    // So far all edges are the same.  shorter cycles are "more" canonical.
    if (lengthA < lengthB)
        return 0;
    if (lengthB < lengthA)
        return 1;
    // All edges the same, lengths the same.  Same cycles.
    return 2; 
}


inline void CycleDecompSearcher::Edge::colour(unsigned newColour) {
    colours[used++] = newColour;
}

inline void CycleDecompSearcher::Edge::unColour() {
    colours[--used] = 0;
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

CycleDecompSearcher::Automorphism::Automorphism(const NIsomorphism * iso,
        const Edge *edges, const unsigned int _nEdges) : nEdges(_nEdges)  {
    edgeMap = new signed int[2*(nEdges+1)];
    realEdgeMap = edgeMap+nEdges;
    for (unsigned int i=0; i < nEdges;i++) {
        unsigned int startFace = edges[i].ends[0]->face;
        unsigned int startTet = edges[i].ends[0]->tet->index;
        unsigned int endFace = edges[i].ends[1]->face;
        unsigned int endTet = edges[i].ends[1]->tet->index;

        unsigned int newStartFace = iso->facePerm(startTet)[startFace];
        unsigned int newStartTet = iso->tetImage(startTet);
        unsigned int newEndFace = iso->facePerm(endTet)[endFace];
        unsigned int newEndTet = iso->tetImage(endTet);
        signed int iIndex = edges[i].index; 
        for(unsigned int j=0; j< nEdges;j++) {
            if (( newStartTet == edges[j].ends[0]->tet->index) &&
                ( newStartFace == edges[j].ends[0]->face) &&
                ( newEndTet == edges[j].ends[1]->tet->index) &&
                ( newEndFace == edges[j].ends[1]->face)) {
                // Edge parity stays the same.
                signed int jIndex = edges[j].index;
                realEdgeMap[iIndex] = jIndex;
                realEdgeMap[std::ptrdiff_t(-iIndex)] = -jIndex;
                break;
            }
            if (( newEndTet == edges[j].ends[0]->tet->index) &&
                ( newEndFace == edges[j].ends[0]->face) &&
                ( newStartTet == edges[j].ends[1]->tet->index) &&
                ( newStartFace == edges[j].ends[1]->face)) {
                signed int jIndex = edges[j].index;
                realEdgeMap[iIndex] = -jIndex;
                realEdgeMap[std::ptrdiff_t(-iIndex)] = jIndex;
                break;
            }
        }
    }
    //std::cout << " New Auto " << std::endl;
    //for(unsigned int i=1; i<= nEdges; i++) {
    //    std::cout << "auto["<<i<<"] = " << realEdgeMap[i] << std::endl;
    //}
}

CycleDecompSearcher::Automorphism::~Automorphism() {
    delete[] edgeMap;
//    for (unsigned int i=0; i < nTets;i++) 
//      delete[] newInts[i];
//    delete[] newInts;
}

signed int inline CycleDecompSearcher::Automorphism::operator [] (const signed int in) {
  //return edgeMap[in+nEdges];
  return realEdgeMap[(std::ptrdiff_t)(in)];
}

//void CycleDecompSearcher::Automorphism::tetAndInt(
//        unsigned int *newTet, unsigned int *newInternal, 
//        unsigned int oldTet, unsigned int oldInternal) {
//  *newTet = newTets[oldTet];
//  *newInternal = newInts[oldTet][oldInternal];
//  return;
//}

/* vim: set ts=4 sw=4 expandtab: */
