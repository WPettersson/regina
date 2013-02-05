
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
  
    // Any minimal triangulation that is not a triangulation of one of RP^3,
    // S^3 or L(3,1) has precisely one vertex.  
    // If we force the number of edges in the triangulation to be n+1 we get
    // \chi = 0 (Euler's formula).
    if (nTets >= 3 && minimal) 
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
    parityArray = new unsigned int[nCycles+1];
    parityArrayCount = new unsigned int[nCycles+1];
    for( unsigned int i=0; i < nCycles+1; i++ ) {
        cycleLengths[i] = 0;
        parityArray[i] = 0;
        parityArrayCount[i] = 0;
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

    // Generate the list of face pairing automorphisms if necessary.
    if (autos == 0) {
        autos = new NFacePairing::IsoList();
        pairing->findAutomorphisms(const_cast<NFacePairing::IsoList&>(*autos));
    }

    if (autos) {
        nAutos = 0;
        automorphisms = new Automorphism*[autos->size()];
        for( NFacePairing::IsoList::const_iterator it = autos->begin();
                it != autos->end(); it++) {
            automorphisms[nAutos] = new Automorphism( (*it),edges,nEdges, nCycles);
            nAutos++;
        }
    } else {
        automorphisms = 0;
    }
}

CycleDecompSearcher::~CycleDecompSearcher() {
    for( unsigned int i=0; i < nCycles+1; i++ ) {
        delete[] cycles[i];
    }
    if (automorphisms != 0) {
        for( unsigned int i=0; i < nAutos; i++ ) {
            delete automorphisms[i];
        }
        delete[] automorphisms;
    }
    delete[] cycles;
    delete[] parityArray;
    delete[] parityArrayCount;
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
    Tetrahedron * tet = nextEdge->ends[0]->tet;
    tet->used++;

    nextColour++;
    nextEdge->colour(nextColour);
    // Try each unused internal edge to start this cycle.
    for(unsigned int i=0; i < 3; i++) {
        unsigned int iEdge = faceEdges[nextEdge->ends[0]->face][i];
        if (tet->internalEdges[iEdge] != 0)
            continue;

        assert(cycleLengths[nextColour] == 0);

        tet->internalEdges[iEdge] = nextColour;
                
        // Note that edgeVertex[a] gives the two vertices on edge a
        // We want the two corresponding faces, so do edgeVertex[5-a]
        unsigned int startFace = NEdge::edgeVertex[5-iEdge][0];
        EdgeEnd *start = tet->externalEdgeEnd[startFace];
        
        
        if ( start == nextEdge->ends[0] ) {
            startFace = NEdge::edgeVertex[5-iEdge][1];
            start = tet->externalEdgeEnd[startFace];
        }


        nextEdge->ends[0]->map[iEdge] = nextEdge->used;
        edgesLeft--;
        EdgeEnd *nextEdgeEnd;
        nextEdgeEnd = nextEdge->ends[1];
        

        signed int dir = nextEdge->index;
        
        
        cycles[nextColour][cycleLengths[nextColour]]=dir;
        cycleLengths[nextColour]++;
        


        // Minimal triangulations with >= 3 tets won't have degree one edges.
        if ((nTets < 3) || (! minimal_) ) {
            // Try to complete the cycle
            if (nextEdgeEnd == start) {
                nextEdgeEnd->map[iEdge] = nextEdgeEnd->edge->used;
                bool goodGluing = finishTet(tet);
                if (goodGluing && checkColourOk() && isCanonical()) {
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
                unFinishTet(tet);
                nextEdgeEnd->map[iEdge] = 0;
            }
        }

        bool goodGluing = finishTet(tet);
        if (goodGluing) {
            // Try to find more paths.
            nextPath(start, iEdge, nextEdgeEnd);
        }

        unFinishTet(tet);

        cycles[nextColour][cycleLengths[nextColour]]=0;
        cycleLengths[nextColour]--;
        
        edgesLeft++;
        nextEdge->ends[0]->map[iEdge] = 0;
        tet->internalEdges[iEdge] = 0;
    }
    tet->used--;
    nextEdge->unColour(); 
    nextColour--;
}

bool CycleDecompSearcher::checkColourOk() {
    // The tests below only apply if we are looking for minimal triangulations.
    if (! minimal_) 
        return true;
    // Tests only work if we have 3 or more tetrahedron.
    if (nTets < 3)
        return true;
    // No short cycles
    if ( cycleLengths[nextColour] < 4) 
        return false;
    return true;
}

void CycleDecompSearcher::nextPath(EdgeEnd *start, unsigned int firstEdge, 
        EdgeEnd *now) {
    unsigned int nextInternal;
    EdgeEnd *nextEnd, *outEnd;
    Edge *nextEdge;
    Tetrahedron *nextTet = now->tet;
    // Count the number of edges left. We are working on cycle number 
    // nextColour. We need 4 edges for each of the (nCycles - nextColour) 
    // cycles after this one. Note that this only applies if we are looking for
    // a minimal triangulation (else we won't actually know how many cycles we
    // are looking for.
   
    if ((nTets >= 3) && (minimal_)) {
        if (edgesLeft < 4*(nCycles - nextColour)) 
            return;
    } else { 
        // We don't actually know how many cycles we might have, so we just
        // have to keep using edges til they run out.
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
                    if (cycles[nextColour][j] == -dir) {
                        dir=0;
                        break;  
                    }
                }
                if (dir == 0)  {
                    continue;
                }
            }
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
            bool goodGluing = finishTet(nextTet);
            if (goodGluing && checkColourOk() && isCanonical()) {
                if ( (edgesLeft == 0) && checkComplete()) {
                    use_(this, useArgs_);
                } else {
                    // At most nCycles cycles.  Remember nextColour starts at
                    // 1, and colourOnLowestEdge() increments it by one.
                    if (nextColour < nCycles)
                        colourLowestEdge();
                }
            } 
            unFinishTet(nextTet);
            start->map[firstEdge] = 0;
        }
        // Try to find more paths.
        bool goodGluing = finishTet(nextTet);
        if (goodGluing) {
            nextPath(start, firstEdge, nextEnd);
        } 
        unFinishTet(nextTet);

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


// Return true if there is nothing wrong.
bool CycleDecompSearcher::finishTet(Tetrahedron *tet) {
    if (nTets <= 3) {
        return true;
    }
    if ( tet->used < 6 ) {
        return true;
    }
    bool goodGluing = true;
    //std::cout << "Finished tet " << tet->index << std::endl;
    //std::cout << "Internal: " << tet->internalEdges[0] << 
    //            " " << tet->internalEdges[1]  <<
    //            " " << tet->internalEdges[2]  <<
    //            " " << tet->internalEdges[3]  <<
    //            " " << tet->internalEdges[4]  <<
    //            " " << tet->internalEdges[5]  <<
    //            std::endl;
    // tet has just been fully completed.  Iterate over all
    // internal edges and check whether edges are above or below.
    for (unsigned int i=0; i< 6;i++) {
        // Find the ends for this actual internal edge.
        
        // Note that edgeVertex[a] gives the two vertices on edge a
        // We want the two corresponding faces, so do edgeVertex[5-a]
        unsigned int endA = NEdge::edgeVertex[5-i][0];
        unsigned int endB = NEdge::edgeVertex[5-i][1];
        // Convert the objects into a numerical value
        EdgeEnd *aa = tet->externalEdgeEnd[endA];
        EdgeEnd *bb = tet->externalEdgeEnd[endB];
        // Get the colour
        unsigned int col = tet->internalEdges[i];
        assert(col>0);
        //std::cout << "i=" << i << std::endl;
        //std::cout << "endA=" << endA << std::endl;
        //std::cout << "endB=" << endB << std::endl;
        unsigned int a=0;
        unsigned int b=0;
        unsigned int c=0;
        unsigned int d=0;
        //if (col == 2) {
        //    std::cout << "OnTet "<< tet->index << " with internal " << i << std::endl;
        //    std::cout << "Faces are " << endA << " and " << endB << std::endl;
        //    std::cout << "Count on col " << col << " is " << parityArrayCount[col]  << std::endl;
        //}
        //unsigned int result = 0;
        for (unsigned int j=0; j < 3; j++) {
            unsigned int e = faceEdges[endA][j];
            // Don't map the current cycle anywhere, only the two edges
            // around it.
            if (e == i) {
                continue;
            }
            // Find the resultant edge.
            signed int eb = edgeParity[i][e];
            //std::cout << "i " << i << " e:eb " << e <<":"<<eb<< std::endl;
            assert(eb>=0);


            // If we have started a new cycle on this tet, then we may not have
            // set up the map[] array on the very last edge.  If not, since
            // it's the last edge to be added, we know it must be the "third"
            // one used.
            unsigned int aAdjust = aa->map[e];
            if (aAdjust == 0) {
                aAdjust = 3;
            }
            unsigned int bAdjust = bb->map[eb];
            if (bAdjust == 0) {
                bAdjust = 3;
            }

            unsigned int valA = 3*(aa->edge->index) + aAdjust; 
            unsigned int valB = 3*(bb->edge->index) + bAdjust; 

            if (a==0) {
                a=valA;
                b=valB;
                //if (col == 2) {
                //    std::cout << "e: " << e << " eb: " << eb <<  std::endl;
                //    std::cout << "a: " << a << " b: " << b <<  std::endl;
                //}
            } else {
                c=valA;
                d=valB;
                //if (col == 2) {
                //    std::cout << "e: " << e << " eb: " << eb <<  std::endl;
                //    std::cout << "c: " << c << " d: " << d <<  std::endl;
                //}
            }
            
            //unsigned int res = ufJoin(col,valA,valB);
            //if (result == 0) {
            //    result = res;
            //} else {
            //    if (result == res) {
            //        std::cout << "Bad triangulation, res = " << res << std::endl;
            //        std::cout << "Completed tet " << tet->index << std::endl;
            //        std::cout << "Completed colour " << col << std::endl;
            //        std::cout << "Edges joined = " << aa->edge->index << ", " <<
            //            bb->edge->index << std::endl;
            //        std::cout << "A,B = " << valA << " , " << valB << std::endl;
            //        std::cout << "Next colour " << nextColour << std::endl;
            //        std::cout << "Arr " << std::endl;
            //        for( unsigned int j=0; j< 3*(nEdges+1); j++) {
            //            std::cout << parityArrays[col][j] << " ";
            //        }
            //        std::cout << std::endl;
            //        dumpData(std::cout);
            //        goodGluing = false;
            //    }
            //}
        }
        assert(a!=0);
        assert(b!=0);
        assert(c!=0);
        assert(d!=0);
        
        parityArrayCount[col]+=1;

        // a,b,c,d have been assigned such that (a,c) <-> (b,d)
        // Check whether parity has changed.
        if ( ((a<c) && (d<b)) || ( (c<a) && (b<d))) {
            parityArray[col]+=1;
            //if (col == 2) {
            //    std::cout << "Flip on col " << col << " with " << a << " " << b << " " << c << " " << d << std::endl;
            //    std::cout << "Flips on col " << col << " is " << parityArray[col]  << std::endl;
            //}
        }
        //std::cout << "Col: " << col << " Len: " << cycleLengths[col] << " Count: " << parityArrayCount[col] << std::endl;
        //if (( col < nextColour ) && (parityArrayCount[col] == cycleLengths[col]) ) {
        //    std::cout << "Finished parity for colour " << col << std::endl;
        //    std::cout << "Parity is " << parityArray[col] << std::endl;
        //}

        // Now check to see if this results in a 2-sided projective
        // plane.  In other words, see if the "positive" and "negative" edges
        // swap an odd number of times. Note that we can only do this
        // if we have completed the current cycle, and also found all
        // possible "flips".
        if (( col < nextColour ) && (parityArrayCount[col] == cycleLengths[col]) && ( parityArray[col]%2 == 1)) {
            //std::cout << "Bad triangulation" << std::endl;
            //std::cout << "Completed tet " << tet->index << std::endl;
            //std::cout << "Completed colour " << col << std::endl;
            //std::cout << "Next colour " << nextColour << std::endl;
            //dumpData(std::cout);
            //std::cout << "Done" << std::endl;
            goodGluing=false;
        }
    }
    return goodGluing;
}

void CycleDecompSearcher::unFinishTet(Tetrahedron *tet) {
    if (nTets <= 3) {
        return;
    }
    if ( tet->used < 6) {
        return;
    }
    // tet has just been fully completed.  Iterate over all
    // internal edges and remove parity checks.
    for (unsigned int i=0; i< 6;i++) {
        // Find the ends for this actual internal edge.

        // Note that edgeVertex[a] gives the two vertices on edge a
        // We want the two corresponding faces, so do edgeVertex[5-a]
        unsigned int endA = NEdge::edgeVertex[5-i][0];
        unsigned int endB = NEdge::edgeVertex[5-i][1];
        // Get the colour
        unsigned int col = tet->internalEdges[i];
        assert(col>0);
        EdgeEnd *aa = tet->externalEdgeEnd[endA];
        EdgeEnd *bb = tet->externalEdgeEnd[endB];
        unsigned int a=0;
        unsigned int b=0;
        unsigned int c=0;
        unsigned int d=0;
        for (unsigned int j=0; j < 3; j++) {
            unsigned int e = faceEdges[endA][j];
            // Don't map the current cycle anywhere, only the two edges
            // around it.
            if (e == i) {
                continue;
            }
            // Find the resultant edge.
            signed int eb = edgeParity[i][e];
            assert(eb>=0);
            

            // If we have started a new cycle on this tet, then we may not have
            // set up the map[] array on the very last edge.  If not, since
            // it's the last edge to be added, we know it must be the "third"
            // one used.
            unsigned int aAdjust = aa->map[e];
            if (aAdjust == 0) {
                aAdjust = 3;
            }
            unsigned int bAdjust = bb->map[eb];
            if (bAdjust == 0) {
                bAdjust = 3;
            }

            // Convert the objects into a numerical value
            unsigned int valA = 3*(aa->edge->index) + aAdjust; 
            unsigned int valB = 3*(bb->edge->index) + bAdjust; 
            //ufUnJoin(col,valA,valB);
            
            if (a==0) {
                a=valA;
                b=valB;
            } else {
                c=valA;
                d=valB;
            }
        }
        
        parityArrayCount[col]-=1;
        //if (col == 2) {
        //    std::cout << "Undo on col " << col << " with " << a << " " << b << " " << c << " " << d << std::endl;
        //}
        // a,b,c,d have been assigned such that (a,c) <-> (b,d)
        // Check whether parity has changed, and undo the changes
        if ( ((a<c) && (b>d)) || ( (a>c) && (b<d))) {
            parityArray[col]-=1;
            //if (col == 2) {
            //    std::cout << "unFlip on col " << col << " with " << a << " " << b << " " << c << " " << d << std::endl;
            //    std::cout << "Flips on col " << col << " is " << parityArray[col]  << std::endl;
            //}
        }
    }
}
 
unsigned int CycleDecompSearcher::ufJoin(unsigned int col, unsigned int A, unsigned int B) {
    // Ensure A<B
    if ( B < A ) {
        return ufJoin(col, B, A);
    }
    // if A==B (which is possible) we don't need to update the array.  Still
    // find the "smallest" entry though.
    unsigned int *arr = parityArray;
    assert(1==0);
    if ( A < B ) {
        arr[B] = A;
    }
    unsigned int temp = A;
    while ( arr[temp] > 0) {
        temp = arr[temp];
    }
    assert(temp>0);
    return temp;
}


void CycleDecompSearcher::ufUnJoin(unsigned int col, unsigned int A, unsigned int B) {
    // Ensure A<B
    if ( B < A ) {
        ufUnJoin(col, B, A);
    } else {
        unsigned int *arr = parityArray;
        assert(1==0);
        arr[B] = 0;
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

void CycleDecompSearcher::runSearch(long maxDepth, 
        std::vector<ThreeCycle*> *threeCycleEdges, unsigned int threeCyclesDone) {
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


    if ( threeCycleEdges == 0 ) {
        // Find all appropriate places where cycles of length 3 can be placed
        threeCycleEdges = new std::vector<ThreeCycle*>;
        for(unsigned int i=0; i < nEdges; i++) {
            Edge *e = &edges[i]; // edges start at 1 and go to nEdges
            unsigned int f1,f2;
            if (e->ends[0]->tet == e->ends[1]->tet) {
                // e is a loop. Check two other edges on this tet are parallel
                Tetrahedron *t = e->ends[0]->tet;
                Tetrahedron *t2 = 0;
                Tetrahedron *t3 = 0;
                // Break out of this cycle once we've found t3 (which means
                // we've also found t2)
                for(unsigned int j=0; (t3 == 0) && (j < 4); j++) {
                    if (t2 == 0) {
                        // Assign t2
                        t2 =(t->externalEdgeEnd[j]->edge->otherEnd(t->externalEdgeEnd[j]))->tet; 
                        // If t2 == t, then this edge is actually the "loop"
                        // edge.
                        if (t2 == t) {
                            t2 = 0;
                        } else {
                            f1 = j;
                        }
                    } else {
                        // Assign t3
                        t3 = (t->externalEdgeEnd[j]->edge->otherEnd(t->externalEdgeEnd[j]))->tet; 
                        // If t3 == t, then this edge is actually the "loop"
                        // edge.
                        if (t3 == t) {
                            t3 = 0;
                        } else {
                            f2 = j;
                        }
                    }
                }
                if ((t2 != 0) && (t2 == t3)) {
                    threeCycleEdges->push_back(new ThreeCycle(e,f1,f2));
                }
            }
        }
    }
    
    if (threeCyclesDone < threeCycleEdges->size()) {
        // Try the search with this 3-cycle not coloured.
        runSearch(maxDepth,threeCycleEdges,threeCyclesDone+1);
        // Colour 3 cycle


        // Extract all three edges
        Edge *e = (*threeCycleEdges)[threeCyclesDone]->loop;
        unsigned int o1 = (*threeCycleEdges)[threeCyclesDone]->otherFaces[0];
        unsigned int o2 = (*threeCycleEdges)[threeCyclesDone]->otherFaces[1];

        unsigned int i;
        Tetrahedron *t = e->ends[0]->tet;
        EdgeEnd *ee2 = t->externalEdgeEnd[o1];
        Edge *e2 = ee2->edge;
        EdgeEnd *e2o = e2->otherEnd(ee2);
        Tetrahedron *t2 = e2o->tet;
        EdgeEnd *ee3 = t->externalEdgeEnd[o2];
        Edge *e3 = ee3->edge;
        EdgeEnd *e3o = e3->otherEnd(ee3);

        nextColour++;
        cycleLengths[nextColour]=0;
        // Find smallest index edge
        if (e->index < e2->index) {
            if ( e->index < e3->index) {
                // e is smallest
                if ( e2->index < e3->index ) {
                    // e,e2,e3
                    e->colour(nextColour);
                    edgesLeft--;
                    i = 5-NEdge::edgeNumber[e->ends[1]->face][o1];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    signed dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    ee2->map[i] = e2->used;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e2o->map[i] = e2->used;
                    dir = e2->index;
                    if (e2->ends[0] == e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e3->colour(nextColour);
                    edgesLeft--;
                    e3o->map[i] = e3->used;
                    i = 5-NEdge::edgeNumber[o2][e->ends[0]->face];
                    t->internalEdges[i] = nextColour;
                    ee3->map[i] = e3->used;
                    e->ends[0]->map[i] = e->used; 
                    dir = e3->index;
                    if (e3->ends[0] != e3o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                } else {
                    // e,e3,e2
                    e->colour(nextColour);
                    edgesLeft--;
                    
                    i = 5-NEdge::edgeNumber[e->ends[1]->face][o2];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    signed dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e3->colour(nextColour);
                    edgesLeft--;
                    ee3->map[i] = e3->used;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e3o->map[i] = e3->used;
                    dir = e3->index;
                    if (e3->ends[0] == e3o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    e2o->map[i] = e2->used;
                    i = 5-NEdge::edgeNumber[o1][e->ends[0]->face];
                    t->internalEdges[i] = nextColour;
                    ee2->map[i] = e2->used;
                    e->ends[0]->map[i] = e->used; 
                    dir = e2->index;
                    if (e2->ends[0] != e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                }
            } else { 
                // e3 is smallest.
                
                if (e3->ends[0] == ee3) { 
                    // e3,e2,e
                    e3->colour(nextColour);
                    edgesLeft--;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e3o->map[i] = e3->used;
                    signed dir = e3->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    e2o->map[i] = e2->used;
                    
                    i = 5-NEdge::edgeNumber[e->ends[0]->face][o1];
                    t->internalEdges[i] = nextColour;
                    ee2->map[i] = e2->used;
                    dir = e2->index;
                    if (e2->ends[1] == e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[o2][e->ends[1]->face];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    ee3->map[i] = e3->used;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                } else {
                    // e3,e,e2
                    e3->colour(nextColour);
                    edgesLeft--;
                    i = 5-NEdge::edgeNumber[o2][e->ends[0]->face];
                    t->internalEdges[i] = nextColour;
                    ee3->map[i] = e3->used;
                    // 3-cycles never appear with shared edges in common 
                    // so we can know these values.
                    signed dir = e3->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[e->ends[1]->face][o1];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    ee2->map[i] = e2->used;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e2o->map[i] = e2->used;
                    e3o->map[i] = e3->used;
                    dir = e2->index;
                    if (e2->ends[0] == e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                }
            }
        } else { 
            // e2 < e
            if ( e2->index < e3->index) {
                // e2 smallest
                if (ee2 == e2->ends[0]) { 
                    // e2,e3,e
                    
                    e2->colour(nextColour);
                    edgesLeft--;
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e2o->map[i] = e2->used;
                    signed dir = e2->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                    
                    e3->colour(nextColour);
                    edgesLeft--;
                    e3o->map[i] = e3->used;
                    i = 5-NEdge::edgeNumber[o2][e->ends[0]->face];
                    t->internalEdges[i] = nextColour;
                    ee3->map[i] = e3->used;
                    dir = e3->index;
                    if (e3->ends[0] != e3o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[e->ends[1]->face][o1];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    ee2->map[i] = e2->used;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                } else {
                    // e2,e,e3
                    
                    e2->colour(nextColour);
                    edgesLeft--;
                    i = 5-NEdge::edgeNumber[e->ends[0]->face][o1];
                    t->internalEdges[i] = nextColour;
                    ee2->map[i] = e2->used;
                    signed dir = e2->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                    
                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[o2][e->ends[1]->face];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e3->colour(nextColour);
                    edgesLeft--;
                    ee3->map[i] = e3->used;
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e3o->map[i] = e3->used;
                    e2o->map[i] = e2->used;
                    dir = e3->index;
                    if (e3->ends[0] == e3o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                }
            } else { 
                // e3 is smallest.
                    
                if (e3->ends[0] == ee3) { 
                    // e3,e2,e
                    e3->colour(nextColour);
                    edgesLeft--;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e3o->map[i] = e3->used;
                    signed dir = e3->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    e2o->map[i] = e2->used;
                    
                    i = 5-NEdge::edgeNumber[e->ends[0]->face][o1];
                    t->internalEdges[i] = nextColour;
                    ee2->map[i] = e2->used;
                    dir = e2->index;
                    if (e2->ends[1] == e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[o2][e->ends[1]->face];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    ee3->map[i] = 1;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                } else {
                    // e3,e,e2
                    e3->colour(nextColour);
                    edgesLeft--;
                    i = 5-NEdge::edgeNumber[o2][e->ends[0]->face];
                    t->internalEdges[i] = nextColour;
                    ee3->map[i] = e3->used;
                    signed dir = e3->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e->colour(nextColour);
                    edgesLeft--;
                    e->ends[0]->map[i] = e->used; 
                    
                    i = 5-NEdge::edgeNumber[e->ends[1]->face][o1];
                    t->internalEdges[i] = nextColour;
                    e->ends[1]->map[i] = e->used;
                    dir = e->index;
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;

                    e2->colour(nextColour);
                    edgesLeft--;
                    ee2->map[i] = e2->used;
                    
                    i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
                    t2->internalEdges[i] = nextColour;
                    e2o->map[i] = e2->used;
                    e3o->map[i] = e3->used;
                    dir = e2->index;
                    if (e2->ends[0] == e2o) {
                        dir = -dir;
                    }
                    cycles[nextColour][cycleLengths[nextColour]] = dir;
                    cycleLengths[nextColour]++;
                }
            }
        }

        if (isCanonical()) {
            runSearch(maxDepth,threeCycleEdges,threeCyclesDone+1);
        }
        
        // Uncolour the 3 edges
        e->unColour();
        e2->unColour();
        e3->unColour();
        edgesLeft+=3;


        // Clear internal edges and edge-end mappings

        // Note that below we "clear" 5 internal edges, yet we know only
        // 3 were used. The two extras we know cannot have been used in any
        // cycles since each 3-cycle in our decomposition cannot share the same
        // main tetrahedra, so clearing them should not have any effect
        i = 5-NEdge::edgeNumber[e->ends[0]->face][o1];
        t->internalEdges[i] = 0;
        e->ends[1]->map[i] = 0;
        ee2->map[i] = 0;
        
        i = 5-NEdge::edgeNumber[e->ends[0]->face][o2];
        t->internalEdges[i] = 0;
        ee3->map[i] = 0;
        e->ends[0]->map[i] = 0; 
       
        i = 5-NEdge::edgeNumber[e->ends[1]->face][o1];
        t->internalEdges[i] = 0;
        ee3->map[i] = 0;
        e->ends[0]->map[i] = 0; 
       
        i = 5-NEdge::edgeNumber[e->ends[1]->face][o2];
        t->internalEdges[i] = 0;
        ee3->map[i] = 0;
        e->ends[0]->map[i] = 0; 
       
        i = 5-NEdge::edgeNumber[e2o->face][e3o->face];
        t2->internalEdges[i] = 0;
        e2o->map[i] = 0;
        e3o->map[i] = 0;



        cycleLengths[nextColour]=0;
        nextColour--;
    } else {
        colourLowestEdge();
    }
    // This function gets called recursively. Ensure that we "complete" the 
    // search only if we are completely finished.
    if ( threeCyclesDone == 0) {
        for(std::vector<ThreeCycle*>::iterator iter = threeCycleEdges->begin();
                iter != threeCycleEdges->end(); iter++ ) {
            delete *iter;
        }
        delete threeCycleEdges;
        use_(0, useArgs_);
    }
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
    delete[] simp;
    return ans;
}

void CycleDecompSearcher::dumpData(std::ostream& out) const {
    //NCompactSearcher::dumpData(out);

    unsigned int i;
    unsigned int j;
    for (i = 1; i <= nCycles && cycleLengths[i] > 0; i++) {
        out << i << ": ";
        for(j = 0; j < cycleLengths[i]; j++) {
            out << cycles[i][j];
            if (j != cycleLengths[i]-1) 
                out << ", ";
        }
        out << std::endl;
    }
    //for (i = 0; i < nEdges; i++) {
    //    out << i << ": " << edges[i].colours[0]
    //        << ", " << edges[i].colours[1]
    //        << ", " << edges[i].colours[2] << std::endl;
    //}

    for (i = 0; i < nTets; i++) {
        for(j = 0; j < 6; j++) {
            out << tets[i].internalEdges[j] << " ";
        }
        if (i != (nTets-1))
            out << ": ";
    }
    out << std::endl;
    out << "------------------" << std::endl;
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
    if ( automorphisms == 0 ) {
        return true;
    }
    unsigned int autoNo;
    bool thisAuto;
    unsigned int *cycleList[nextColour+1];
    unsigned int cycleListLengths[nextColour+1];
    unsigned int offset[nCycles];
    //for(unsigned int i=1; i<=nextColour; i++) {
    //    cycleList[i] = new unsigned int[3*nEdges];//cycleLengths[i]];
    //    //std::cout << "cycleList["<<i<<"] has size " << cycleLengths[i];
    //    //std::cout << " At " << cycleList[i] << std::endl;
    //}
    bool debug = false;
    if (debug) {
        std::cout << "CANON Test" << std::endl;
    }
    //if ( nextColour == 5 && 
    //     (cycleLengths[1] == 1) &&
    //     (cycleLengths[2] == 8) &&
    //     (cycleLengths[3] == 9) &&
    //     (cycleLengths[4] == 5) &&
    //     (cycleLengths[5] == 1) && 
    //     (cycles[1][0] == 1) &&
    //     (cycles[2][0] == 1) &&
    //     (cycles[2][1] == 2) &&
    //     (cycles[2][2] == 4) &&
    //     (cycles[2][3] == -3) &&
    //     (cycles[2][4] == -1) &&
    //     (cycles[2][5] == 3) &&
    //     (cycles[2][6] == -5) &&
    //     (cycles[2][7] == -2) &&
    //     (cycles[3][0] == 2) &&
    //     (cycles[3][1] == 6) &&
    //     (cycles[3][2] == -7) &&
    //     (cycles[3][3] == -5) &&
    //     (cycles[3][4] == 4) &&
    //     (cycles[3][5] == 7) &&
    //     (cycles[3][6] == 8) &&
    //     (cycles[3][7] == -7) &&
    //     (cycles[3][8] == -3) &&
    //     (cycles[4][0] == 4) &&
    //     (cycles[4][1] == -5) &&
    //     (cycles[4][2] == 6) &&
    //     (cycles[4][3] == 8) &&
    //     (cycles[4][4] == -6) &&
    //     (cycles[5][0] == 8)
    //     ) {
    //         debug=true;
    //}

    for(autoNo=0; autoNo < nAutos; autoNo++) {
        Automorphism *theAuto = automorphisms[autoNo];
        // As long as thisAuto is true, we keep looking at this automorphism.
        // If it becomes false, it means that this automorphism results in a
        // less canonical representation.
        thisAuto=true;
        // Generate new cycle lists
        for(unsigned int i=1; i<nextColour; i++) {
            cycleList[i] = theAuto->cycles[i];
            cycleListLengths[i] = theAuto->cycleLength[i];
            offset[i] = theAuto->offset[i];
        }
        offset[nextColour]=0;
        cycleList[nextColour] = theAuto->cycles[nextColour];
        // offset[i] is going to hold the current "starting" place for the
        // cycle.

        unsigned int nextEdge = 255;
        // nextEdge will store the "second" edge in the cycle, if offset[i]
        // is the first edge.  Note that the "second" edge may be the one
        // before the first edge, if the first edge is negative (since
        // swapping signs also reverses direction).
        // A value of 255 means something is broken.
        
        unsigned int newEdge = (*automorphisms[autoNo])[cycles[nextColour][0]];
        theAuto->cycles[nextColour][0] = newEdge;
        unsigned int min = newEdge;

        bool checkNextPair=false;
        bool setNextEdge=false;
        if ( newEdge %2 == 0) {
            // first edge +ve
            setNextEdge = true;
        } else {
            nextEdge =(*automorphisms[autoNo])[cycles[nextColour][cycleLengths[nextColour]-1]];
        }
        if (debug) std::cout << "Pos: 0   Offset: 0   Start: " << min << " Next: " << nextEdge << std::endl; 
        for(unsigned int j=1; j < cycleLengths[nextColour]; j++) {
            if (debug) std::cout << "Pos: "<<j<<"   Offset: "<<offset[nextColour]<<"   Start: " << min << " Next: " << nextEdge << std::endl; 
            newEdge = (*automorphisms[autoNo])[cycles[nextColour][j]];

            // If checkNextPair is true, that means the last edge we
            // checked was equal-smallest.  We check this new edge against
            // the "next" edge from the current smallest, and if we have a
            // smaller edge, update the offset.
            if (checkNextPair) {
                if (newEdge < nextEdge) {
                    offset[nextColour]=j-1;
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
                if ((min%2 == 1) && ( newEdge + 1 == min )) {
                    // Current lowest is negative, new one is positive.
                    // Need to check "next".
                    if ( j == cycleLengths[nextColour]-1) {
                        // "next" edge is first
                        unsigned int newNext = cycleList[nextColour][0];
                        if  ( newNext < nextEdge) {
                            offset[nextColour]=j;
                            // no need to change nextEdge
                        }

                    } else {
                        // We need to check the next edge against the
                        // "next edge" from the current minimum, so
                        // make this note.
                        checkNextPair = true;
                    }
                } else {
                    min = newEdge;
                    offset[nextColour] = j;
                    if (min%2 ==1) {
                        // Negative
                            
                        // If i%2 == 0, then 1 - 2*(i%2) == 1 
                        // If i%2 == 1, then 1 - 2*(i%2) == -1
                        //
                        // So the line below adds one (if i%2 ==0) or
                        // subtracts one (if i%2 == 1), which flips the
                        // sign.
                        nextEdge = cycleList[nextColour][j-1] + (1- 2*(cycleList[nextColour][j-1]%2));
                    } else {
                        // Positive
                        setNextEdge = true;
                    }
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
                        unsigned int newNext = cycleList[nextColour][j-1] + (1 - 2*(cycleList[nextColour][j-1]%2));

                        if (newNext < nextEdge) {
                                offset[nextColour]=j;
                                nextEdge = newNext;
                                // min stays the same.
                        }
                    } else { 
                            // min and current edges are positive. 

                        if ( j == cycleLengths[nextColour]-1) {
                            // "next" edge is first
                            unsigned int newNext = cycleList[nextColour][0];
                            if  ( newNext < nextEdge) {
                                offset[nextColour]=j;
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
                        unsigned int newNext = cycleList[nextColour][j-1] + (1 - 2*(cycleList[nextColour][j-1]%2));
                        if ( newNext < nextEdge) {
                            min = newEdge;
                            nextEdge = newNext;
                            offset[nextColour]=j;
                        }
                    }
                }
            }
            theAuto->cycles[nextColour][j] = newEdge;
        }
        theAuto->cycleLength[nextColour] = cycleLengths[nextColour];
        theAuto->offset[nextColour] = offset[nextColour];
        cycleList[nextColour] = theAuto->cycles[nextColour];;
        cycleListLengths[nextColour] = cycleLengths[nextColour];
        if (false && debug) {
            std::cout << "Before sorting" << std::endl;
            for(unsigned int k=1; k<=nextColour; k++) {
                for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                    signed int val = cycleList[k][l];
                    if ( l == offset[k] ) {
                        std::cout << "(";
                    }
                    if (val % 2 == 1) {
                        val = - ( val-1)/2;
                    } else {
                        val = val/2;
                    }
                    std::cout << val;
                    if ( l == offset[k] ) {
                        std::cout << ")";
                    }
                    std::cout << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "---" << std::endl;
        }

        // Sort cycleList based on values of cycleList[i][offset[i]] 

        // Note that 3-cycles must come first in this implementation so we
        // must sort 3-cycles first.
        
        
        unsigned int i=1; 
        // Can probably do this better.  This is actually one loop split into
        // two. I decided to split it to ease readability.
        for (;thisAuto &&i <= nextColour && (cycleListLengths[i] == 3) ; i++) {
            // Reverse bubble sort order so position[0] is definitely lowest
            // after one run through.  We can then check to see if this is more or
            // less canonical after 1 run through.
            for(unsigned int j=nextColour; j > i; j--) {
                if ( cycleListLengths[j] != 3 ) {
                    continue;
                }
                // Since we're working on 3-cycles (and >= 5 tets) , we know 
                // the 'first/lowest' edges can never be the same so we only
                // have to compare them.
                if ( cycleList[j][offset[j]] < cycleList[j-1][offset[j-1]]) {
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
            // Check the 3 edges
            signed int add = (cycleList[i][offset[i]] % 2 == 0) ? 1 : -1;
            unsigned int off = offset[i];
            for (unsigned int j=0;j<3;j++) {

                signed int cycleVal = 2*cycles[i][j];

                if (cycleVal < 0) {
                    cycleVal = (-cycleVal) + 1;
                }

                // Get the value of the edge in the automorphism
                unsigned newE = cycleList[i][off];

                // If the offset is negative, flip sign on this edge.
                if (add == -1) {
                    // If i%2 == 0, then 1 - 2*(i%2) == 1 
                    // If i%2 == 1, then 1 - 2*(i%2) == -1
                    //
                    // So the line below adds one (if i%2 ==0) or
                    // subtracts one (if i%2 == 1), which flips the
                    // sign.
                    newE += (1 - 2*(newE%2));
                }
                if ( cycleVal > newE ) {
                    return false;
                }
                if ( cycleVal < newE ) {
                    thisAuto=false;
                    break;
                }
                if  ( add == 1) {
                    off+=1;
                    if (off == cycleListLengths[i] ) {
                        off = 0;
                    }
                } else {
                    // add == -1
                    off-=1;
                    if (off == -1) {
                        off = cycleListLengths[i]-1;
                    }
                }
            }
        }
        for(; thisAuto && i <= nextColour ; i++ ) {
            // Do some sorting!
            assert(cycleListLengths[i] > 3);
            

            // Reverse bubble sort order so position[0] is definitely lowest
            // after one run through.  We can then check to see if this is more or
            // less canonical after 1 run through.
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
            if (debug) {
                std::cout << "Sorted " << i << " time(s)" << std::endl;
                for(unsigned int k=1; k<=nextColour; k++) {
                    for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                        signed int val = cycleList[k][l];
                        if ( l == offset[k] ) {
                            std::cout << "(";
                        }
                        if (val % 2 == 1) {
                            val = - ( val-1)/2;
                        } else {
                            val = val/2;
                        }
                        std::cout << val;
                        if ( l == offset[k] ) {
                            std::cout << ")";
                        }
                        std::cout << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "---" << std::endl;
            }
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
                    if (debug) {
                        std::cout << "Checking" << std::endl;
                        for(unsigned int k=1; k<=nextColour; k++) {
                            for(unsigned int l=0; l < cycleLengths[k]; l++) {
                                std::cout << cycles[k][l] << " ";
                            }
                            std::cout << std::endl;
                        }
                        std::cout << "---" << std::endl;
                        std::cout << " bigger on " << j << std::endl;
                        for(unsigned int k=1; k<=nextColour; k++) {
                            for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                                signed int val = cycleList[k][l];
                                if ( l == offset[k] ) {
                                    std::cout << "(";
                                }
                                if (val % 2 == 1) {
                                    val = - ( val-1)/2;
                                } else {
                                    val = val/2;
                                }
                                std::cout << val;
                                if ( l == offset[k] ) {
                                    std::cout << ")";
                                }
                                std::cout << " ";
                            }
                            std::cout << std::endl;
                        }
                        std::cout << "---" << std::endl;
                    }
                    //for(unsigned int k=1; k<=nextColour; k++)
                    //    delete[] cycleList[k];
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
                if (debug) {
                    std::cout << "Beaten by length at cycle " << i << std::endl;
                    for(unsigned int k=1; k<=nextColour; k++) {
                        for(unsigned int l=0; l < cycleListLengths[k]; l++) {
                            signed int val = cycleList[k][l];
                            if ( l == offset[k] ) {
                                std::cout << "(";
                            }
                            if (val % 2 == 1) {
                                val = - (val-1)/2;
                            } else {
                                val = val/2;
                            }
                            std::cout << val;
                            if ( l == offset[k] ) {
                                std::cout << ")";
                            }
                            std::cout << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << "---" << std::endl;
                }
                //for(unsigned int k=1; k<=nextColour; k++)
                //    delete[] cycleList[k];
                return false;
            }
        }
    }
    //for(unsigned int i=1; i<=nextColour; i++)
    //    delete[] cycleList[i];
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
    assert(0 <= used);
    assert(used < 3);
    colours[used++] = newColour;
}

inline void CycleDecompSearcher::Edge::unColour() {
    assert(0 < used);
    assert(used <= 3);
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
        const Edge *edges, const unsigned int _nEdges, 
        const unsigned int _nCycles) : nEdges(_nEdges), nCycles(_nCycles)  {
    edgeMap = new unsigned int[2*(nEdges+1)];
    realEdgeMap = edgeMap+nEdges;
    cycles = new unsigned int*[nCycles+1];
    cycleLength = new unsigned int[nCycles+1];
    offset = new unsigned int[nCycles+1];
    for(unsigned int i=0;i< nCycles+1;i++) {
        cycles[i] = new unsigned int[3*nEdges];
    }
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
                unsigned int jIndex = 2*edges[j].index;
                realEdgeMap[iIndex] = jIndex;
                realEdgeMap[std::ptrdiff_t(-iIndex)] = jIndex+1;
                break;
            }
            if (( newEndTet == edges[j].ends[0]->tet->index) &&
                ( newEndFace == edges[j].ends[0]->face) &&
                ( newStartTet == edges[j].ends[1]->tet->index) &&
                ( newStartFace == edges[j].ends[1]->face)) {
                unsigned int jIndex = 2*edges[j].index;
                realEdgeMap[std::ptrdiff_t(-iIndex)] = jIndex;
                realEdgeMap[iIndex] = jIndex+1;
                break;
            }
        }
    }
}

CycleDecompSearcher::Automorphism::~Automorphism() {
    for(unsigned int i=0;i< nCycles+1;i++) {
        delete[] cycles[i];
    }
    delete[] cycles;
    delete[] cycleLength;
    delete[] offset;
    delete[] edgeMap;
}

unsigned int inline CycleDecompSearcher::Automorphism::operator [] (const signed int in) {
  //return edgeMap[in+nEdges];
  return realEdgeMap[(std::ptrdiff_t)(in)];
}

CycleDecompSearcher::ThreeCycle::ThreeCycle(Edge *e, 
        unsigned int f1, unsigned int f2) {
    loop = e;
    otherFaces[0] = f1;
    otherFaces[1] = f2;
}

/* vim: set ts=4 sw=4 expandtab: */
