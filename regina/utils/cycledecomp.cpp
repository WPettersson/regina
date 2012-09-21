
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


    if ( autos ) {
        nAutos = 0;
        automorphisms = new Automorphism*[autos->size()];
        for( NFacePairing::IsoList::const_iterator it = autos->begin();
                it != autos->end(); it++) {
            automorphisms[nAutos] = new Automorphism( (*it),edges,nEdges);
            nAutos++;
        }
        //std::cout << "Original " << std::endl;

        //cycleLengths[1] = 3;
        //cycles[1][0] = -2;
        //cycles[1][1] = 3;
        //cycles[1][2] = 5;
        //cycles[1][3] = 7;
        //std::cout << "1: ";
        //for(unsigned int i=0; i < cycleLengths[1]; i++) {
        //    std::cout << cycles[1][i];
        //    if (i < (cycleLengths[1]-1))
        //        std::cout << ", ";
        //}
        //std::cout << std::endl;
        //std::cout << "Automorphisms " << std::endl;
        //
        //signed int *newCycle;
        //newCycle = new signed int[cycleLengths[1]];
        //signed int *internalEdges;
        //internalEdges = new signed int[6*nTets];
        //for (unsigned int autoCounter=0; autoCounter< nAutos; autoCounter++) {
        //    for(unsigned int i=1; i <= nEdges; i++) 
        //        std::cout << (*automorphisms[autoCounter])[i] << " ";
        //    std::cout << "::\t";
        //    for(unsigned int i=0; i < cycleLengths[1]; i++) {
        //        newCycle[i] = (*automorphisms[autoCounter])[cycles[1][i]];
        //        std::cout << newCycle[i];
        //        if (i < (cycleLengths[1]-1))
        //            std::cout << ", ";
        //    }
        //    //std::cout << std::endl;
        //    std::cout << "::\t";
        //    for(unsigned int i=0; i < 6*nTets; i++) 
        //        internalEdges[i]=0;

        //    for(unsigned int i=1; i < cycleLengths[1]; i++) {
        //        signed int edge = newCycle[i];
        //        signed int lastEdge = newCycle[i-1];
        //        unsigned int tet;
        //        unsigned int oneFace;
        //        unsigned int otherFace;
        //        if (lastEdge < 0) {
        //            tet = edges[(-lastEdge)-1].ends[0]->tet->index;
        //            oneFace = edges[-(lastEdge-1)].ends[0]->face;
        //        } else {
        //            tet = edges[lastEdge-1].ends[1]->tet->index;
        //            oneFace = edges[lastEdge-1].ends[1]->face;
        //        }
        //        if (edge < 0) {
        //            otherFace = edges[(-edge)-1].ends[1]->face;
        //        } else {
        //            otherFace = edges[edge-1].ends[0]->face;
        //        }
        //        assert((0 <= tet) && (tet < 6));
        //        // NEdge::edgeNumber[oneFace][otherFace] gives the edge joining
        //        // vertices oneFace and otherFace (as vertices).  We want the
        //        // exact opposite edge, so we do 5 - this value.
        //        internalEdges[(6*tet)+(5-NEdge::edgeNumber[oneFace][otherFace])] = 1;
        //    }
        //    for(unsigned int i=0; i < 6*nTets; i++) {
        //      std::cout << internalEdges[i] << " ";
        //      if ( (i%6) == 5)
        //          std::cout << ": ";
        //    }
        //    std::cout << std::endl;
        //}
        //delete[] internalEdges;
        //delete[] newCycle;
        //cycleLengths[1] = 0;
        //cycles[1][0] = 0;
        //cycles[1][1] = 0;
        //cycles[1][2] = 0;
        //cycles[1][3] = 0;
    } else {
        automorphisms = 0;
    }
}

CycleDecompSearcher::~CycleDecompSearcher() {
    for( unsigned int i=0; i < nTets+1; i++ ) {
        delete[] cycles[i];
    }
    for( unsigned int i=0; i < nAutos; i++ ) {
        delete[] automorphisms[i];
    }
    delete[] automorphisms;
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
    // Work out which edge should be last (used to map automorphisms).
    // If ends[0] of this edge is last visited, then we're going backwards
    // along this edge.
    if ( start->edge->ends[0] == start )
        lastEdge = -1*(start->edge->index);
    else
        lastEdge = start->edge->index;


    unsigned int outFace = NEdge::edgeVertex[5-edge][1];
    EdgeEnd *outEdgeEnd = tets[tet].externalEdgeEnd[outFace];
    Edge *nextEdge = outEdgeEnd->edge;
    nextEdge->colour(nextColour);
    outEdgeEnd->map[edge] = nextEdge->used;
    edgesLeft--;
    EdgeEnd *nextEdgeEnd = nextEdge->otherEnd(outEdgeEnd);
    
    std::cout << "t"<<tet<<"i"<<edge<<"l"<<lastEdge<<"n"<<nextEdge->index << std::endl;
    
    // An edge has two ends, ends[0] and ends[1]
    // A "forward" direction along an edge is from ends[0] to ends[1]
    // and "backwards" is the reverse. Forwards is denoted by a positive
    // edge number and backwards by a negative edge number.
    //
    // If outEdgeEnd == ends[1] that means we're using the edge 
    // ends[1] -> ends[0]
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
      
        if (isCanonical()) {
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
            // Try to find more paths.
            nextPath(start, firstEdge, nextEnd);
        }

        

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

bool CycleDecompSearcher::isCanonical() {
    unsigned int *internalEdges = new unsigned int[6*nTets];
    unsigned int *mapping = new unsigned int[nextColour+1]; // Colours start from 1
    unsigned int autoNo, cycleNo;
    unsigned int start;
    Tetrahedron *tet;
    unsigned int tetIndex;
    bool iso;
    for(autoNo=0; autoNo < nAutos; autoNo++) {
        iso = true; 
        for(unsigned int i=0; i < 6*nTets; i++) 
            internalEdges[i]=0;
        for(cycleNo = 1; cycleNo <= nextColour; cycleNo++) {
            mapping[cycleNo] = 0;
            for(unsigned int i=0; i < cycleLengths[cycleNo]; i++) {
                signed int edge;
                if ( i == 0) {
                    if ( cycleNo == nextColour ) {
                        // Last cycle
                        edge = lastEdge;
                    } else {
                        edge = cycles[cycleNo][cycleLengths[cycleNo]-1];
                    }
                } else {
                    edge = cycles[cycleNo][i-1];
                }
                signed int nextEdge = cycles[cycleNo][i];
                unsigned int tet;
                unsigned int oneFace;
                unsigned int otherFace;
                if (edge < 0) {
                    tet = edges[(-edge)-1].ends[0]->tet->index;
                    oneFace = edges[-(edge)-1].ends[0]->face;
                } else {
                    tet = edges[edge-1].ends[1]->tet->index;
                    oneFace = edges[edge-1].ends[1]->face;
                }
                if (nextEdge < 0) {
                    otherFace = edges[(-nextEdge)-1].ends[1]->face;
                } else {
                    otherFace = edges[nextEdge-1].ends[0]->face;
                }
                // NEdge::edgeNumber[oneFace][otherFace] gives the edge joining
                // vertices oneFace and otherFace (as vertices).  We want the
                // exact opposite edge, so we do 5 - this value.
                if (!((6*tet)+(5-NEdge::edgeNumber[oneFace][otherFace]) < 6*nTets)) {
                    std::cout << "Blah";
                }
                assert((6*tet)+(5-NEdge::edgeNumber[oneFace][otherFace]) < 6*nTets);
                assert(0 <= (6*tet)+(5-NEdge::edgeNumber[oneFace][otherFace]));

                internalEdges[(6*tet)+(5-NEdge::edgeNumber[oneFace][otherFace])] = cycleNo;
            }
        }
        //for(unsigned int i=0; i < 6*nTets; i++) {
        //  std::cout << internalEdges[i] << " ";
        //  if ( (i%6) == 5)
        //      std::cout << ": ";
        //}
        //std::cout << std::endl;
        // Check to see if the first "nextColour-1" cycles can be matched up
        // perfectly. Note that we assume by induction that these cycles are
        // all canonical already.
        // Note that this loop terminates if iso=false (no longer isomorphic).
        for ( tetIndex=0; iso && tetIndex < nTets; tetIndex++ ) {
            for( unsigned int i=0; i < 6; i++) {
                unsigned int currentNumber = tets[tetIndex].internalEdges[i];
                unsigned int isoNumber = internalEdges[(6*tetIndex)+i];
                // If it's the current colour we're looking at, skip it.
                if (currentNumber == nextColour)
                    continue;
                if (isoNumber == nextColour)
                    continue;
                // If edge uncoloured in current graph, then it must be
                // uncoloured in isomorphism.
                if (currentNumber == 0) {
                    if ( isoNumber == 0) 
                        continue;
                    else {
                        iso = false;
                        break;
                    }
                }
                // So edge is coloured.  If no mapping is set up, and the
                // corresponding edge in the iso is coloured, set up a mapping.
                // If corresponding edge not coloured, break out
                // Else, just check that the corresponding edge makes 
                if ( mapping[currentNumber] == 0 ) {
                    if (isoNumber == 0) {
                        iso = false;
                        break;
                    } else {
                        mapping[currentNumber] = isoNumber;
                    }
                } else {
                    if ( mapping[currentNumber] != isoNumber ) {
                        iso = false;
                        break;
                    }
                }
            }
        }
        if (! iso) 
            continue;
        // Is isomorphic so far, check the last cycle
        // First find the "start" of the cycle so we can keep that fixed.
        // Remember, the start of the cycle is the "lowest free edge" on the
        // tetrahedra with the most used edges (and we shouldn't cound an edge
        // as used if it has nextColour on it, as that is part of the current
        // cycle.
        tetIndex=0;
        unsigned int mostUsed=0;
        unsigned int currentUsed;
        for( unsigned int i=0; i < nTets; i++) {
            currentUsed=0; 
            for( unsigned int j=0; j < 6; j++) {
                if ((tets[i].internalEdges[j] != 0) && 
                    (tets[i].internalEdges[j] != nextColour)) {
                    currentUsed+=1;
                }
            }
            if ((currentUsed < 6) && (currentUsed > mostUsed)) {
                mostUsed = currentUsed;
                tetIndex=i;
            }
        }
            
        tet = &(tets[tetIndex]);
        start = 0;
        // Note that since this tetrahedra had the most used edges, we must
        // have started colouring on it, so it must have an edge coloured
        // nextColour.
        while ( tet->internalEdges[start] != nextColour) {
            start++;
            assert(start < 6);
        }
        // Check the rest of the tetrahedra. 
        // I internalEdges[] is nextColour before
        // tet[start]->internalEdges then this current colouring is not
        // canonical. 
        //
        // If tet[start]->internalEdges is nextColour before
        // internalEdges, then this current colouring is definitely more
        // canonical than this isomorphism.
        
        // If possiblyMoreCanonical is set, that means the isomorphic graph is 
        // possibly more canonical, so keep checking.
        bool possiblyMoreCanonical=true;
        for( unsigned int i=0; possiblyMoreCanonical && i < nTets; i++) {
            for( unsigned int j=0; j < 6; j++) {
                if  ( internalEdges[(6*i)+j] == nextColour) {
                    if ( tets[i].internalEdges[j] != nextColour) 
                        return false;
                }
                if ( tets[i].internalEdges[j] == nextColour) {
                    if ( internalEdges[(6*i)+j] != nextColour) {
                        possiblyMoreCanonical=false;
                        break;
                    }
                }
            }
        }
    }
    delete[] internalEdges;
    delete[] mapping;
    return true;
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

CycleDecompSearcher::Automorphism::Automorphism(const NIsomorphism * iso,
        const Edge *edges, const unsigned int _nEdges) {
    nEdges = _nEdges;
    edgeMap = new signed int[nEdges+1];
    for (unsigned int i=0; i< nEdges;i++) {
        unsigned int startFace = edges[i].ends[0]->face;
        unsigned int startTet = edges[i].ends[0]->tet->index;
        unsigned int endFace = edges[i].ends[1]->face;
        unsigned int endTet = edges[i].ends[1]->tet->index;

        unsigned int newStartFace = iso->facePerm(startTet)[startFace];
        unsigned int newStartTet = iso->tetImage(startTet);
        unsigned int newEndFace = iso->facePerm(endTet)[endFace];
        unsigned int newEndTet = iso->tetImage(endTet);
        
        for(unsigned int j=0; j< nEdges;j++) {
            if (( newStartTet == edges[j].ends[0]->tet->index) &&
                ( newStartFace == edges[j].ends[0]->face) &&
                ( newEndTet == edges[j].ends[1]->tet->index) &&
                ( newEndFace == edges[j].ends[1]->face)) {
                assert((0 <= j) && (j < nEdges));
                edgeMap[i+1] = j+1;
                break;
            }
            if (( newEndTet == edges[j].ends[0]->tet->index) &&
                ( newEndFace == edges[j].ends[0]->face) &&
                ( newStartTet == edges[j].ends[1]->tet->index) &&
                ( newStartFace == edges[j].ends[1]->face)) {
                assert((0 <= j) && (j < nEdges));
                edgeMap[i+1] = -(j+1);
                break;
            }
        }
    }
}

signed int CycleDecompSearcher::Automorphism::operator [] (const signed int in) {
  signed int sign = ( in < 0 ? -1 : 1);
  return sign * edgeMap[ abs(in) ];
}
