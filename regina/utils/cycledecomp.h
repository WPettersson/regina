
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2011, Ben Burton                                   *
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

/*! \file cycledecomp.h
 *  \brief Supports searching face pairings of tetrahedra via cycle 
 *  decompositions of the face pairing graph.
 */

#ifndef __CYCLEDECOMPSEARCHER_H
#ifndef __DOXYGEN
#define __CYCLEDECOMPSEARCHER_H
#endif

#include "regina-core.h"
#include "census/ngluingperms.h"
#include "census/nfacepairing.h"
#include "utilities/nqitmask.h"


class CycleDecompSearcher;
/**
 * A routine used to do arbitrary processing upon a particular set of
 * tetrahedron gluing permutations.  Such routines are used to process
 * permutation sets found when running NGluingPermSearcher::findAllPerms().
 *
 * The first parameter passed will be a set of gluing permutations
 * (in fact it will be of the subclass NGluingPermSearcher in order to
 * support partial searches as well as full searches).  This set of
 * gluing permutations must not be deallocated by this routine, since it
 * may be used again later by the caller.  The second parameter may contain
 * arbitrary data as passed to either NGluingPerms::findAllPerms() or
 * the NGluingPermSearcher class constructor.
 *
 * Note that the first parameter passed might be \c null to signal that
 * gluing permutation generation has finished.
 */
typedef void (*UseCycles)(const CycleDecompSearcher*, void*);


/**
 * A gluing permutation search class that turns the permutation search problem
 * into a cycle decomposition problem. 
 *
 * \ifacespython Not present.
 */
class CycleDecompSearcher {
    private:
        static const unsigned faceEdges[4][3];
            /**< Lists the edges contained within a given face in lexicographic
             *   order. */
        static const unsigned otherFace[4][3];
            /**< Lists the other faces that are adjacent to a given face based
             *   upon the edge they have in common.  The index into faceEdges[]
             *   is used as an index into otherFace. For example,
             *   otherface[i][j] is the face that contains edge faceEdges[j]
             *   and is adjacent to face j. */
        static const signed edgeParity[6][6];
            /**< Lists which edges should be identified, in terms of being
             *   "above" or "below" some given edge.  For example, if edge i is
             *   considered to be "above" edge j, then edge edgeParity[j][i]
             *   should also be "above" edge j. */
 
        static const signed otherVert[4][6];
            /**< Identifies the vertex on the face of a tetrahedron that is 
             *   not incident upon a given edge.  For example, the vertex which
             *   is on face i but is not on edge j is otherVertex[i][j]. */

        UseCycles use_;
            /**< A routine to call each time a gluing permutation set is
                 found during the search. */
        void* useArgs_;
            /**< Additional user-supplied data to be passed as the second
                 argument to the \a use_ routine. */
        const regina::NFacePairing *pairing_;
            /**< The face pairing graph. */

        class EdgeEnd;

        class Tetrahedron {
            /**< A tetrahedron, as represented in the face pairing graph. */
            public:
                Tetrahedron();
                    /**< Constructor. */
                unsigned int used;
                    /**< The number of internal edges used. */
                unsigned int index;
                    /**< An identifying number for this tetrahedron. */
                unsigned int internalEdges[6];
                    /**< The 6 internal edges of the tetrahedron. */
                EdgeEnd* externalEdgeEnd[4];
                    /**< The EdgeEnds which are attached to each face of the
                     *   tetrahedron. */

        };

        class Edge {
            /**< Represents a single edge of the face pairingr graph. */
            public:
                Edge(); 
                    /**< Constructor. */
                void colour(unsigned int newColour);
                    /**< Adds the "newColour" cycle to this edge. */
                void unColour();
                    /**< Removes the most recent cycle from this edge. */
                EdgeEnd* otherEnd(EdgeEnd* one);
                    /**< Given EdgeEnd one belonging to this Edge, return
                     *   the other EdgeEnd. */

                signed int colours[3];
                    /**< The 3 cycles an edge may have. */
                int used;
                    /**< The number of cycles currently using this edge. */
                //unsigned int endTet[2];
                    /**< The tetrahedron on either end  of this edge. */
                EdgeEnd* ends[2];
                    /**< The edge-ends of this edge. */
                unsigned int index;
                    /**< The index of this edge, such that 
                     *   Edge == edges[Edge.index]. */
        };

        class EdgeEnd {
            /**< The "end" of an edge in the face pairing graph. This end
             *   contains information about how the edge is attached to the
             *   tetrahedron. */
            public:
                EdgeEnd();
                    /**< Constructor. */
                Tetrahedron* tet;
                    /**< The tetrahedron attached here. */
                Edge *edge;
                    /**< The edge containing this end. */
                unsigned int face;
                    /**< The face of the tetrahedron which this edge is
                     *   attached to. */
                signed int map[6];
                    /**< A mapping of the edges of the face to the cycles on
                     *   the edge. */
        };

        class Automorphism {
            /**< This class represents an automorphism of the face pairing
             *   graph, in terms of a mapping between the edges of the graph
             *   and the edges of the automorphism. */
            public:
                Automorphism( const regina::NIsomorphism * iso, 
                              const Edge* edges,
                              const unsigned int _nEdges);
                    /**< Constructor from regina type isomorphism. */
                ~Automorphism();
                    /**< Destructor. */
                unsigned int operator [] ( const signed int initial);
                    /**< Returns the equivalent edge under the automorphism in
                     *   question. */
                void tetAndInt(unsigned int *newTet, unsigned int *newInternal,
                               unsigned int oldTet, unsigned int oldInternal);
                    /**< As above, returns the equivalent tet/internal_edge
                     *   pair under this automorphism. */
            private:
                unsigned int nEdges;
                    /**< Number of edges. */
                unsigned int *edgeMap;
                    /**< An array storing the relationship between edges.  
                     *   These are stored such that for edge e,
                     *   edgeMap[nEdge + e] gives the mapping.
                     *   That is, negative edge mappings are stored in the
                     *   range edgeMap[0] .. edgeMap[nEdge-1] and positive edge
                     *   mappings are stored in the range edgeMap[nEdge+1] ..
                     *   edgeMap[2*nEdge].*/
                unsigned int *realEdgeMap;
                    /**< An array storing the relationship between edges. 
                     *   This is simply a pointer to edgeMap[nEdge] so pointer
                     *   arithmetic is simpler when accessing mappings. */
                //unsigned int *newTets;
                    /**< An array storing how tetrahedra change under this
                     *   automorphism. */
                //unsigned int **newInts;
                    /**< An array of arrays, storing how internal edges change
                     *   under this  automorphism. For example, newInts[i][j]
                     *   will return a value from 0 to 5 inclusive, telling us
                     *   that internal edge j, of tetrahedron i, will be mapped
                     *   to internal edge "newInts[i][j]" in some tetrahedron. */
        };

        Tetrahedron *tets;
            /**< The tetrahedron representations in the face pairing graph. */
        unsigned int nTets;
            /**< The number of tetrahedra. */
        Edge *edges;
            /**< The edges of the face pairing graph. */
        unsigned int *firstEdge;
            /**< The first internal edge to be used 
             *   in each cycle. */
        unsigned int *firstOtherFace;
            /**< The first face to be used 
             *   in each cycle, after using the internal edge firstEdge.
             *   That is, follow the firstEdge throgh a tetrahedron, go out
             *   through the face NEdge::edgeVertex[5-firstEdge][0], and follow
             *   that new edge.  This new edge goes into face firstOtherFace on
             *   the new tetrahedron. */
        unsigned int nCycles;
            /**< The maximum number of cycles to search for. */
        unsigned int nEdges;
            /**< The number of edges. */
        EdgeEnd *ends;
            /**< The ends of the edges. */
        unsigned int nEnds;
            /**< The number of edge ends. */
        unsigned int nextColour;
            /**< The next "colour" to be used to mark out a cycle in the 
             *   face pairing graph. */
        unsigned int edgesLeft;
            /**< How many edges have not been coloured yet. */
        signed int ** cycles;
            /**< Contains a list of all edge numbers used in each cycle.
             *   cycles[x][y] denotes the y'th edge in cycle number x. 
             *   Some edge numbers are negative, this indicates following
             *   the edge in the "reverse" direction (end[1] to end[0]).
             *   These are the spine codes as Matveev uses them. */
        unsigned int * cycleLengths;
            /**< The length of each cycle as stored in the array above. */
        unsigned int nAutos;
            /**< The number of automorphisms of the graph. */
        Automorphism** automorphisms;
            /**< Store the automorphisms as maps from the edges in the
             *   graph to edges in the automorphism of the graph. */
        unsigned int * isAutoOk;
            /**< An array of size nAutos. Each entry corresponds an
             *   automorphism, such if automorphism i is not valid when 
             *   adding edge j, isAutoOk[i] = j.
             *   The value of j here is purely a number indicating how many
             *   edges have been added so far. */
        bool orientable;
            /**< Whether we are searching for orientable manifolds. */
        bool minimal_;
            /**< Is the search only looking for minimal triangulations. */

        void colourLowestEdge();
            /**< Start colouring a cycle on the lowest edge. */
       
        bool checkColourOk();
            /**< Checks whether the cycle indicated by colour nextColour
             *   is a valid cycle. */

        bool isCanonical();
            /**< Determine whether the current internal edge allocation, as
             *   The parameter otherOption, if non-zero, */

        void  nextPath(EdgeEnd *start, unsigned int firstEdge, EdgeEnd *now);
            /**< Tries all choices for continuing a cycle from *now.
             *   Also checks to see if a cycle can be completed after using
             *   one of these choices. */

        bool checkComplete();
            /**< Checks to see if a decomposition has been found. This is done
             *   by checking to see that each Edge has 3 cycles on it. */
  
        unsigned int findTetWithMostInternalEdgesUsed();
            /**< Finds the tetrahedra with the most used internal edges. */

        unsigned int compareCycles(unsigned int *cycleListA, unsigned int *cycleListB, 
                unsigned int lengthA, unsigned int lengthB,
                unsigned int offsetA, unsigned int offsetB); 
            /**< Check to see which of two cycles is more canonical.  Return values are
             *   0 if A < B
             *   1 if B > A
             *   2 if A == B
             *   */
    public:
        static const char dataTag_;
            /**< A character used to identify this class when reading
                 and writing tagged data in text format. */


    public:
        /**
         * Creates a new search manager for use when (i) only closed prime
         * minimal P2-irreducible triangulations are required, and (ii) the
         * given face pairing has order at least three.  Note that other
         * unwanted triangulations may still be produced (e.g.,
         * non-prime or non-minimal triangulations), but there will be
         * far fewer of these than when using the NGluingPermSearcher
         * class directly.
         *
         * For details on how a search manager is used, see the
         * NGluingPermSearcher documentation.  Note in particular that
         * this class will be automatically used by
         * NGluingPermSearcher::findAllPerms() if possible, so there is
         * often no need for an end user to instantiate this class
         * directly.
         *
         * All constructor arguments are the same as for the
         * NGluingPermSearcher constructor, though some arguments (such as
         * \a finiteOnly and \a whichPurge) are not needed here since they
         * are already implied by the specialised search context.
         *
         * \pre The given face pairing is connected, i.e., it is possible
         * to reach any tetrahedron from any other tetrahedron via a
         * series of matched face pairs.
         * \pre The given face pairing is in canonical form as described
         * by NFacePairing::isCanonical().  Note that all face pairings
         * constructed by NFacePairing::findAllPairings() are of this form.
         * \pre The given face pairing has no boundary faces and has at
         * least three tetrahedra.
         */
        CycleDecompSearcher(const regina::NFacePairing* pairing,
                const regina::NFacePairing::IsoList* autos,
                bool orientableOnly, UseCycles use, void* useArgs = 0,
                bool minimal = true);


        /**
         * Create a triangulation from the data the search manager has 
         * found. */
        regina::NTriangulation* triangulate() const;

        /**
         * Initialises a new search manager based on data read from the
         * given input stream.  This may be a new search or a partially
         * completed search.
         *
         * This routine reads data in the format written by dumpData().
         * If you wish to read data whose precise class is unknown,
         * consider using dumpTaggedData() and readTaggedData() instead.
         *
         * If the data found in the input stream is invalid or incorrectly
         * formatted, the routine inputError() will return \c true but
         * the contents of this object will be otherwise undefined.
         *
         * The arguments \a use and \a useArgs are the same as for the
         * NGluingPermSearcher constructor.
         *
         * \warning The data format is liable to change between Regina
         * releases.  Data in this format should be used on a short-term
         * temporary basis only.
         *
         * @param in the input stream from which to read.
         */
        //NCycleDecompSearcher(std::istream& in,
        //    UseGluingPerms use, void* useArgs = 0);

        /**
         * Destroys this search manager and all supporting data
         * structures.
         */
        ~CycleDecompSearcher();

        // Overridden methods:
        virtual void dumpData(std::ostream& out) const;
        virtual void runSearch(long maxDepth = -1);


    protected:
        // Overridden methods:
        virtual char dataTag() const;

};

inline char CycleDecompSearcher::dataTag() const {
    return CycleDecompSearcher::dataTag_;
}

inline CycleDecompSearcher::Edge::Edge() {
    used=0;
}

inline CycleDecompSearcher::EdgeEnd::EdgeEnd() {
    for(unsigned int i=0; i<6;i++) 
        map[i]=0;
}

inline CycleDecompSearcher::Tetrahedron::Tetrahedron() {
    used=0;
    for(unsigned int i=0; i<6;i++) 
        internalEdges[i]=0;
}

#endif

/* vim: set ts=4 sw=4 expandtab: */
