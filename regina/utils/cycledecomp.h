
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

#ifndef __NGLUINGPERMSEARCHER_H
#ifndef __DOXYGEN
#define __NGLUINGPERMSEARCHER_H
#endif

#include "regina-core.h"
#include "census/ngluingperms.h"
#include "utilities/nqitmask.h"

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
typedef void (*UseGluingPerms)(const NGluingPermSearcher*, void*);


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

        class Edge {
            /**< Represents a single edge of the face pairingr graph. */
            public:
                signed colours[3];
                    /**< The 3 cycles an edge may have. */
                int used;
                    /**< The number of cycles currently using this edge. */
                unsigned endTet[2];
                    /**< The tetrahedron on either end  of this edge. */
                EdgeEnd* ends[2];
                    /**< The edge-ends of this edge. */
        };

        class EdgeEnd {
            /**< The "end" of an edge in the face pairing graph. This end
             *   contains information about how the edge is attached to the
             *   tetrahedron. */
            public:
                Tetrahedron* tet;
                    /**< The tetrahedron attached here. */
                Edge *edge;
                    /**< The edge containing this end. */
                unsigned face;
                    /**< The face of the tetrahedron which this edge is
                     *   attached to. */
                signed map[6];
                    /**< A mapping of the edges of the face to the cycles on
                     *   the edge. */

        };

        class Tetrahedron {
            /**< A tetrahedron, as represented in the face pairing graph. */
            public:
                int used;
                    /**< The number of internal edges used. */
                unsigned internalEdges[6];
                    /**< The 6 internal edges of the tetrahedron. */
                EdgeEnd* externalEdgeEnd[4];
                    /**< The EdgeEnds which are attached to each face of the
                     *   tetrahedron. */

        };

        Tetrahedron[] tets;
            /**< The tetrahedron representations in the face pairing graph. */
        unsigned nTets;
            /**< The number of tetrahedra. */
        Edge[] edges;
            /**< The edges of the face pairing graph. */
        unsigned nEdges;
            /**< The number of edges. */
        EdgeEnd[] ends;
            /**< The ends of the edges. */
        unsigned nEnds;
            /**< The number of edge ends. */
        unsigned nextColour;
            /**< The next "colour" to be used to mark out a cycle in the 
             *   face pairing graph. */
        unsigned** cycles;
            /**< Contains a list of all edge numbers used in each cycle.
             *   cycles[x][y] denotes the y'th edge in cycle number x. 
             *   These are the spine codes as Matveev uses them. */
        unsigned* cycleLengths;
            /**< The length of each cycle as stored in the array above. */

        

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
        NCycleDecompSearcher(const NFacePairing* pairing,
                const NFacePairing::IsoList* autos,
                bool orientableOnly, UseGluingPerms use, void* useArgs = 0);

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
        NCycleDecompSearcher(std::istream& in,
            UseGluingPerms use, void* useArgs = 0);

        /**
         * Destroys this search manager and all supporting data
         * structures.
         */
        virtual ~NCycleDecompSearcher();

        // Overridden methods:
        virtual void dumpData(std::ostream& out) const;
        virtual void runSearch(long maxDepth = -1);

    protected:
        // Overridden methods:
        virtual char dataTag() const;

};

inline char NCycleDecompSearcher::dataTag() const {
    return NCycleDecompSearcher::dataTag_;
}

inline NCycleDecompSearcher::Edge::Edge() {
    used=0;
}

inline NCycleDecompSearcher::Tetrahedron::Tetrahedron() {
    used=0;
    for(unsigned i=0; i<6;i++) 
        internalEdges[i]=0;
}

#endif

