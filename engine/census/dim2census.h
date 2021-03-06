
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
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

/*! \file census/dim2census.h
 *  \brief Deals with forming a census of all 2-manifold triangulations
 *  of a given size.
 */

#ifndef __DIM2CENSUS_H
#ifndef __DOXYGEN
#define __DIM2CENSUS_H
#endif

#include "regina-core.h"
#include "dim2/dim2edgepairing.h"
#include "utilities/nbooleans.h"

namespace regina {

class Dim2GluingPerms;
class Dim2GluingPermSearcher;
class NPacket;

template <int> class Triangulation;
typedef Triangulation<2> Dim2Triangulation;

/**
 * \weakgroup census
 * @{
 */

/**
 * A utility class used to search for triangulations across one or more
 * 2-manifold census databases.
 */
class REGINA_API Dim2Census {
    public:
        /**
         * A routine used to determine whether a particular 2-manifold
         * triangulation should be included in a census.  Routines of this
         * type are used by Dim2Census::formCensus().
         *
         * The first parameter passed should be a triangulation currently
         * under consideration.
         * The second parameter may contain arbitrary data as passed to
         * Dim2Census::formCensus().
         *
         * The return value should be \c true if the triangulation passed
         * should be included in the census, or \c false otherwise.
         *
         * \deprecated The Dim2Census enumeration facilities are on their
         * way out of Regina, and in the future the Dim2Census class will be
         * used purely for census lookups.  If you wish to build a census
         * yourself, you should call Dim2EdgePairing::findAllPairings() and
         * Dim2GluingPermSearcher::findAllPerms() directly.
         */
        REGINA_DEPRECATED typedef bool (*AcceptTriangulation)(
            Dim2Triangulation*, void*);

    private:
        NPacket* parent_;
            /**< The argument passed to formCensus(). */
        NBoolSet orientability_;
            /**< The argument passed to formCensus(). */

        bool (*sieve_)(Dim2Triangulation*, void*);
            /**< The arbitrary constraint function to run triangulations
                 through.  This is of type AcceptTriangulation, but because
                 that typedef is deprecated we do not use it here. */
        void* sieveArgs_;
            /**< The second argument to pass to function \a sieve. */

        unsigned long whichSoln_;
            /**< The number of the solution we are up to. */

    public:
        /**
         * Deprecated routine that fills the given packet with all
         * triangulations in a census of 2-manifold triangulations satisfying
         * the given constraints.  Each triangulation in the census will
         * appear as a child of the given packet.
         *
         * This routine will conduct a census of all valid triangulations
         * containing a given number of triangles.  All such triangulations
         * are included in the census up to combinatorial isomorphism; given
         * any isomorphism class, exactly one representative will appear in
         * the census.
         *
         * The census can be optionally restricted to only include
         * triangulations satisfying further constraints (such as
         * orientability); see the individual parameter
         * descriptions for further details.  In particular, parameter
         * \a sieve can be used to impose arbitrary restrictions that are
         * not hard-coded into this class.
         *
         * Note that if constraints may be imposed using the hard-coded
         * parameters (such as orientability), it is
         * generally better to do this than to use the arbitrary
         * constraint parameter \a sieve.  Hard-coded parameters will be
         * tested earlier, and some (such as orientability) can be
         * incorporated directly into the census algorithm to give a vast
         * performance increase.
         *
         * Note that this routine should only be used if the census
         * contains a small enough total number of triangulations to
         * avoid any memory disasters.
         *
         * \deprecated The Dim2Census enumeration facilities are on their
         * way out of Regina, and in the future the Dim2Census class will be
         * used purely for census lookups.  To perform the kind of
         * enumeration that is described here, you should call
         * Dim2EdgePairing::findAllPairings() and
         * Dim2GluingPermSearcher::findAllPerms() directly.
         *
         * \ifacespython Parameters \a sieve and \a sieveArgs
         * are not present (and will be treated as 0).
         *
         * @param parent the packet beneath which members of the census will
         * be placed.
         * @param nTriangles the number of triangles in each triangulation
         * in the census.
         * @param orientability determines whether to include orientable
         * and/or non-orientable triangulations.  The set should contain \c true
         * if orientable triangulations are to be included, and should contain
         * \c false if non-orientable triangulations are to be included.
         * @param boundary determines whether to include triangulations with
         * and/or without boundary edges.  The set should contain \c true
         * if triangulations with boundary edges are to be included, and
         * should contain \c false if triangulations with only internal
         * edges are to be included.
         * @param nBdryEdges specifies the precise number of boundary edges
         * that should be present in the triangulations produced.  If this
         * parameter is negative, it is ignored and no additional restriction
         * is imposed.  See the documentation for routine
         * Dim2EdgePairing::findAllPairings() for details regarding this
         * parameter and how it interacts with parameter \a boundary.
         * @param sieve an additional constraint function that may be
         * used to exclude certain triangulations from the census.  If
         * this parameter is non-zero, each triangulation produced (after
         * passing all other criteria) will be passed through this
         * function.  If this function returns \c true then the triangulation
         * will be included in the census; otherwise it will not.
         * When this function is called, the first (triangulation)
         * argument will be a triangulation under consideration for
         * inclusion in the census.  The second argument will be
         * parameter \a sieveArgs as passed to formCensus().
         * Parameter \a sieve may be passed as \c null (in which case no
         * additional constraint function will be used).
         * @param sieveArgs the pointer to pass as the final parameter
         * for the function \a sieve which will be called upon each
         * triangulation found.  If \a sieve is \c null then \a sieveArgs
         * will be ignored.
         * @return the number of triangulations produced in the census.
         */
        REGINA_DEPRECATED static unsigned long formCensus(NPacket* parent,
            unsigned nTriangles, NBoolSet orientability, NBoolSet boundary,
            int nBdryEdges, AcceptTriangulation sieve = 0,
            void* sieveArgs = 0);

        /**
         * Deprecated routine that fills the given packet with all
         * triangulations in a partial census of 2-manifold triangulations
         * satisfying the given constraints.  Each triangulation in the
         * partial census will appear as a child of the given packet.
         *
         * This routine will conduct a census of all valid 2-manifold
         * triangulations that are modelled by the given triangle edge
         * pairing.  All such triangulations are included in the census up to
         * combinatorial isomorphism; given any isomorphism class, exactly
         * one representative will appear in the census.
         *
         * The census can be optionally restricted to only include
         * triangulations satisfying further constraints (such as
         * orientability); see the individual parameter
         * descriptions for further details.  In particular, parameter
         * \a sieve can be used to impose arbitrary restrictions that are
         * not hard-coded into this class.
         *
         * Note that if constraints may be imposed using the hard-coded
         * parameters (such as orientability), it is
         * generally better to do this than to use the arbitrary
         * constraint parameter \a sieve.  Hard-coded parameters will be
         * tested earlier, and some (such as orientability) can be
         * incorporated directly into the census algorithm to give a vast
         * performance increase.
         *
         * Note that this routine should only be used if the partial census
         * contains a small enough total number of triangulations to
         * avoid any memory disasters.
         *
         * The partial census will run in the current thread.  This
         * routine will only return once the partial census is complete.
         *
         * \deprecated The Dim2Census enumeration facilities are on their
         * way out of Regina, and in the future the Dim2Census class will be
         * used purely for census lookups (if at all).  To perform the kind of
         * enumeration that is described here, you should call
         * Dim2GluingPermSearcher::findAllPerms() directly.
         *
         * \pre The given edge pairing is connected, i.e., it is possible
         * to reach any triangle from any other triangle via a
         * series of matched edge pairs.
         * \pre The given edge pairing is in canonical form as described
         * by Dim2EdgePairing::isCanonical().  Note that all facet pairings
         * constructed by Dim2EdgePairing::findAllPairings() are of this form.
         *
         * \ifacespython Parameters \a sieve and \a sieveArgs
         * are not present (and will be treated as 0).
         *
         * @param pairing the edge pairing that
         * triangulations in this partial census must be modelled by.
         * @param parent the packet beneath which members of the partial
         * census will be placed.
         * @param orientability determines whether to include orientable
         * and/or non-orientable triangulations.  The set should contain \c true
         * if orientable triangulations are to be included, and should contain
         * \c false if non-orientable triangulations are to be included.
         * @param sieve an additional constraint function that may be
         * used to exclude certain triangulations from the census.  If
         * this parameter is non-zero, each triangulation produced (after
         * passing all other criteria) will be passed through this
         * function.  If this function returns \c true then the triangulation
         * will be included in the census; otherwise it will not.
         * When this function is called, the first (triangulation)
         * argument will be a triangulation under consideration for
         * inclusion in the census.  The second argument will be
         * parameter \a sieveArgs as passed to formPartialCensus().
         * Parameter \a sieve may be passed as \c null (in which case no
         * additional constraint function will be used).
         * @param sieveArgs the pointer to pass as the final parameter
         * for the function \a sieve which will be called upon each
         * triangulation found.  If \a sieve is \c null then \a sieveArgs
         * will be ignored.
         * @return the number of triangulations produced in the partial census.
         */
        REGINA_DEPRECATED static unsigned long formPartialCensus(
            const Dim2EdgePairing* pairing, NPacket* parent,
            NBoolSet orientability, AcceptTriangulation sieve = 0,
            void* sieveArgs = 0);

    private:
        /**
         * Creates a new structure to hold the given information.
         * All parameters not explained are taken directly from
         * formCensus().
         *
         * \deprecated The Dim2Census enumeration facilities are on their
         * way out of Regina, and in the future the Dim2Census class will be
         * used purely for census lookups.  If you wish to build a census
         * yourself, you should call Dim2EdgePairing::findAllPairings() and
         * Dim2GluingPermSearcher::findAllPerms() directly.
         */
        REGINA_DEPRECATED Dim2Census(NPacket* parent,
            const NBoolSet& orientability,
            AcceptTriangulation sieve, void* sieveArgs);

        /**
         * Called when a particular triangle edge pairing has been
         * found.  This routine hooks up the edge pairing generation with
         * the gluing permutation generation.
         *
         * @param pairing the edge pairing that has been found.
         * @param autos the set of all automorphisms of the given edge
         * pairing.
         * @param census the census currently being generated;
         * this must really be of class Dim2Census.
         */
        static void foundEdgePairing(const Dim2EdgePairing* pairing,
            const Dim2EdgePairing::IsoList* autos, void* census);

        /**
         * Called when a particular set of gluing permutations has been
         * found.  This routine generates the corresponding triangulation
         * and decides whether it really belongs in the census.
         *
         * \pre The given set of gluing permutations is complete, i.e.,
         * it is not just a partially-complete search state.  Partial
         * searches are formed by calling Dim2GluingPermSearcher::runSearch()
         * with a non-negative depth argument.
         *
         * @param perms the gluing permutation set that has been found.
         * @param census the census currently being generated;
         * this must really be of class Dim2Census.
         */
        static void foundGluingPerms(const Dim2GluingPermSearcher* perms,
            void* census);
};

/*@}*/

} // namespace regina

#endif

