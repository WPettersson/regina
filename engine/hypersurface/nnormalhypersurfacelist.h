
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

/*! \file hypersurface/nnormalhypersurfacelist.h
 *  \brief Contains a packet representing a collection of normal
 *  hypersurfaces in a 4-manifold triangulation.
 */

#ifndef __NNORMALHYPERSURFACELIST_H
#ifndef __DOXYGEN
#define __NNORMALHYPERSURFACELIST_H
#endif

#include <algorithm>
#include <iterator>
#include <vector>
#include "regina-core.h"
#include "hypersurface/hypercoords.h"
#include "hypersurface/hyperflags.h"
#include "hypersurface/nnormalhypersurface.h"
#include "packet/npacket.h"
#include "utilities/memutils.h"

namespace regina {

class NMatrixInt;
class NNormalHypersurface;
class NNormalHypersurfaceList;
class NProgressTracker;
class NXMLNormalHypersurfaceListReader;
class NXMLPacketReader;

template <int> class Triangulation;
// We *must* declare the specialisation here; otherwise this header has
// the side-effect of instantiating Trianglation<4> using the generic template.
template <> class Triangulation<4>;
typedef Triangulation<4> Dim4Triangulation;

/**
 * \weakgroup hypersurface
 * @{
 */

#ifndef __DOXYGEN // Doxygen complains about undocumented specialisations.
template <>
struct PacketInfo<PACKET_NORMALHYPERSURFACELIST> {
    typedef NNormalHypersurfaceList Class;
    inline static const char* name() {
        return "Normal Hypersurface List";
    }
};
#endif

/**
 * A packet representing a collection of normal hypersurfaces in a 4-manifold
 * triangulation.  Such a packet must always be a child packet of the
 * triangulation from which the surfaces were obtained.  If this triangulation
 * changes, the information contained in this packet will become invalid.
 *
 * See the NNormalHypersurfaceVector class notes for details of what to do
 * when introducing a new coordinate system.
 *
 * Normal hypersurface lists should be created using the routine enumerate().
 */
class REGINA_API NNormalHypersurfaceList : public NPacket {
    REGINA_PACKET(NNormalHypersurfaceList, PACKET_NORMALHYPERSURFACELIST)

    public:
        class VectorIterator;

    protected:
        std::vector<NNormalHypersurface*> surfaces_;
            /**< Contains the normal hypersurfaces stored in this packet. */
        HyperCoords coords_;
            /**< Stores which coordinate system is being
                 used by the normal hypersurfaces in this packet. */
        HyperList which_;
            /**< Indicates which normal hypersurfaces these represent
                 within the underlying triangulation. */
        HyperAlg algorithm_;
            /**< Stores the details of the enumeration algorithm that
                 was used to generate this list.  This might not be the
                 same as the \a algorithmHints flag passed to the
                 corresponding enumeration routine (e.g., if invalid or
                 inappropriate flags were passed). */

    public:
        /**
         * Destroys this list and all the hypersurfaces within.
         */
        virtual ~NNormalHypersurfaceList();

        /**
         * Enumerates all vertex normal hypersurfaces in the given
         * triangulation using the given coordinate system.
         * These vertex normal hypersurfaces will be stored in a new normal
         * hypersurface list.  Their representations will
         * use the smallest possible integer coordinates.
         * The option is offered to find only embedded normal hypersurfaces
         * or to also include immersed and singular normal hypersurfaces.
         *
         * The normal hypersurface list that is created will be inserted as the
         * last child of the given triangulation.  This triangulation \b must
         * remain the parent of this normal hypersurface list, and must not
         * change while this normal hypersurface list remains in existence.
         *
         * If a progress tracker is passed, the normal hypersurface
         * enumeration will take place in a new thread and this routine
         * will return immediately.  If the user cancels the operation
         * from another thread, then the normal hypersurface list will \e not
         * be inserted into the packet tree (but the caller of this
         * routine will still need to delete it).  Regarding progress
         * tracking, this routine will declare and work through a series
         * of stages whose combined weights sum to 1; typically this
         * means that the given tracker must not have been used before.
         *
         * If no progress tracker is passed, the enumeration will run
         * in the current thread and this routine will return only when
         * the enumeration is complete.  Note that this enumeration can
         * be extremely slow for larger triangulations.
         *
         * @param owner the triangulation upon which this list of normal
         * hypersurfaces will be based.
         * @param coords the coordinate system to be used.
         * @param embeddedOnly \c true if only embedded normal hypersurfaces
         * are to be produced, or \c false if immersed and singular
         * normal hypersurfaces are also to be produced; this defaults to
         * \c true.
         * @param tracker a progress tracker through which progress will
         * be reported, or 0 if no progress reporting is required.
         * @return the newly created normal hypersurface list.  Note that if
         * a progress tracker is passed then this list may not be completely
         * filled when this routine returns.  If a progress tracker is
         * passed and a new thread could not be started, this routine
         * returns 0 (and no normal hypersurface list is created).
         */
        static NNormalHypersurfaceList* enumerate(Dim4Triangulation* owner,
            HyperCoords coords,
            HyperList which = HS_LIST_DEFAULT,
            HyperAlg algHints = HS_ALG_DEFAULT,
            NProgressTracker* tracker = 0); // TODO

        /**
         * Returns the coordinate system being used by the
         * hypersurfaces stored in this set.
         *
         * @return the coordinate system used.
         */
        HyperCoords coords() const;
        /**
         * Returns details of which normal hypersurfaces this list represents
         * within the underlying triangulation.
         *
         * This may not be the same HyperList that was passed to enumerate().
         * In particular, default values will have been explicitly
         * filled in (such as HS_VERTEX and/or HS_EMBEDDED_ONLY), and
         * invalid and/or redundant values will have been removed.
         *
         * @return details of what this list represents.
         */
        HyperList which() const;
        /**
         * Returns details of the algorithm that was used to enumerate
         * this list.
         *
         * These may not be the same HyperAlg flags that were passed to
         * enumerate().  In particular, default values will have been
         * explicitly filled in, invalid and/or redundant values will have
         * been removed, and unavailable and/or unsupported combinations
         * of algorithm flags will be replaced with whatever algorithm
         * was actually used.
         *
         * @return details of the algorithm used to enumerate this list.
         */
        HyperAlg algorithm() const;
        /**
         * Returns whether this set is known to contain only embedded normal
         * hypersurfaces.
         *
         * If this returns \c false, it does not guarantee that immersed
         * and/or singular hypersurfaces are present; it merely indicates
         * that they were not deliberately excluded (for instance, the
         * prism constraints were not enforced).
         *
         * @return \c true if this list was constructed to contain only
         * properly embedded hypersurfaces, or \c false otherwise.
         */
        bool isEmbeddedOnly() const;
        /**
         * Returns the triangulation in which these normal hypersurfaces live.
         *
         * @return the triangulation in which these hypersurfaces live.
         */
        Dim4Triangulation* triangulation() const;

        /**
         * Returns the number of hypersurfaces stored in this list.
         *
         * @return the number of hypersurfaces.
         */
        size_t size() const;
        /**
         * Returns the hypersurface at the requested index in this list.
         *
         * @param index the index of the requested hypersurface in this list;
         * this must be between 0 and size()-1 inclusive.
         *
         * @return the normal hypersurface at the requested index in this list.
         */
        const NNormalHypersurface* hypersurface(size_t index) const;

        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeTextLong(std::ostream& out) const;
        static NXMLPacketReader* xmlReader(NPacket* parent,
            NXMLTreeResolver& resolver);
        virtual bool dependsOnParent() const;

        /**
         * Returns a newly created matrix containing the matching
         * equations that were used to create this normal hypersurface list.
         * The destruction of this matrix is the responsibility of the
         * caller of this routine.  Multiple calls to this routine will
         * result in the construction of multiple matrices.  This
         * routine in fact merely calls makeMatchingEquations() with the
         * appropriate parameters.
         *
         * The format of the matrix is identical to that returned by
         * makeMatchingEquations().
         * 
         * @return the matching equations used to create this normal
         * hypersurface list.
         */
        NMatrixInt* recreateMatchingEquations() const;

        /**
         * An iterator that gives access to the raw vectors for hypersurfaces
         * in this list, pointing to the beginning of this hypersurface list.
         *
         * \ifacespython Not present.
         *
         * @return an iterator at the beginning of this hypersurface list.
         */
        VectorIterator beginVectors() const;

        /**
         * An iterator that gives access to the raw vectors for hypersurfaces
         * in this list, pointing past the end of this hypersurface list.
         * This iterator is not dereferenceable.
         *
         * \ifacespython Not present.
         *
         * @return an iterator past the end of this hypersurface list.
         */
        VectorIterator endVectors() const;

        /**
         * A bidirectional iterator that runs through the raw vectors for
         * hypersurfaces in this list.
         *
         * \ifacespython Not present.
         */
        class VectorIterator : public std::iterator<
                std::bidirectional_iterator_tag,
                const NNormalHypersurfaceVector*> {
            private:
                std::vector<NNormalHypersurface*>::const_iterator it_;
                    /**< An iterator into the underlying list of
                         hypersurfaces. */

            public:
                /**
                 * Creates a new uninitialised iterator.
                 */
                VectorIterator();

                /**
                 * Creates a copy of the given iterator.
                 *
                 * @param cloneMe the iterator to clone.
                 */
                VectorIterator(const VectorIterator& cloneMe);

                /**
                 * Makes this a copy of the given iterator.
                 *
                 * @param cloneMe the iterator to clone.
                 * @return a reference to this iterator.
                 */
                VectorIterator& operator = (const VectorIterator& cloneMe);

                /**
                 * Compares this with the given operator for equality.
                 *
                 * @param other the iterator to compare this with.
                 * @return \c true if the iterators point to the same
                 * element of the same normal surface list, or \c false
                 * if they do not.
                 */
                bool operator == (const VectorIterator& other) const;

                /**
                 * Compares this with the given operator for inequality.
                 *
                 * @param other the iterator to compare this with.
                 * @return \c false if the iterators point to the same
                 * element of the same normal surface list, or \c true
                 * if they do not.
                 */
                bool operator != (const VectorIterator& other) const;

                /**
                 * Returns the raw vector for the normal hypersurface that this
                 * iterator is currently pointing to.
                 *
                 * \pre This iterator is dereferenceable (in particular,
                 * it is not past-the-end).
                 *
                 * @return the corresponding normal hypersurface vector.
                 */
                const NNormalHypersurfaceVector* operator *() const;

                /**
                 * The preincrement operator.
                 *
                 * @return a reference to this iterator after the increment.
                 */
                VectorIterator& operator ++();

                /**
                 * The postincrement operator.
                 *
                 * @return a copy of this iterator before the
                 * increment took place.
                 */
                VectorIterator operator ++(int);

                /**
                 * The predecrement operator.
                 *
                 * @return a reference to this iterator after the decrement.
                 */
                VectorIterator& operator --();

                /**
                 * The postdecrement operator.
                 *
                 * @return a copy of this iterator before the
                 * decrement took place.
                 */
                VectorIterator operator --(int);

            private:
                /**
                 * Initialise a new vector iterator using an iterator for
                 * the internal list of normal hypersurfaces.
                 */
                VectorIterator(
                    const std::vector<NNormalHypersurface*>::const_iterator& i);

            friend class NNormalHypersurfaceList;
        };

    protected:
        /**
         * Creates an empty list of normal hypersurfaces with the given
         * parameters.
         *
         * @param coords the coordinate system to be used for filling
         * this list.
         * @param which indicates which normal hypersurfaces these will
         * represent within the underlying triangulation.
         * @param algorithm details of the enumeration algorithm that
         * will be used to fill this list.
         */
        NNormalHypersurfaceList(HyperCoords coords, HyperList which,
            HyperAlg algorithm);

        virtual NPacket* internalClonePacket(NPacket* parent) const;
        virtual void writeXMLPacketData(std::ostream& out) const;

        /**
         * An output iterator used to insert hypersurfaces into an
         * NNormalHypersurfaceList.
         *
         * Objects of type <tt>NNormalHypersurface*</tt> and
         * <tt>NNormalHypersurfaceVector*</tt> can be assigned to this
         * iterator.  In the latter case, a surrounding NNormalHypersurface
         * will be automatically created.
         */
        struct HypersurfaceInserter : public std::iterator<
                std::output_iterator_tag, NNormalHypersurfaceVector*> {
            NNormalHypersurfaceList* list_;
                /**< The list into which hypersurfaces will be inserted. */
            Dim4Triangulation* owner_;
                /**< The triangulation in which the hypersurfaces to be
                 *   inserted are contained. */

            /**
             * Creates a new output iterator.  The member variables of
             * this iterator will be initialised according to the
             * parameters passed to this constructor.
             *
             * @param list the list into which hypersurfaces will be inserted.
             * @param owner the triangulation in which the hypersurfaces
             * to be inserted are contained.
             */
            HypersurfaceInserter(NNormalHypersurfaceList& list,
                Dim4Triangulation* owner);
            /**
             * Creates a new output iterator that is a clone of the
             * given iterator.
             *
             * @param cloneMe the output iterator to clone.
             */
            HypersurfaceInserter(const HypersurfaceInserter& cloneMe);

            /**
             * Sets this iterator to be a clone of the given output iterator.
             *
             * @param cloneMe the output iterator to clone.
             * @return this output iterator.
             */
            HypersurfaceInserter& operator = (
                const HypersurfaceInserter& cloneMe);

            /**
             * Appends a normal hypersurface to the end of the appropriate
             * surface list.
             *
             * The given hypersurface will be deallocated with the other
             * hypersurfaces in this list when the list is eventually
             * destroyed.
             *
             * @param hypersurface the normal hypersurface to insert.
             * @return this output iterator.
             */
            HypersurfaceInserter& operator = (
                NNormalHypersurface* hypersurface);
            /**
             * Appends the normal hypersurface corresponding to the given
             * vector to the end of the appropriate hypersurface list.
             *
             * The given vector will be owned by the newly created
             * normal hypersurface and will be deallocated with the other
             * hypersurfaces in this list when the list is eventually
             * destroyed.
             *
             * @param vector the vector of the normal hypersurface to insert.
             * @return this output iterator.
             */
            HypersurfaceInserter& operator = (
                NNormalHypersurfaceVector* vector);

            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            HypersurfaceInserter& operator *();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            HypersurfaceInserter& operator ++();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            HypersurfaceInserter& operator ++(int);
        };

    private:
        /**
         * A functor that performs all normal hypersurface enumeration.
         */
        class Enumerator {
            private:
                NNormalHypersurfaceList* list_;
                    /**< The hypersurface list to be filled. */
                Dim4Triangulation* triang_;
                    /**< The triangulation in which these hypersurfaces lie. */
                NProgressTracker* tracker_;
                    /**< The progress tracker through which progress is
                         reported and cancellation requests are accepted,
                         or 0 if no progress tracker is in use. */

            public:
                /**
                 * Creates a new functor with the given parameters.
                 *
                 * @param list the hypersurface list to be filled.
                 * @param triang the triangulation in which these
                 * hypersurfaces lie.
                 * @param tracker the progress tracker to use for
                 * progress reporting and cancellation polling, or 0 if
                 * these capabilities are not required.
                 */
                Enumerator(NNormalHypersurfaceList* list,
                    Dim4Triangulation* triang, NProgressTracker* tracker);

                /**
                 * Performs the real enumeration work, in a setting
                 * where the underlying coordinate system is
                 * a compile-time constant.
                 *
                 * We assume here that neither list_->which_ nor
                 * list_->algorithm_ have been sanity-checked.
                 *
                 * This routine fills \a list_ with surfaces, and then once
                 * this is finished it inserts \a list_ into the packet
                 * tree as a child of \a triang_.
                 *
                 * \tparam Coords an instance of the HyperInfo<> template class.
                 */
                template <typename Coords>
                void operator() ();

            private:
                /**
                 * The enumeration code for enumerating vertex hypersurfaces.
                 * This is internal to operator().
                 *
                 * We assume that the flag set which_ is set correctly,
                 * and we do not alter it here.
                 * We make no assumptions about the state of algorithm_,
                 * and we set this during the course of this routine.
                 *
                 * This routine only fills \a list_ with hypersurfaces.
                 * It does not make any adjustments to the structure of
                 * the packet tree.
                 *
                 * If \a tracker_ is non-null, this routine will declare and
                 * work through a series of tracker stages whose
                 * combined weights sum to 1.  It will not, however,
                 * call NProgressTracker::setFinished().
                 */
                template <typename Coords>
                void fillVertex();

                /**
                 * The enumeration code for enumerating fundamental
                 * hypersurfaces.  This is internal to operator().
                 *
                 * We assume that the flag set which_ is set correctly,
                 * and we do not alter it here.
                 * We make no assumptions about the state of algorithm_,
                 * and we set this during the course of this routine.
                 *
                 * This routine only fills \a list_ with surfaces.
                 * It does not make any adjustments to the structure of
                 * the packet tree.
                 *
                 * If \a tracker_ is non-null, this routine declare and
                 * work through a series of tracker stages whose
                 * combined weights sum to 1.  It will not, however,
                 * call NProgressTracker::setFinished().
                 */
                template <typename Coords>
                void fillFundamental();

                /**
                 * The enumeration code for enumerating vertex surfaces
                 * using the double description method.
                 * This is internal to fillVertex().
                 *
                 * This routine assumes that \a algorithm_ has been set
                 * correctly, and does not alter it.
                 *
                 * If \a tracker_ is non-null, this routine assumes that
                 * an appropriate tracker stage has already been
                 * declared, and works through that stage only.
                 *
                 * \pre The underlying triangulation is non-empty.
                 */
                template <typename Coords>
                void fillVertexDD();

                /**
                 * The enumeration code for enumerating fundamental surfaces
                 * using the primal method.
                 * This is internal to fillFundamental().
                 *
                 * This routine assumes nothing about the state of the
                 * \a algorithm_ flag set, and sets it appropriately.
                 *
                 * If \a tracker_ is non-null, this routine will declare and
                 * work through a series of tracker stages whose
                 * combined weights sum to 1.  It will not, however,
                 * call NProgressTracker::setFinished().
                 *
                 * \pre The underlying triangulation is non-empty.
                 */
                template <typename Coords>
                void fillFundamentalPrimal();

                /**
                 * The enumeration code for enumerating fundamental surfaces
                 * using the dual method.
                 * This is internal to fillFundamental().
                 *
                 * This routine assumes nothing about the state of the
                 * \a algorithm_ flag set, and sets it appropriately.
                 *
                 * If \a tracker_ is non-null, this routine will declare and
                 * work through a series of tracker stages whose
                 * combined weights sum to 1.  It will not, however,
                 * call NProgressTracker::setFinished().
                 *
                 * \pre The underlying triangulation is non-empty.
                 */
                template <typename Coords>
                void fillFundamentalDual();
        };

    friend class regina::NXMLNormalHypersurfaceListReader;
};

/**
 * Returns a new normal hypersurface vector of the appropriate length for the
 * given triangulation and the given coordinate system.
 * All elements of this vector will be initialised to zero.
 *
 * The new vector will be of the subclass of NNormalHypersurfaceVector
 * corresponding to the given coordinate system.  The caller
 * of this routine is responsible for destroying the new vector.
 *
 * \ifacespython Not present.
 *
 * @param triangulation the triangulation upon which the underlying
 * coordinate system is based.
 * @param coords the coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalHypersurfaceList.
 * @return a new zero vector of the correct class and length.
 */
REGINA_API NNormalHypersurfaceVector* makeZeroVector(
    const Dim4Triangulation* triangulation, HyperCoords coords);
/**
 * Creates a new set of normal hypersurface matching equations for the
 * given triangulation using the given coordinate system.
 * The returned matrix will be newly allocated and its destruction will
 * be the responsibility of the caller of this routine.
 *
 * Each equation will be represented as a row of the matrix.
 * Each column of the matrix represents a coordinate in the given
 * coordinate system.
 *
 * @param triangulation the triangulation upon which these matching equations
 * will be based.
 * @param coords the coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalHypersurfaceList.
 * @return a newly allocated set of matching equations.
 */
REGINA_API NMatrixInt* makeMatchingEquations(
    const Dim4Triangulation* triangulation, HyperCoords coords);
/**
 * Creates a new set of validity constraints representing the condition that
 * normal hypersurfaces be embedded.  The validity constraints will be expressed
 * relative to the given coordinate system.
 *
 * \ifacespython Not present.
 *
 * @param triangulation the triangulation upon which these validity constraints
 * will be based.
 * @param coords the coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalHypersurfaceList.
 * @return a newly allocated set of constraints.
 */
REGINA_API NEnumConstraintList* makeEmbeddedConstraints(
    const Dim4Triangulation* triangulation, HyperCoords coords);

/*@}*/

// Inline functions for NNormalHypersurfaceList

inline NNormalHypersurfaceList::~NNormalHypersurfaceList() {
    for_each(surfaces_.begin(), surfaces_.end(),
        FuncDelete<NNormalHypersurface>());
}

inline HyperCoords NNormalHypersurfaceList::coords() const {
    return coords_;
}

inline HyperList NNormalHypersurfaceList::which() const {
    return which_;
}

inline HyperAlg NNormalHypersurfaceList::algorithm() const {
    return algorithm_;
}

inline bool NNormalHypersurfaceList::isEmbeddedOnly() const {
    return which_.has(HS_EMBEDDED_ONLY);
}

inline size_t NNormalHypersurfaceList::size() const {
    return surfaces_.size();
}

inline const NNormalHypersurface* NNormalHypersurfaceList::hypersurface(
        size_t index) const {
    return surfaces_[index];
}

inline bool NNormalHypersurfaceList::dependsOnParent() const {
    return true;
}

inline NMatrixInt* NNormalHypersurfaceList::recreateMatchingEquations() const {
    return makeMatchingEquations(triangulation(), coords_);
}

inline NNormalHypersurfaceList::VectorIterator::VectorIterator() {
}

inline NNormalHypersurfaceList::VectorIterator::VectorIterator(
        const NNormalHypersurfaceList::VectorIterator& cloneMe) :
        it_(cloneMe.it_) {
}

inline NNormalHypersurfaceList::VectorIterator&
        NNormalHypersurfaceList::VectorIterator::operator =(
        const NNormalHypersurfaceList::VectorIterator& cloneMe) {
    it_ = cloneMe.it_;
    return *this;
}

inline bool NNormalHypersurfaceList::VectorIterator::operator ==(
        const NNormalHypersurfaceList::VectorIterator& other) const {
    return (it_ == other.it_);
}

inline bool NNormalHypersurfaceList::VectorIterator::operator !=(
        const NNormalHypersurfaceList::VectorIterator& other) const {
    return (it_ != other.it_);
}

inline const NNormalHypersurfaceVector*
        NNormalHypersurfaceList::VectorIterator::operator *() const {
    return (*it_)->rawVector();
}

inline NNormalHypersurfaceList::VectorIterator&
        NNormalHypersurfaceList::VectorIterator::operator ++() {
    ++it_;
    return *this;
}

inline NNormalHypersurfaceList::VectorIterator
        NNormalHypersurfaceList::VectorIterator::operator ++(int) {
    return NNormalHypersurfaceList::VectorIterator(it_++);
}

inline NNormalHypersurfaceList::VectorIterator&
        NNormalHypersurfaceList::VectorIterator::operator --() {
    --it_;
    return *this;
}

inline NNormalHypersurfaceList::VectorIterator
        NNormalHypersurfaceList::VectorIterator::operator --(int) {
    return NNormalHypersurfaceList::VectorIterator(it_--);
}

inline NNormalHypersurfaceList::VectorIterator::VectorIterator(
        const std::vector<NNormalHypersurface*>::const_iterator& i) : it_(i) {
}

inline NNormalHypersurfaceList::VectorIterator
        NNormalHypersurfaceList::beginVectors() const {
    return VectorIterator(surfaces_.begin());
}

inline NNormalHypersurfaceList::VectorIterator
        NNormalHypersurfaceList::endVectors() const {
    return VectorIterator(surfaces_.end());
}

inline NNormalHypersurfaceList::HypersurfaceInserter::HypersurfaceInserter(
        NNormalHypersurfaceList& list, Dim4Triangulation* owner) :
        list_(&list), owner_(owner) {
}

inline NNormalHypersurfaceList::HypersurfaceInserter::HypersurfaceInserter(
        const HypersurfaceInserter& cloneMe) : list_(cloneMe.list_),
        owner_(cloneMe.owner_) {
}


inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator =(
        const HypersurfaceInserter& cloneMe) {
    list_ = cloneMe.list_;
    owner_ = cloneMe.owner_;
    return *this;
}

inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator =(
        NNormalHypersurface* surface) {
    list_->surfaces_.push_back(surface);
    return *this;
}

inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator =(
        NNormalHypersurfaceVector* vector) {
    list_->surfaces_.push_back(new NNormalHypersurface(owner_, vector));
    return *this;
}

inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator *() {
    return *this;
}

inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator ++() {
    return *this;
}

inline NNormalHypersurfaceList::HypersurfaceInserter&
        NNormalHypersurfaceList::HypersurfaceInserter::operator ++(int) {
    return *this;
}

inline NNormalHypersurfaceList::NNormalHypersurfaceList(HyperCoords coords,
        HyperList which, HyperAlg algorithm) :
        coords_(coords), which_(which), algorithm_(algorithm) {
}

inline NNormalHypersurfaceList::Enumerator::Enumerator(
        NNormalHypersurfaceList* list, Dim4Triangulation* triang,
        NProgressTracker* tracker) :
        list_(list), triang_(triang), tracker_(tracker) {
}

} // namespace regina

#endif

