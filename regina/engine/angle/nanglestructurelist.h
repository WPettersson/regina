
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

/*! \file angle/nanglestructurelist.h
 *  \brief Contains a packet representing a collection of angle
 *  structures on a triangulation.
 */

#ifndef __NANGLESTRUCTURELIST_H
#ifndef __DOXYGEN
#define __NANGLESTRUCTURELIST_H
#endif

#include <algorithm>
#include <iterator>
#include <vector>
#include "regina-core.h"
#include "angle/nanglestructure.h"
#include "packet/npacket.h"
#include "utilities/memutils.h"
#include "utilities/nproperty.h"
#include "utilities/nthread.h"

namespace regina {

class NProgressManager;
class NXMLPacketReader;
class NXMLAngleStructureListReader;

/**
 * \weakgroup angle
 * @{
 */

/**
 * A packet representing a collection of angle structures on a triangulation.
 * Such a packet must always be a child packet of the triangulation on
 * which the angle structures lie.  If this triangulation changes, the
 * information contained in this packet will become invalid.
 *
 * Angle structure lists should be created using the routine enumerate(),
 * which is new as of Regina 3.95.
 *
 * \testpart
 */
class REGINA_API NAngleStructureList :
        public NPacket, public NFilePropertyReader {
    public:
        static const int packetType;

    private:
        std::vector<NAngleStructure*> structures;
            /**< Contains the angle structures stored in this packet. */
        bool tautOnly_;
            /**< Stores whether we are only interested in taut structures.
                 This is an option selected by the user before enumeration
                 takes place. */

        mutable NProperty<bool> doesSpanStrict;
            /**< Does the convex span of this list include a strict
                 angle structure?
                 This is determined by looking at the output angle structues
                 after enumeration has taken place. */
        mutable NProperty<bool> doesSpanTaut;
            /**< Does this list include a taut structure?
                 This is determined by looking at the output angle structues
                 after enumeration has taken place. */

    public:
        /**
         * Destroys this list and all the angle structures within.
         */
        virtual ~NAngleStructureList();

        /**
         * Returns the triangulation on which these angle structures
         * lie.
         *
         * @return the corresponding triangulation.
         */
        NTriangulation* getTriangulation() const;

        /**
         * Returns whether this list was produced by enumerating taut angle
         * structures only.
         *
         * @return \c true if this list was produced by enumerating
         * taut angle structures only, or \c false if the enumeration
         * procedure allowed for any angle structures.
         */
        bool isTautOnly() const;

        /**
         * Returns the number of angle structures stored in this list.
         *
         * @return the number of angle structures.
         */
        unsigned long getNumberOfStructures() const;
        /**
         * Returns the angle structure at the requested index in this
         * list.
         *
         * @param index the index of the requested angle structure in
         * this list; this must be between 0 and
         * getNumberOfStructures()-1 inclusive.
         * @return the angle structure at the requested index.
         */
        const NAngleStructure* getStructure(unsigned long index) const;

        /**
         * Determines whether any convex combination of the angle
         * structures in this list is a strict angle structure.
         * See NAngleStructure::isStrict() for details on strict angle
         * structures.
         *
         * @return \c true if and only if a strict angle structure can
         * be produced.
         */
        bool spansStrict() const;
        /**
         * Determines whether any angle structure in this list is a
         * taut structure.  Because taut structures always appear as
         * vertices of the angle structure solution space, this routine
         * is equivalent to testing whether any convex combination of
         * the angle structures in this list is a taut structure.
         *
         * See NAngleStructure::isTaut() for details on taut
         * structures.
         *
         * @return \c true if and only if a taut structure can be produced.
         */
        bool spansTaut() const;

        /**
         * Determines whether any convex combination of the angle
         * structures in this list is a strict angle structure.
         *
         * \deprecated This routine will be removed in a future version
         * of Regina.  Users should switch to the identical routine
         * spansStrict() instead.
         *
         * @return \c true if and only if a strict angle structure can
         * be produced.
         */
        bool allowsStrict() const;

        /**
         * Determines whether any angle structure in this list is a
         * taut structure.
         *
         * \deprecated This routine will be removed in a future version
         * of Regina.  Users should switch to the identical routine
         * spansTaut() instead.
         *
         * @return \c true if and only if a taut angle structure can
         * be produced.
         */
        bool allowsTaut() const;

        /**
         * Enumerates all angle structures on the given triangulation.
         * A list containing all vertices of the angle structure solution
         * space will be returned.
         *
         * The option is offered to find only taut structures (which are
         * considerably faster to enumerate) instead of enumerating all
         * vertex angle structures.  See the \a tautOnly argument below.
         *
         * The angle structure list that is created will be inserted as the
         * last child of the given triangulation.  This triangulation \b must
         * remain the parent of this angle structure list, and must not
         * change while this angle structure list remains in existence.
         *
         * If a progress manager is passed, the angle structure
         * enumeration will take place in a new thread and this routine
         * will return immediately.  The NProgress object assigned to
         * this progress manager is guaranteed to be of the class
         * NProgressNumber.
         *
         * If no progress manager is passed, the enumeration will run
         * in the current thread and this routine will return only when
         * the enumeration is complete.  Note that this enumeration can
         * be extremely slow for larger triangulations.
         *
         * @param owner the triangulation for which the vertex
         * angle structures will be enumerated.
         * @param tautOnly \c true if only taut structures are to be
         * enuemrated, or \c false if we should enumerate all vertices
         * of the angle structure solution space; this defaults to \c false.
         * @param manager a progress manager through which progress will
         * be reported, or 0 if no progress reporting is required.  If
         * non-zero, \a manager must point to a progress manager for
         * which NProgressManager::isStarted() is still \c false.
         * @return the newly created angle structure list.  Note that if
         * a progress manager is passed then this list may not be completely
         * filled when this routine returns.  If a progress manager is
         * passed and a new thread could not be started, this routine
         * returns 0 (and no angle structure list is created).
         */
        static NAngleStructureList* enumerate(NTriangulation* owner,
            bool tautOnly = false, NProgressManager* manager = 0);

        virtual int getPacketType() const;
        virtual std::string getPacketTypeName() const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeTextLong(std::ostream& out) const;
        static NXMLPacketReader* getXMLReader(NPacket* parent);
        virtual void writePacket(NFile& out) const;
        static NAngleStructureList* readPacket(NFile& in, NPacket* parent);
        virtual bool dependsOnParent() const;

    protected:
        /**
         * Creates a new empty angle structure list.
         * All properties are marked as unknown.
         *
         * @param tautOnly \c true if only taut structures are to be
         * enuemrated (when the time comes for enumeration to be performed),
         * or \c false if we should enumerate all vertices of the angle
         * structure solution space.
         */
        NAngleStructureList(bool tautOnly);

        virtual NPacket* internalClonePacket(NPacket* parent) const;
        virtual void writeXMLPacketData(std::ostream& out) const;
        virtual void readIndividualProperty(NFile& infile, unsigned propType);

        /**
         * Calculate whether the convex span of this list includes a
         * strict angle structure.
         */
        void calculateSpanStrict() const;
        /**
         * Calculate whether the convex span of this list includes a
         * taut structure.
         */
        void calculateSpanTaut() const;

        /**
         * An output iterator used to insert angle structures into an
         * NAngleStructureList.
         *
         * Objects of type <tt>NAngleStructure*</tt> and
         * <tt>NAngleStructureVector*</tt> can be assigned to this
         * iterator.  In the latter case, a surrounding NAngleStructure
         * will be automatically created.
         */
        struct StructureInserter : public std::iterator<
                std::output_iterator_tag, NAngleStructureVector*> {
            NAngleStructureList* list;
                /**< The list into which angle structures will be inserted. */
            NTriangulation* owner;
                /**< The triangulation on which the angle structures to
                 *   be inserted lie. */

            /**
             * Creates a new uninitialised output iterator.
             *
             * \warning This iterator must not be used until its
             * structure list and triangulation have been initialised.
             */
            StructureInserter();
            /**
             * Creates a new output iterator.  The member variables of
             * this iterator will be initialised according to the
             * parameters passed to this constructor.
             *
             * @param newList the list into which angle structures will
             * be inserted.
             * @param newOwner the triangulation on which the structures
             * to be inserted lie.
             */
            StructureInserter(NAngleStructureList& newList,
                NTriangulation* newOwner);
            /**
             * Creates a new output iterator that is a clone of the
             * given iterator.
             *
             * @param cloneMe the output iterator to clone.
             */
            StructureInserter(const StructureInserter& cloneMe);

            /**
             * Sets this iterator to be a clone of the given output iterator.
             *
             * @param cloneMe the output iterator to clone.
             * @return this output iterator.
             */
            StructureInserter& operator =(const StructureInserter& cloneMe);

            /**
             * Appends an angle structure to the end of the appropriate
             * structure list.
             *
             * The given angle structure will be deallocated with the
             * other angle structures in this list when the list is
             * eventually destroyed.
             *
             * @param structure the angle structure to insert.
             * @return this output iterator.
             */
            StructureInserter& operator =(NAngleStructure* structure);
            /**
             * Appends the angle structure corresponding to the given
             * vector to the end of the appropriate structure list.
             *
             * The given vector will be owned by the newly created
             * angle structure and will be deallocated with the
             * other angle structures in this list when the list is
             * eventually destroyed.
             *
             * @param vector the vector of the angle structure to insert.
             * @return this output iterator.
             */
            StructureInserter& operator =(NAngleStructureVector* vector);

            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            StructureInserter& operator *();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            StructureInserter& operator ++();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            StructureInserter& operator ++(int);
        };

    private:
        /**
         * A thread class that actually performs the angle structure
         * enumeration.
         */
        class Enumerator : public NThread {
            private:
                NAngleStructureList* list;
                    /**< The angle structure list to be filled. */
                NTriangulation* triang;
                    /**< The triangulation upon which this angle
                         structure list will be based. */
                NProgressManager* manager;
                    /**< The progress manager through which progress is
                         reported, or 0 if no progress manager is in use. */

            public:
                /**
                 * Creates a new enumerator thread with the given
                 * parameters.
                 *
                 * @param newList the angle structure list to be filled.
                 * @param useTriang the triangulation upon which this
                 * angle structure list will be based.
                 * @param useManager the progress manager to use for
                 * progress reporting, or 0 if progress reporting is not
                 * required.
                 */
                Enumerator(NAngleStructureList* newList,
                    NTriangulation* useTriang, NProgressManager* useManager);

                void* run(void*);
        };

    friend class regina::NXMLAngleStructureListReader;
};

/*@}*/

// Inline functions for NAngleStructureList

inline NAngleStructureList::~NAngleStructureList() {
    for_each(structures.begin(), structures.end(),
        FuncDelete<NAngleStructure>());
}

inline bool NAngleStructureList::isTautOnly() const {
    return tautOnly_;
}

inline unsigned long NAngleStructureList::getNumberOfStructures() const {
    return structures.size();
}

inline const NAngleStructure* NAngleStructureList::getStructure(
        unsigned long index) const {
    return structures[index];
}

inline bool NAngleStructureList::spansStrict() const {
    if (! doesSpanStrict.known())
        calculateSpanStrict();
    return doesSpanStrict.value();
}

inline bool NAngleStructureList::spansTaut() const {
    if (! doesSpanTaut.known())
        calculateSpanTaut();
    return doesSpanTaut.value();
}

inline bool NAngleStructureList::allowsStrict() const {
    return spansStrict();
}

inline bool NAngleStructureList::allowsTaut() const {
    return spansTaut();
}

inline bool NAngleStructureList::dependsOnParent() const {
    return true;
}

inline NAngleStructureList::NAngleStructureList(bool tautOnly) :
        tautOnly_(tautOnly) {
}

inline NAngleStructureList::StructureInserter::StructureInserter() : list(0),
        owner(0) {
}

inline NAngleStructureList::StructureInserter::StructureInserter(
        NAngleStructureList& newList, NTriangulation* newOwner) :
        list(&newList), owner(newOwner) {
}

inline NAngleStructureList::StructureInserter::StructureInserter(
        const StructureInserter& cloneMe) : list(cloneMe.list),
        owner(cloneMe.owner) {
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator =(
        const StructureInserter& cloneMe) {
    list = cloneMe.list;
    owner = cloneMe.owner;
    return *this;
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator =(
        NAngleStructure* structure) {
    list->structures.push_back(structure);
    return *this;
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator =(
        NAngleStructureVector* vector) {
    list->structures.push_back(new NAngleStructure(owner, vector));
    return *this;
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator *() {
    return *this;
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator ++() {
    return *this;
}

inline NAngleStructureList::StructureInserter&
        NAngleStructureList::StructureInserter::operator ++(int) {
    return *this;
}

inline NAngleStructureList::Enumerator::Enumerator(NAngleStructureList* newList,
        NTriangulation* useTriang, NProgressManager* useManager) :
        list(newList), triang(useTriang), manager(useManager) {
}

} // namespace regina

#endif

