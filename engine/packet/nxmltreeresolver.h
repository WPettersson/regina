
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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

/* end stub */

/*! \file packet/nxmltreeresolver.h
 *  \brief Support for resolving dangling packet references after a
 *  complete packet tree has been read from file.
 */

#ifndef __NXMLTREERESOLVER_H
#ifndef __DOXYGEN
#define __NXMLTREERESOLVER_H
#endif

#include "regina-core.h"
#include <list>
#include <map>

namespace regina {

class NPacket;

/**
 * \weakgroup packet
 * @{
 */

class NXMLTreeResolver;

/**
 * An individual task for resolving dangling packet references after an
 * XML data file has been read.
 *
 * See the NXMLTreeResolver class notes for an overview of how dangling
 * references and related issues are resolved, and the role that
 * NXMLTreeResolutionTask plays in this process.
 *
 * Specifically, if an individual NXMLPacketReader cannot
 * completely flesh out the internal data for a packet as the packet is
 * being read, it should construct a new NXMLTreeResolutionTask and
 * queue it to the master NXMLTreeResolver.  The NXMLTreeResolver will
 * then call resolve() for each queued task after the complete data file
 * has been read, at which point whatever information was missing when the
 * packet was initially read should now be available.
 *
 * Each packet reader that requires this machinery should subclass
 * NXMLTreeResolutionTask, and override resolve() to perform whatever
 * "fleshing out" procedure is required for its particular type of packet.
 */
class REGINA_API NXMLTreeResolutionTask {
    public:
        /**
         * A default construct that does nothing.
         */
        virtual ~NXMLTreeResolutionTask();
        /**
         * Called by NXMLTreeResolver after the entire data file has
         * been read.  Subclasses should override this routine to
         * perform whatever "fleshing out" is necessary for a packet
         * whose internal data is not yet complete.
         *
         * @param resolver the master resolver managing the resolution
         * process, as outlined in the NXMLTreeResolver class notes.
         */
        virtual void resolve(const NXMLTreeResolver& resolver) = 0;
};

/**
 * Provides a mechanism to resolve dangling packet references after a
 * complete packet tree has been read from an XML data file.
 *
 * There are situations in which, when reading an XML data file, the data
 * stored in an individual packet cannot be fully constructed until after
 * the entire data file has been read.  For instance, a packet might need to
 * store pointers or references to other packets that could appear later in
 * the packet tree (e.g., a script storing pointers to its variables).
 *
 * This problem is solved by the NXMLTreeResolver class.  The complete
 * process of reading an XML data file works as follows:
 *
 * - The top-level routine managing the file I/O constructs a new
 *   NXMLTreeResolver.  This resolver is then passed to each
 *   NXMLPacketReader in turn as each individual packet is read.
 *
 * - If an NXMLPacketReader is not able to fully flesh out its data
 *   because it requires information that is not yet available, it
 *   should create a new NXMLTreeResolutionTask and queue this task to
 *   the resolver via NXMLTreeResolver::queueTask().
 *
 * - Once the entire packet tree has been read, the top-level file I/O
 *   manager will call NXMLTreeResolver::resolve().  This will run
 *   NXMLTreeResolutionTask::resolve() for each task in turn, whereby any
 *   missing data for individual packets can be resolved.
 *
 * Each task should be an instance of a subclass of NXMLTreeResolutionTask,
 * whose virtual resolve() function is overridden to perform whatever
 * "fleshing out" work is required for the type of packet under consideration.
 */
class REGINA_API NXMLTreeResolver {
    public:
        typedef std::map<std::string, NPacket*> IDMap;
            /**< A type that maps internal IDs from the data file to the
                 corresponding packets.  See ids() for details. */

    private:
        IDMap ids_;
            /**< Maps internal IDs from the data file to the
                 corresponding packets. */
        std::list<NXMLTreeResolutionTask*> tasks_;
            /**< The list of tasks that have been queued for processing. */

    public:
        /**
         * Constructs a resolver with no tasks queued.
         */
        NXMLTreeResolver();
        /**
         * Destroys any tasks that were queued but not performed.
         */
        ~NXMLTreeResolver();

        /**
         * Queues a task for processing.  When the file I/O manager
         * calls resolve(), this will call NXMLTreeResolutionTask::resolve()
         * for each task that has been queued.
         *
         * This object will claim ownership of the given task, and will
         * destroy it after resolve() has been called (or, if resolve()
         * is never called, when this NXMLTreeResolver is destroyed).
         *
         * @param task the task to be queued.
         */
        void queueTask(NXMLTreeResolutionTask* task);
        /**
         * Stores the fact that the given packet is stored in the data
         * file using the given internal ID.  Associations between IDs
         * and packets can be queried through the ids() function.
         * See ids() for further information on internal IDs.
         *
         * This will be called automatically by NXMLPacketReader as it
         * processes packet tags in the data file.  Users and/or subclasses
         * of NXMLPacketReader do not need to call this function themselves.
         *
         * @param id the internal ID of the given packet, as stored in
         * the data file.
         * @param packet the corresponding packet.
         */
        void storeID(const std::string& id, NPacket* packet);

        /**
         * Returns the map from internal IDs to packets, as stored in
         * the data file.
         *
         * Packets in a data file may have individual string IDs stored
         * alongside them, in the \a id attribute of the
         * <tt>&lt;packet&gt;</tt> tag.  These strings are optional,
         * and do not need to be human-readable.
         * Although packets are not required to have IDs, any IDs that \e are
         * stored must be unique (i.e., two different packets cannot
         * share the same ID).
         *
         * Note that IDs read from the data file need not bear any
         * relation to the IDs that are returned from NPacket::internalID(),
         * although this is typically how they are constructed when a
         * file is saved.
         *
         * This map will be fleshed out as the data file is read.  In
         * particular, since each task runs NXMLTreeResolutionTask::resolve()
         * only after the entire tree has been read, tasks may assume that
         * this map contains all IDs that were explicitly stored in the
         * data file.
         *
         * Only packets with IDs will appear in this map (i.e., there may well
         * be packets in the data file that do not appear in this map at all).
         *
         * @return the map from internal file IDs to packets.
         */
        const IDMap& ids() const;

        /**
         * Calls NXMLTreeResolutionTask::resolve() for all queued tasks.
         *
         * The tasks will then be destroyed and removed from the queue
         * (so subsequent calls to resolve() are safe and will do nothing).
         */
        void resolve();
};

/*@}*/

// Inline functions for NXMLTreeResolutionTask

inline NXMLTreeResolutionTask::~NXMLTreeResolutionTask() {
}

// Inline functions for NXMLTreeResolver

inline NXMLTreeResolver::NXMLTreeResolver() {
}

inline NXMLTreeResolver::~NXMLTreeResolver() {
    for (std::list<NXMLTreeResolutionTask*>::iterator it = tasks_.begin();
            it != tasks_.end(); ++it)
        delete *it;
}

inline void NXMLTreeResolver::queueTask(NXMLTreeResolutionTask* task) {
    tasks_.push_back(task);
}

inline void NXMLTreeResolver::storeID(const std::string& id, NPacket* packet) {
    ids_.insert(std::make_pair(id, packet));
}

inline const NXMLTreeResolver::IDMap& NXMLTreeResolver::ids() const {
    return ids_;
}

inline void NXMLTreeResolver::resolve() {
    for (std::list<NXMLTreeResolutionTask*>::iterator it = tasks_.begin();
            it != tasks_.end(); ++it) {
        (*it)->resolve(*this);
        delete *it;
    }
    tasks_.clear();
}

} // namespace regina

#endif
