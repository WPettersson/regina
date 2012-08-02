
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
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

// Regina core includes:
#include "packet/npacket.h"

// UI includes:
#include "packetchooser.h"
#include "packetfilter.h"
#include "packetmanager.h"

#include <algorithm>

using regina::NPacket;

PacketChooser::PacketChooser(regina::NPacket* newSubtree,
        QWidget* parent) :
        QComboBox(parent), subtree(newSubtree), filter(0),
        onAutoUpdate(false), isUpdating(false) {
    setMinimumContentsLength(30);
    setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    fill(false, 0);
}

PacketChooser::PacketChooser(regina::NPacket* newSubtree,
        PacketFilter* newFilter, QWidget* parent) :
        QComboBox(parent), subtree(newSubtree), filter(newFilter),
        onAutoUpdate(false), isUpdating(false) {
    setMinimumContentsLength(30);
    setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    fill(false, 0);
}

PacketChooser::PacketChooser(regina::NPacket* newSubtree,
        PacketFilter* newFilter, bool allowNone,
        regina::NPacket* initialSelection, QWidget* parent) :
        QComboBox(parent), subtree(newSubtree), filter(newFilter),
        onAutoUpdate(false), isUpdating(false) {
    setMinimumContentsLength(30);
    setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    fill(allowNone, initialSelection);
}

PacketChooser::~PacketChooser() {
    if (filter)
        delete filter;
}

NPacket* PacketChooser::selectedPacket() {
    if (count() == 0)
        return 0;
    int curr = currentIndex();
    return (curr < 0 ? 0 : packets[curr]);
}

void PacketChooser::selectPacket(regina::NPacket* packet) {
    int index = 0;
    std::vector<regina::NPacket*>::const_iterator it = packets.begin();
    while (it != packets.end()) {
        if ((*it) == packet) {
            setCurrentIndex(index);
            return;
        }
        ++it;
        ++index;
    }

    // Not found.
    if (! packets.empty())
        setCurrentIndex(0);
    return;
}

void PacketChooser::setAutoUpdate(bool shouldAutoUpdate) {
    if (onAutoUpdate == shouldAutoUpdate)
        return;

    onAutoUpdate = shouldAutoUpdate;
    if (onAutoUpdate)
        for (std::vector<regina::NPacket*>::iterator it = packets.begin();
                it != packets.end(); it++)
            if (*it)
                (*it)->listen(this);
    else
        unregisterFromAllPackets();
}

void PacketChooser::packetWasRenamed(regina::NPacket* renamed) {
    // Just rename the item that was changed.
    std::vector<regina::NPacket*>::iterator it = std::find(
        packets.begin(), packets.end(), renamed);

    if (it != packets.end()) {
        // This may trigger a refreshContents(), but that's okay since
        // we're at the end of the routine.
        int index = it - packets.begin();
        setItemIcon(index, PacketManager::icon(renamed));
        setItemText(index, renamed->getPacketLabel().c_str());
    }
}

void PacketChooser::packetToBeDestroyed(regina::NPacket* toDestroy) {
    // Just remove the item that is being destroyed.
    std::vector<regina::NPacket*>::iterator it = std::find(
        packets.begin(), packets.end(), toDestroy);

    if (it != packets.end()) {
        // Make sure the call to removeItem() comes last since it could
        // trigger a refreshContents().
        long destroyIndex = it - packets.begin();
        long currIndex = currentIndex();

        packets.erase(it);
        if (destroyIndex == currIndex) {
            // We know count() > 0 since currentItem() exists.
            setCurrentIndex(0);

            // If the item to destroy *is* 0, this should just fall through
            // to whatever's next when we remove it from the chooser.
        } else if (destroyIndex < currIndex) {
            // The selected item has moved up the list.
            setCurrentIndex(currIndex - 1);
        }

        removeItem(destroyIndex);

        // Don't bother unlistening; this will happen in the packet
        // destructor anyway.
    }
}

void PacketChooser::refreshContents() {
    if (isUpdating)
        return;

    // Don't change anything if we don't have to.
    if (verify())
        return;

    // Do a straight empty-and-refill.
    // There are more streamlined ways of doing it, but we'll come to it
    // when we come to it.
    isUpdating = true;

    // Remember how it used to look.
    NPacket* remember = selectedPacket();
    bool allowNone = ((! packets.empty()) && (! packets[0]));

    // Empty the combo box.
    // Empty from the end in case it's stored as a vector deep inside.
    if (onAutoUpdate)
        unregisterFromAllPackets();

    clear();
    packets.clear();

    // Fill it again.
    fill(allowNone, remember);

    isUpdating = false;
}

void PacketChooser::fill(bool allowNone, NPacket* select) {
    // Insert the None entry if appropriate.
    if (allowNone) {
        addItem(tr("<None>"));
        packets.push_back(0);

        if (select == 0)
            setCurrentIndex(0);
    }

    // Insert the regular packets.
    regina::NPacket* p = subtree;
    while (p && subtree->isGrandparentOf(p)) {
        if ((! filter) || (filter->accept(p))) {
            addItem(PacketManager::icon(p), p->getPacketLabel().c_str());
            packets.push_back(p);
            if (onAutoUpdate)
                p->listen(this);
            if (p == select)
                setCurrentIndex(count() - 1);
        }
        p = p->nextTreePacket();
    }
}

bool PacketChooser::verify() {
    regina::NPacket* p = subtree;
    std::vector<regina::NPacket*>::const_iterator it = packets.begin();

    // Ignore the None entry if it exists.
    if (it != packets.end() && (! *it))
        it++;

    // Now match the packets up one by one.
    while (it != packets.end() || p != 0) {
        // Are we ignoring this packet?
        if (p && filter && ! filter->accept(p)) {
            p = p->nextTreePacket();
            continue;
        }

        // Out of packets?
        if (it == packets.end())
            return false;
        // Out of combo box entries?
        if (! p)
            return false;
        // Mismatched entries?
        if (p != *it)
            return false;

        // All good.
        it++;
        p = p->nextTreePacket();
    }

    return true;
}

