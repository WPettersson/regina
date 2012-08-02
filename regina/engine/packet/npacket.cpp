
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

#include <set>
#include <sstream>
#include "engine.h"
#include "packet/npacket.h"
#include "packet/npacketlistener.h"
#include "utilities/xmlutils.h"

namespace regina {

NPacket::~NPacket() {
    inDestructor = true;

    // Orphan this packet before doing anything else.
    // The destructor can lead to callbacks for packet listeners, which
    // might in turn involve tree traversal.  It can't be good for
    // anyone to start querying packets whose destructors are already
    // being carried out.
    if (treeParent)
        makeOrphan();

    // Destroy all descendants.
    // Note that the NPacket destructor now orphans the packet as well.
    while(firstTreeChild)
        delete firstTreeChild;

    // Fire a packet event and unregister all listeners.
    fireDestructionEvent();
}

void NPacket::setPacketLabel(const std::string& newLabel) {
    fireEvent(&NPacketListener::packetToBeRenamed);
    packetLabel = newLabel;
    fireEvent(&NPacketListener::packetWasRenamed);
}

bool NPacket::listen(NPacketListener* listener) {
    if (! listeners.get())
        listeners.reset(new std::set<NPacketListener*>());

    listener->packets.insert(this);
    return listeners->insert(listener).second;
}

bool NPacket::unlisten(NPacketListener* listener) {
    if (! listeners.get())
        return false;

    listener->packets.erase(this);
    return listeners->erase(listener);
}

NPacket* NPacket::getTreeMatriarch() const {
    NPacket* p = const_cast<NPacket*>(this);
    while (p->treeParent)
        p = p->treeParent;
    return p;
}

void NPacket::insertChildFirst(NPacket* child) {
    fireEvent(&NPacketListener::childToBeAdded, child);

    child->treeParent = this;
    child->prevTreeSibling = 0;
    child->nextTreeSibling = firstTreeChild;

    if (firstTreeChild) {
        firstTreeChild->prevTreeSibling = child;
        firstTreeChild = child;
    } else {
        firstTreeChild = child;
        lastTreeChild = child;
    }

    fireEvent(&NPacketListener::childWasAdded, child);
}

void NPacket::insertChildLast(NPacket* child) {
    fireEvent(&NPacketListener::childToBeAdded, child);

    child->treeParent = this;
    child->prevTreeSibling = lastTreeChild;
    child->nextTreeSibling = 0;

    if (lastTreeChild) {
        lastTreeChild->nextTreeSibling = child;
        lastTreeChild = child;
    } else {
        firstTreeChild = child;
        lastTreeChild = child;
    }

    fireEvent(&NPacketListener::childWasAdded, child);
}

void NPacket::insertChildAfter(NPacket* newChild, NPacket* prevChild) {
    fireEvent(&NPacketListener::childToBeAdded, newChild);

    if (prevChild == 0)
        insertChildFirst(newChild);
    else {
        newChild->treeParent = this;
        newChild->nextTreeSibling = prevChild->nextTreeSibling;
        newChild->prevTreeSibling = prevChild;
        prevChild->nextTreeSibling = newChild;
        if (newChild->nextTreeSibling)
            newChild->nextTreeSibling->prevTreeSibling = newChild;
        else
            lastTreeChild = newChild;
    }

    fireEvent(&NPacketListener::childWasAdded, newChild);
}

void NPacket::makeOrphan() {
    if (! treeParent)
        return;
    NPacket* oldParent = treeParent;

    oldParent->fireEvent(&NPacketListener::childToBeRemoved,
        this, inDestructor);

    if (treeParent->firstTreeChild == this)
        treeParent->firstTreeChild = nextTreeSibling;
    else
        prevTreeSibling->nextTreeSibling = nextTreeSibling;

    if (treeParent->lastTreeChild == this)
        treeParent->lastTreeChild = prevTreeSibling;
    else
        nextTreeSibling->prevTreeSibling = prevTreeSibling;

    treeParent = 0;

    oldParent->fireEvent(&NPacketListener::childWasRemoved,
        this, inDestructor);
}

void NPacket::reparent(NPacket* newParent, bool first) {
    if (treeParent)
        makeOrphan();

    if (first)
        newParent->insertChildFirst(this);
    else
        newParent->insertChildLast(this);
}

void NPacket::moveUp(unsigned steps) {
    if (steps == 0 || ! prevTreeSibling)
        return;

    // This packet is not the first packet in the child list.
    treeParent->fireEvent(&NPacketListener::childrenToBeReordered);

    NPacket* prev = prevTreeSibling;
    while (prev && steps) {
        prev = prev->prevTreeSibling;
        steps--;
    }

    // Pull us out of the tree.
    if (nextTreeSibling)
        nextTreeSibling->prevTreeSibling = prevTreeSibling;
    else
        treeParent->lastTreeChild = prevTreeSibling;

    prevTreeSibling->nextTreeSibling = nextTreeSibling;

    // Reinsert ourselves into the tree.
    prevTreeSibling = prev;
    nextTreeSibling =
        (prev ? prev->nextTreeSibling : treeParent->firstTreeChild);
    nextTreeSibling->prevTreeSibling = this;
    if (prev)
        prev->nextTreeSibling = this;
    else
        treeParent->firstTreeChild = this;

    treeParent->fireEvent(&NPacketListener::childrenWereReordered);
}

void NPacket::moveDown(unsigned steps) {
    if (steps == 0 || ! nextTreeSibling)
        return;

    // This packet is not the last packet in the child list.
    treeParent->fireEvent(&NPacketListener::childrenToBeReordered);

    NPacket* next = nextTreeSibling;
    while (next && steps) {
        next = next->nextTreeSibling;
        steps--;
    }

    // Pull us out of the tree.
    if (prevTreeSibling)
        prevTreeSibling->nextTreeSibling = nextTreeSibling;
    else
        treeParent->firstTreeChild = nextTreeSibling;

    nextTreeSibling->prevTreeSibling = prevTreeSibling;

    // Reinsert ourselves into the tree.
    nextTreeSibling = next;
    prevTreeSibling =
        (next ? next->prevTreeSibling : treeParent->lastTreeChild);
    prevTreeSibling->nextTreeSibling = this;
    if (next)
        next->prevTreeSibling = this;
    else
        treeParent->lastTreeChild = this;

    treeParent->fireEvent(&NPacketListener::childrenWereReordered);
}

void NPacket::moveToFirst() {
    if (! prevTreeSibling)
        return;

    // This packet is not the first packet in the child list.
    treeParent->fireEvent(&NPacketListener::childrenToBeReordered);

    // Pull us out of the tree.
    if (nextTreeSibling)
        nextTreeSibling->prevTreeSibling = prevTreeSibling;
    else
        treeParent->lastTreeChild = prevTreeSibling;

    prevTreeSibling->nextTreeSibling = nextTreeSibling;

    // Reinsert ourselves into the tree.
    treeParent->firstTreeChild->prevTreeSibling = this;
    nextTreeSibling = treeParent->firstTreeChild;
    prevTreeSibling = 0;
    treeParent->firstTreeChild = this;

    treeParent->fireEvent(&NPacketListener::childrenWereReordered);
}

void NPacket::moveToLast() {
    if (! nextTreeSibling)
        return;

    // This packet is not the last packet in the child list.
    treeParent->fireEvent(&NPacketListener::childrenToBeReordered);

    // Pull us out of the tree.
    if (prevTreeSibling)
        prevTreeSibling->nextTreeSibling = nextTreeSibling;
    else
        treeParent->firstTreeChild = nextTreeSibling;

    nextTreeSibling->prevTreeSibling = prevTreeSibling;

    // Reinsert ourselves into the tree.
    treeParent->lastTreeChild->nextTreeSibling = this;
    prevTreeSibling = treeParent->lastTreeChild;
    nextTreeSibling = 0;
    treeParent->lastTreeChild = this;

    treeParent->fireEvent(&NPacketListener::childrenWereReordered);
}

void NPacket::sortChildren() {
    // Run through the packets from largest to smallest, moving each to
    // the beginning of the child list in turn.
    NPacket* endpoint = 0;
    NPacket* current;
    NPacket* largest;

    fireEvent(&NPacketListener::childrenToBeReordered);
    while (1) {
        // Put current at the beginning of the clump of yet-unsorted children.
        if (! endpoint)
            current = firstTreeChild;
        else
            current = endpoint->nextTreeSibling;
        if (! current)
            break;

        // Find the largest amongst the yet-unsorted children.
        largest = current;
        current = current->nextTreeSibling;
        while (current) {
            if (current->getPacketLabel() > largest->getPacketLabel())
                largest = current;
            current = current->nextTreeSibling;
        }

        // Move the largest to the front of the list.
        if (firstTreeChild != largest) {
            // We know that largest has a previous sibling.
            largest->prevTreeSibling->nextTreeSibling =
                largest->nextTreeSibling;

            if (largest->nextTreeSibling)
                largest->nextTreeSibling->prevTreeSibling =
                    largest->prevTreeSibling;
            else
                lastTreeChild = largest->prevTreeSibling;

            firstTreeChild->prevTreeSibling = largest;
            largest->nextTreeSibling = firstTreeChild;
            largest->prevTreeSibling = 0;
            firstTreeChild = largest;
        }

        if (! endpoint)
            endpoint = largest;
    }

    fireEvent(&NPacketListener::childrenWereReordered);
}

void NPacket::swapWithNextSibling() {
    if (! nextTreeSibling)
        return;

    treeParent->fireEvent(&NPacketListener::childrenToBeReordered);

    if (prevTreeSibling)
        prevTreeSibling->nextTreeSibling = nextTreeSibling;
    else
        treeParent->firstTreeChild = nextTreeSibling;

    if (nextTreeSibling->nextTreeSibling)
        nextTreeSibling->nextTreeSibling->prevTreeSibling = this;
    else
        treeParent->lastTreeChild = this;

    NPacket* other = nextTreeSibling;

    nextTreeSibling = other->nextTreeSibling;
    other->prevTreeSibling = prevTreeSibling;
    prevTreeSibling = other;
    other->nextTreeSibling = this;

    treeParent->fireEvent(&NPacketListener::childrenWereReordered);
}

NPacket* NPacket::nextTreePacket() {
    if (firstTreeChild)
        return firstTreeChild;
    if (nextTreeSibling)
        return nextTreeSibling;
    NPacket* tmp = treeParent;
    while (tmp) {
        if (tmp->nextTreeSibling)
            return tmp->nextTreeSibling;
        tmp = tmp->treeParent;
    }
    return 0;
}

const NPacket* NPacket::nextTreePacket() const {
    if (firstTreeChild)
        return firstTreeChild;
    if (nextTreeSibling)
        return nextTreeSibling;
    NPacket* tmp = treeParent;
    while (tmp) {
        if (tmp->nextTreeSibling)
            return tmp->nextTreeSibling;
        tmp = tmp->treeParent;
    }
    return 0;
}

NPacket* NPacket::firstTreePacket(const std::string& type) {
    if (getPacketTypeName() == type)
        return this;
    return nextTreePacket(type);
}

const NPacket* NPacket::firstTreePacket(const std::string& type) const {
    if (getPacketTypeName() == type)
        return this;
    return nextTreePacket(type);
}

NPacket* NPacket::nextTreePacket(const std::string& type) {
    NPacket* ans = nextTreePacket();
    while (ans) {
        if (ans->getPacketTypeName() == type)
            return ans;
        ans = ans->nextTreePacket();
    }
    return 0;
}

const NPacket* NPacket::nextTreePacket(const std::string& type) const {
    const NPacket* ans = nextTreePacket();
    while (ans) {
        if (ans->getPacketTypeName() == type)
            return ans;
        ans = ans->nextTreePacket();
    }
    return 0;
}

NPacket* NPacket::findPacketLabel(const std::string& label) {
    if (packetLabel == label)
        return this;
    NPacket* tmp = firstTreeChild;
    NPacket* ans;
    while (tmp) {
        ans = tmp->findPacketLabel(label);
        if (ans)
            return ans;
        tmp = tmp->nextTreeSibling;
    }
    return 0;
}

const NPacket* NPacket::findPacketLabel(const std::string& label) const {
    if (packetLabel == label)
        return this;
    NPacket* tmp = firstTreeChild;
    NPacket* ans;
    while (tmp) {
        ans = tmp->findPacketLabel(label);
        if (ans)
            return ans;
        tmp = tmp->nextTreeSibling;
    }
    return 0;
}

unsigned NPacket::levelsDownTo(const NPacket* descendant) const {
    unsigned levels = 0;
    while (descendant != this) {
        descendant = descendant->treeParent;
        levels++;
    }
    return levels;
}

bool NPacket::isGrandparentOf(const NPacket* descendant) const {
    while (descendant) {
        if (descendant == this)
            return true;
        descendant = descendant->treeParent;
    }
    return false;
}

unsigned long NPacket::getNumberOfChildren() const {
    unsigned long tot = 0;
    for (NPacket* tmp = firstTreeChild; tmp; tmp = tmp->nextTreeSibling)
        tot++;
    return tot;
}

unsigned long NPacket::getTotalTreeSize() const {
    unsigned long tot = 1;
    for (NPacket* tmp = firstTreeChild; tmp; tmp = tmp->nextTreeSibling)
        tot += tmp->getTotalTreeSize();
    return tot;
}

bool NPacket::isPacketEditable() const {
    NPacket* tmp = firstTreeChild;
    while (tmp) {
        if (tmp->dependsOnParent())
            return false;
        tmp = tmp->nextTreeSibling;
    }
    return true;
}

NPacket* NPacket::clone(bool cloneDescendants, bool end) const {
    if (treeParent == 0)
        return 0;
    NPacket* ans = internalClonePacket(treeParent);
    ans->setPacketLabel(makeUniqueLabel(packetLabel + " - clone"));
    if (end)
        treeParent->insertChildLast(ans);
    else
        treeParent->insertChildAfter(ans, const_cast<NPacket*>(this));
    if (cloneDescendants)
        internalCloneDescendants(ans);
    return ans;
}

void NPacket::internalCloneDescendants(NPacket* parent) const {
    NPacket* child = firstTreeChild;
    NPacket* clone;
    while (child) {
        clone = child->internalClonePacket(parent);
        clone->setPacketLabel(makeUniqueLabel(child->packetLabel
            + " - clone"));
        parent->insertChildLast(clone);
        child->internalCloneDescendants(clone);
        child = child->nextTreeSibling;
    }
}

std::string NPacket::makeUniqueLabel(const std::string& base) const {
    const NPacket* topLevel = this;
    while (topLevel->treeParent)
        topLevel = topLevel->treeParent;

    if (! topLevel->findPacketLabel(base))
        return base;

    std::string ans;
    unsigned long extraInt = 2;
    while(1) {
        std::ostringstream out;
        out << ' ' << extraInt;
        ans = base + out.str();
        if (! topLevel->findPacketLabel(ans))
            return ans;
        else
            extraInt++;
    }
    return "";
}

bool NPacket::makeUniqueLabels(NPacket* reference) {
    NPacket* tree[3];
    if (reference) {
        tree[0] = reference;
        tree[1] = this;
        tree[2] = 0;
    } else {
        tree[0] = this;
        tree[1] = 0;
    }

    std::set<std::string> labels;

    int whichTree;
    NPacket* p;
    std::string label, newLabel;
    unsigned long extraInt;
    bool changed = false;
    for (whichTree = 0; tree[whichTree]; whichTree++)
        for (p = tree[whichTree]; p; p = p->nextTreePacket()) {
            label = p->getPacketLabel();
            if (! labels.insert(label).second) {
                extraInt = 1;
                do {
                    extraInt++;
                    std::ostringstream out;
                    out << ' ' << extraInt;
                    newLabel = label + out.str();
                } while (! labels.insert(newLabel).second);

                p->setPacketLabel(newLabel);
                changed = true;
            }
        }

    return changed;
}

bool NPacket::addTag(const std::string& tag) {
    fireEvent(&NPacketListener::packetToBeRenamed);
    if (! tags.get())
        tags.reset(new std::set<std::string>());

    bool ans = tags->insert(tag).second;
    fireEvent(&NPacketListener::packetWasRenamed);
    return ans;
}

bool NPacket::removeTag(const std::string& tag) {
    if (! tags.get())
        return false;

    fireEvent(&NPacketListener::packetToBeRenamed);
    bool ans = tags->erase(tag);
    fireEvent(&NPacketListener::packetWasRenamed);
    return ans;
}

void NPacket::removeAllTags() {
    if (tags.get() && ! tags->empty()) {
        fireEvent(&NPacketListener::packetToBeRenamed);
        tags->clear();
        fireEvent(&NPacketListener::packetWasRenamed);
    }
}

void NPacket::writeXMLFile(std::ostream& out) const {
    // Write the XML header.
    out << "<?xml version=\"1.0\"?>\n";

    // Write the regina data opening tag including engine version.
    out << "<reginadata engine=\"" << regina::getVersionString() << "\">\n";

    // Write the packet tree.
    writeXMLPacketTree(out);

    // Write the regina data closing tag.
    out << "</reginadata>\n";
}

void NPacket::fireEvent(void (NPacketListener::*event)(NPacket*)) {
    if (listeners.get()) {
        std::set<NPacketListener*>::const_iterator it = listeners->begin();
        while (it != listeners->end())
            ((*it++)->*event)(this);
    }
}

void NPacket::fireEvent(void (NPacketListener::*event)(NPacket*, NPacket*),
        NPacket* arg2) {
    if (listeners.get()) {
        std::set<NPacketListener*>::const_iterator it = listeners->begin();
        while (it != listeners->end())
            ((*it++)->*event)(this, arg2);
    }
}

void NPacket::fireEvent(void (NPacketListener::*event)(NPacket*,
        NPacket*, bool), NPacket* arg2, bool arg3) {
    if (listeners.get()) {
        std::set<NPacketListener*>::const_iterator it = listeners->begin();
        while (it != listeners->end())
            ((*it++)->*event)(this, arg2, arg3);
    }
}

void NPacket::fireDestructionEvent() {
    if (listeners.get()) {
        std::set<NPacketListener*>::const_iterator it;
        NPacketListener* tmp;
        while (! listeners->empty()) {
            it = listeners->begin();
            tmp = *it;

            // Unregister *before* we fire the event for each listener.
            // If we have a listener that deletes itself (or other listeners),
            // we don't want things to get nasty.
            listeners->erase(it);
            tmp->packets.erase(this);

            tmp->packetToBeDestroyed(this);
        }
    }
}

void NPacket::writeXMLPacketTree(std::ostream& out) const {
    using regina::xml::xmlEncodeSpecialChars;
    using regina::xml::xmlEncodeComment;

    // Write the packet opening tag including packet label and type.
    out << "<packet label=\"" << xmlEncodeSpecialChars(packetLabel) << "\"\n";
    out << "\ttype=\"" << getPacketTypeName() << "\" typeid=\""
        << getPacketType() << "\"\n";
    out << "\tparent=\"";
    if (treeParent)
        out << xmlEncodeSpecialChars(treeParent->packetLabel);
    out << "\">\n";

    // Write the internal packet data.
    writeXMLPacketData(out);

    // Write any packet tags.
    if (tags.get())
        for (std::set<std::string>::const_iterator it = tags->begin();
                it != tags->end(); it++)
            out << "  <tag name=\"" << xmlEncodeSpecialChars(*it) << "\"/>\n";

    // Write the child packets.
    for (NPacket* p = firstTreeChild; p; p = p->nextTreeSibling)
        p->writeXMLPacketTree(out);

    // Write the packet closing tag.
    out << "</packet> <!-- " << xmlEncodeComment(packetLabel)
        << " (" << xmlEncodeComment(getPacketTypeName()) << ") -->\n";
}

} // namespace regina

