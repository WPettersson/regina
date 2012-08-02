
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

/*! \file packetui.h
 *  \brief Provides a basic infrastructure for packet interface components.
 */

#ifndef __PACKETUI_H
#define __PACKETUI_H

#include "packet/npacketlistener.h"

#include <QLinkedList>
#include <QWidget>

class PacketEditIface;
class PacketPane;
class PacketWindow;
class QAction;
class QLabel;
class QMenu;
class QToolButton;
class QTreeWidget;
class ReginaMain;

namespace regina {
    class NPacket;
};

/**
 * A packet-specific interface component for viewing or editing a packet.
 * Generic interface elements such as actions for refreshing or closing
 * a packet interface are not included.
 *
 * Different packet types should use different subclasses of PacketUI for
 * user interaction.  The PacketManager class is responsible for
 * creating an interface component appropriate for a given packet.
 *
 * Subclasses should call setDirty(true) whenever changes are made in
 * the interface.  Likewise, they should call setDirty(false) at the
 * end of their implementations of commit() and refresh().
 * Changes must never be made to the underlying packet except for within
 * the commit() routine.
 *
 * Each packet interface is either in read-write mode or in read-only
 * mode, and may be required to change modes throughout its life span.
 * See setReadWrite() for details.
 *
 * Subclasses will generally wish to override many of the PacketUI virtual
 * functions.
 */
class PacketUI {
    protected:
        /**
         * External components
         */
        PacketPane* enclosingPane;

    public:
        /**
         * Constructor and destructor.
         */
        PacketUI(PacketPane* newEnclosingPane);
        virtual ~PacketUI();

        /**
         * Return the packet that this pane is managing.
         * This routine should always return the same pointer throughout
         * the life of this object.
         */
        virtual regina::NPacket* getPacket() = 0;

        /**
         * Return the entire interface component.
         * This routine should always return the same pointer throughout
         * the life of this object.
         */
        virtual QWidget* getInterface() = 0;

        /**
         * Return details of the interface's interaction with standard
         * edit and clipboard operations.  This may be 0 if there is no
         * such interaction.
         *
         * The default implementation simply returns 0.
         */
        virtual PacketEditIface* getEditIface();

        /**
         * Return a list of actions specific to the particular type of
         * packet handled by this PacketUI subclass.  Such actions might
         * (for instance) query or manipulate packets of a particular
         * type.
         *
         * The default implementation of this routine simply returns an
         * empty list.
         */
        virtual const QLinkedList<QAction*>& getPacketTypeActions();

        /**
         * Return the label of the menu that should contain the actions
         * specific to this particular type of packet.  This label
         * may contain an ampersand (for keyboard shortcuts).
         */
        virtual QString getPacketMenuText() const = 0;

        /**
         * Store any changes currently made in this interface in the
         * underlying packet.
         *
         * Note that if a packet interface wishes to force a commit, it
         * should call enclosingPane->commit() and not this routine.
         * This will ensure that the enclosing pane can keep up to date
         * with what is taking place.
         *
         * No checking needs to be done as to whether the commit should
         * be allowed; this is taken care of by PacketPane::commit().
         *
         * This routine should call setDirty(false) once changes have
         * been made.
         */
        virtual void commit() = 0;

        /**
         * Update this interface to reflect the current contents of the
         * underlying packet.
         *
         * Note that if a packet interface wishes to force a refresh, it
         * should call enclosingPane->refresh() and not this routine.
         * This will ensure that the enclosing pane can keep up to date
         * with what is taking place.
         *
         * This routine should call setDirty(false) once the interface
         * has been updated.
         */
        virtual void refresh() = 0;

        /**
         * Modify this interface to be read-write or read-only according
         * to the given argument.
         *
         * This routine should never be called directly; instead
         * PacketPane::setReadWrite() should be used.
         *
         * If this interface is incapable of editing packets (e.g.,
         * interfaces for packet types that are inherently read-only
         * such as containers), this routine need not do anything.
         */
        virtual void setReadWrite(bool readWrite) = 0;

    protected:
        /**
         * Notifies external interface elements that this interface does
         * (or does not) contain changes that have not yet been committed
         * to the underlying packet.
         *
         * This routine must be called whenever changes are made in the
         * interface, and must also be called at the end of the
         * implementations of commit() and refresh().
         *
         * This routine should generally not be overridden by
         * subclasses.  If it is however, any reimplementations must
         * call the parent implementation.
         */
        virtual void setDirty(bool newDirty);

    private:
        /**
         * An empty action list.
         */
        static QLinkedList<QAction*> noActions;
};

/**
 * A packet interface that does not allow changes to be made.
 * Packets types that are inherently read-only (such as containers)
 * will probably want to use a subclass of PacketReadOnlyUI for their
 * interfaces.
 */
class PacketReadOnlyUI : public PacketUI {
    public:
        /**
         * Constructor.
         */
        PacketReadOnlyUI(PacketPane* newEnclosingPane);

        /**
         * An implementation of commit() that does nothing but call
         * setDirty(false).
         */
        virtual void commit();

        /**
         * An implementation of setReadWrite() that does nothing
         * whatsoever.
         */
        virtual void setReadWrite(bool readWrite);
};

/**
 * A packet interface that simply displays a given error message.
 */
class ErrorPacketUI : public PacketReadOnlyUI {
    private:
        regina::NPacket* packet;
        QLabel* label;

    public:
        /**
         * Constructor.
         */
        ErrorPacketUI(regina::NPacket* newPacket,
            PacketPane* newEnclosingPane, const QString& errorMessage);

        /**
         * Implementations of PacketUI virtual functions.
         */
        virtual regina::NPacket* getPacket();
        virtual QWidget* getInterface();
        virtual QString getPacketMenuText() const;
        virtual void refresh();
};

/**
 * A packet interface that should be used for unknown packet types.
 * A simple message is displayed indicating that the packet cannot be
 * viewed.
 */
class DefaultPacketUI : public ErrorPacketUI {
    public:
        /**
         * Constructor.
         */
        DefaultPacketUI(regina::NPacket* newPacket,
            PacketPane* newEnclosingPane);
};

/**
 * A full-featured component through which the user can view or edit a
 * single packet.
 *
 * Packet panes may be either docked within a main window
 * or may be floating freely in their own frames.
 */
class PacketPane : public QWidget, public regina::NPacketListener {
    Q_OBJECT

    private:
        /**
         * External components
         */
        ReginaMain* mainWindow;
        PacketWindow* frame;

        /**
         * Internal components
         */
        QLabel* headerIcon;
        QLabel* headerTitle;
        PacketUI* mainUI;
        QToolButton* dockUndockBtn;

        /**
         * Properties
         */
        bool dirty;
        bool readWrite;
        bool emergencyClosure;
        bool emergencyRefresh;
        bool isCommitting;

        /**
         * Internal actions
         */
        QAction* actCommit;
        QAction* actRefresh;
        QAction* actDockUndock;
        QAction* actClose;

        /**
         * Externally registered edit actions and their sources
         */
        QAction* editCut;
        QAction* editCopy;
        QAction* editPaste;

    public:
        /**
         * Constructs a new packet pane, managed by the given main window,
         * that views or edits the given packet.
         *
         * An appropriate internal interface component will be selected
         * by way of the PacketManager class.
         */
        PacketPane(ReginaMain* newMainWindow, regina::NPacket* newPacket,
            QWidget* parent = 0);
        ~PacketPane();

        /**
         * Query components and actions.
         */
        regina::NPacket* getPacket();
        ReginaMain* getMainWindow();
        PacketUI* getUI();

        /**
         * Set this pane to support or not support docking as appropriate.
         */
        void supportDock(bool shouldSupport);

        /**
         * Fill the given menu with the internal packet actions.
         * The menu must already be in a menu bar; otherwise Windows
         * platforms get horribly confused.
         * It is assumed that the given menu is empty.
         */
        void fillPacketTypeMenu(QMenu* menu);

        /**
         * Does this packet pane contain any changes that have not yet
         * been committed?
         */
        bool isDirty();

        /**
         * Signals that there are (or are not) changes in this interface
         * that have not yet been committed.  External interface
         * components will be updated accordingly.
         */
        void setDirty(bool newDirty);

        /**
         * Is this packet pane currently in read-write (as opposed to
         * read-only) mode?
         */
        bool isReadWrite();

        /**
         * Attempts to put this pane into read-write or read-only mode
         * as signalled by the \a allowReadWrite parameter.
         *
         * If \a allowReadWrite is \c true but nevertheless the pane
         * cannot be put into read-write mode, i.e., if
         * NPacket::isPacketEditable() returns \c false or the
         * underlying file is in read-only mode, then this routine will
         * do nothing and return \c false.  Otherwise this routine will
         * set the read-write status as requested and return \c true.
         */
        bool setReadWrite(bool allowReadWrite);

        /**
         * Are we allowed to close this packet pane?
         *
         * If this routine returns \c true, the caller of this routine
         * must ensure that the packet pane is actually closed (since in
         * this case queryClose() will call ReginaMain::isClosing()).
         */
        bool queryClose();

        /**
         * Registers or deregisters standard editor actions to operate
         * on this packet interface.
         *
         * Registered actions will be connected to appropriate edit
         * operations in this interface, and will be enabled and disabled
         * over time according to the current status of the internal UI
         * components.
         *
         * When the actions are deregistered, these relationships will be
         * broken and the actions will be left in a disabled state.
         *
         * Only one set of editor actions may be registered at a time.
         * Any attempt to register a new set of actions will
         * automatically deregister any previously registered actions.
         *
         * When a packet pane is destroyed, it is not a requirement that
         * currently registered actions be deregistered beforehand,
         * though the final enabled/disabled status of any remaining
         * actions that are not deregistered is not guaranteed.
         */
        void registerEditOperations(QAction* actCut, QAction* actCopy,
            QAction* actPaste);
        void deregisterEditOperations();

        /**
         * NPacketListener overrides.
         */
        void packetWasChanged(regina::NPacket* packet);
        void packetWasRenamed(regina::NPacket* packet);
        void packetToBeDestroyed(regina::NPacket* packet);
        void childWasAdded(regina::NPacket* packet, regina::NPacket* child);
        void childWasRemoved(regina::NPacket* packet, regina::NPacket* child,
            bool inParentDestructor);

    public slots:
        /**
         * Queries the packet and refreshes the interface accordingly.
         * Any uncommitted changes will be lost, though the user will be
         * prompted first.
         */
        void refresh();

        /**
         * Like refresh(), except that the user is never prompted.
         * Any uncommitted changes will be lost.
         */
        void refreshForce();

        /**
         * Commits any changes made in the user interface to the
         * underlying packet.
         *
         * If the commit should not be allowed (i.e., the packet
         * pane is in read-only mode or the underlying packet may not
         * be modified for mathematical reasons), it will not take
         * place.  Instead an appropriate error message will be displayed.
         *
         * Returns \c true if the commit took place (or if there were in
         * fact no changes to commit), or \c false if the commit was not
         * allowed to proceed.
         */
        bool commit();

        /**
         * Ensures that the underlying packet may be modified, and
         * commits any changes made in the user interface.
         *
         * This routine behaves identically to commit(), with one
         * exception.  Even if no changes have been made in the user
         * interface, an error will still be raised if changes to the
         * underlying packet are not allowed.
         *
         * This routine should be used when an operation plans to modify
         * the underlying packet, but wishes to commit any changes made
         * in the UI before proceeding.
         *
         * Returns \c true if changes are allowed to the underlying
         * packet (and therefore any necessary commit took place), or
         * \c false if changes are not allowed.
         */
        bool commitToModify();

        /**
         * Tries to commit any changes made in the user interface to the
         * underlying packet, but allows the user to continue if the
         * commit is not allowed.
         *
         * This routine behaves identically to commit(), except that if
         * the commit is not allowed then the error messages will be
         * replaced by gentle warnings explaining that an old copy of
         * the packet is to be used.  The user will be offered the
         * chance to continue or cancel.
         *
         * This routine should be used when an operation would like to
         * commit current changes made in the UI, but if the commit fails
         * then it can happily work with the last committed copy.
         *
         * Returns \c false if the user chose to cancel the operation,
         * or \c true if either the commit succeeded or the user was
         * happy to continue regardless.
         */
        bool tryCommit();

        /**
         * Closes this packet pane.  The user will be prompted if
         * necessary.
         *
         * For a packet pane that is currently docked, this routine
         * is equivalent to calling ReginaMain::closeDockedPane().
         *
         * Note that all this routine does is delegate the closure
         * operation to whatever component currently owns this packet pane.
         */
        bool close();

        /**
         * Closes this packet pane without prompting the user and
         * without the chance of the closure being cancelled.
         */
        void closeForce();

        /**
         * Docks this packet pane into the main window, if
         * it is not already docked.  If another packet pane is
         * already docked and refuses to be closed, the other pane will
         * be moved into its own freely floating window.
         *
         * This routine is the one and only way to dock a packet pane
         * that is currently floating in its own window.
         *
         * It is assumed that the packet pane is already registered with
         * the main window, and is either already docked or
         * currently floating in its own window.
         */
        void dockPane();

        /**
         * Floats this packet pane in its own top-level window, if it is
         * not already in such a state.
         *
         * This routine is the one and only way to construct a
         * top-level window enclosing a packet pane.
         *
         * It is assumed that the packet pane is already registered with
         * the main window, though it does not matter if the
         * pane is currently floating, docked or parentless.
         *
         * Note that a currently docked packet pane can also be floated by
         * calling ReginaMain::floatDockedPane(), which simply calls
         * this routine.
         */
        void floatPane();

        /**
         * Updates the enabled statuses of various registered editor
         * actions.  These slots are for internal use.
         */
        void updateClipboardActions();

    protected:
        /**
         * Refresh the packet label and icon in the header.
         */
        void refreshHeader();

        /**
         * Allow GUI updates from within a non-GUI thread.
         */
        void customEvent(QEvent* evt);
};

/**
 * A packet-specific interface for opening a packet using an external
 * viewer.
 */
typedef void (*PacketExternalViewer)(regina::NPacket* /* packet */,
    QWidget* /* parentWidget */);

inline PacketUI::PacketUI(PacketPane* newEnclosingPane) :
        enclosingPane(newEnclosingPane) {
}

inline PacketUI::~PacketUI() {
}

inline PacketEditIface* PacketUI::getEditIface() {
    return 0;
}

inline const QLinkedList<QAction*>& PacketUI::getPacketTypeActions() {
    return noActions;
}

inline void PacketUI::setDirty(bool newDirty) {
    enclosingPane->setDirty(newDirty);
}

inline PacketReadOnlyUI::PacketReadOnlyUI(PacketPane* newEnclosingPane) :
        PacketUI(newEnclosingPane) {
}

inline void PacketReadOnlyUI::commit() {
    setDirty(false);
}

inline void PacketReadOnlyUI::setReadWrite(bool) {
}

inline regina::NPacket* PacketPane::getPacket() {
    return mainUI->getPacket();
}

inline ReginaMain* PacketPane::getMainWindow() {
    return mainWindow;
}

inline PacketUI* PacketPane::getUI() {
    return mainUI;
}

inline bool PacketPane::isDirty() {
    return dirty;
}

inline bool PacketPane::isReadWrite() {
    return readWrite;
}

#endif
