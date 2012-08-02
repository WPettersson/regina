
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

/*! \file ntricomposition.h
 *  \brief Provides a combinatorial composition viewer for triangulations.
 */

#ifndef __NTRICOMPOSITION_H
#define __NTRICOMPOSITION_H

#include "packet/npacketlistener.h"

#include "../packettabui.h"

#include <memory>

class PacketChooser;
class PacketEditIface;
class QMenu;
class QPushButton;
class QTreeWidget;
class QTreeWidgetItem;

namespace regina {
    class NIsomorphism;
    class NMatrix2;
    class NPacket;
    class NPerm4;
    class NSatRegion;
    class NStandardTriangulation;
    class NTriangulation;
};

/**
 * A triangulation page for viewing the combinatorial composition.
 */
class NTriCompositionUI : public QObject, public PacketViewerTab,
        public regina::NPacketListener {
    Q_OBJECT

    private:
        /**
         * Describes the type of isomorphism relationship that has been
         * discovered, if any.
         */
        enum IsomorphismType
            { NoRelationship, IsIsomorphic, IsSubcomplex, IsSupercomplex };

        /**
         * Packet details
         */
        regina::NTriangulation* tri;
        regina::NTriangulation* comparingTri;
        std::auto_ptr<regina::NIsomorphism> isomorphism;
        IsomorphismType isoType;

        /**
         * Internal components
         */
        QWidget* ui;
        PacketChooser* isoTest;
        QLabel* isoResult;
        QPushButton* isoView;
        QTreeWidget* details;
        QTreeWidgetItem* components;
        QTreeWidgetItem* lastComponent;
        PacketEditIface* editIface;

    public:
        /**
         * Constructor.
         */
        NTriCompositionUI(regina::NTriangulation* packet,
                PacketTabbedUI* useParentUI);
        ~NTriCompositionUI();

        /**
         * PacketViewerTab overrides.
         */
        regina::NPacket* getPacket();
        QWidget* getInterface();
        PacketEditIface* getEditIface();
        void refresh();
        void editingElsewhere();

        /**
         * NPacketListener overrides.
         */
        void packetToBeDestroyed(regina::NPacket* packet);

    public slots:
        /**
         * Update the isomorphism test panel.
         */
        void updateIsoPanel();

        /**
         * View the isomorphism details.
         */
        void viewIsomorphism();

    private:
        /**
         * Add new items to the list view.
         */
        QTreeWidgetItem* addTopLevelSection(const QString& text);
        QTreeWidgetItem* addComponentSection(const QString& text);

        /**
         * Fill the list view with information.
         */
        void findAugTriSolidTori();
        void findBlockedTriangulations();
        void findL31Pillows();
        void findLayeredChainPairs();
        void findLayeredLensSpaces();
        void findLayeredLoops();
        void findLayeredSolidTori();
        void findPillowSpheres();
        void findPlugTriSolidTori();
        void findSnappedBalls();
        void findSnappedSpheres();
        void findSpiralSolidTori();
        void describeSatRegion(const regina::NSatRegion& region,
            QTreeWidgetItem* parent);

        /**
         * Return string representations of various items.
         */
        static QString edgeString(unsigned long tetIndex, int edge1,
            int edge2);
        static QString edgeString(unsigned long tetIndex,
            const regina::NPerm4& roles, int startPreimage, int endPreimage);
        static QString matrixString(const regina::NMatrix2& matrix);
};

inline PacketEditIface* NTriCompositionUI::getEditIface() {
    return editIface;
}

#endif
