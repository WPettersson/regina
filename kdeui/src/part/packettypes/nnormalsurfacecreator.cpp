
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2003, Ben Burton                                   *
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
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

// Regina core includes:
#include "surfaces/nnormalsurfacelist.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "coordinates.h"
#include "nnormalsurfacecreator.h"

#include <kcombobox.h>
#include <klocale.h>
#include <kmessagebox.h>
#include <qcheckbox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qwhatsthis.h>

using regina::NNormalSurfaceList;

namespace {
    /**
     * IDs corresponding to indices within the coordinate system combo box.
     */
    const int COORD_BOX_STANDARD = 0;
    const int COORD_BOX_AN_STANDARD = 1;
    const int COORD_BOX_QUAD = 2;
}

NNormalSurfaceCreator::NNormalSurfaceCreator() {
    // Set up the basic layout.
    ui = new QWidget();
    QBoxLayout* layout = new QVBoxLayout(ui);

    QBoxLayout* coordArea = new QHBoxLayout(layout, 5);
    QString expln = i18n("Specifies the coordinate system in which the "
        "vertex normal surfaces will be enumerated.");
    QLabel* label = new QLabel(i18n("Coordinate system:"), ui);
    QWhatsThis::add(label, expln);
    coordArea->addWidget(label);
    coords = new KComboBox(ui);
    QWhatsThis::add(coords, expln);
    coordArea->addWidget(coords, 1);

    layout->addSpacing(5);

    embedded = new QCheckBox(i18n("Embedded surfaces only"), ui);
    embedded->setChecked(true);
    QWhatsThis::add(embedded, i18n("Specifies whether only embedded "
        "normal surfaces should be enumerated, or whether all normal "
        "surfaces (embedded, immersed and singular) should be enumerated."));
    layout->addWidget(embedded);

    // Set up the individual coordinate systems.
    // Note that the order in which these options are added to the combo
    // box must correspond precisely to the coordinate system IDs defined
    // at the head of this file.
    coords->insertItem(Coordinates::name(NNormalSurfaceList::STANDARD));
    coords->insertItem(Coordinates::name(NNormalSurfaceList::AN_STANDARD));
    coords->insertItem(Coordinates::name(NNormalSurfaceList::QUAD));
    coords->setCurrentItem(0);
}

QWidget* NNormalSurfaceCreator::getInterface() {
    return ui;
}

regina::NPacket* NNormalSurfaceCreator::createPacket(regina::NPacket* parent,
        QWidget* parentWidget) {
    if (parent->getPacketType() != regina::NTriangulation::packetType) {
        KMessageBox::error(parentWidget, i18n(
            "Normal surface lists can only be created directly beneath "
            "triangulations."));
        return 0;
    }

    /**
     * Translate combo box indices into real coordinate system IDs.
     */
    int realCoordId;
    switch (coords->currentItem()) {
        case COORD_BOX_STANDARD:
            realCoordId = NNormalSurfaceList::STANDARD;
            break;
        case COORD_BOX_AN_STANDARD:
            realCoordId = NNormalSurfaceList::AN_STANDARD;
            break;
        case COORD_BOX_QUAD:
            realCoordId = NNormalSurfaceList::QUAD;
            break;
        default:
            KMessageBox::error(parentWidget, i18n("No coordinate system "
                "has been selected in which to create the new normal "
                "surfaces."));
            return 0;
    }

    return NNormalSurfaceList::enumerate(
        dynamic_cast<regina::NTriangulation*>(parent), realCoordId,
        embedded->isChecked());
}

