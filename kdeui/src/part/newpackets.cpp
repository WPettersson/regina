
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

#include "packet/ncontainer.h"
#include "packet/nscript.h"
#include "packet/ntext.h"
#include "triangulation/ntriangulation.h"

#include "newpacketdialog.h"
#include "packetcreator.h"
#include "packetfilter.h"
#include "packettreeview.h"
#include "reginapart.h"
#include "packettypes/nanglestructurecreator.h"
#include "packettypes/nsurfacefiltercreator.h"

#include <klocale.h>

void ReginaPart::newAngleStructures() {
    newPacket(new NAngleStructureCreator(),
        new SingleTypeFilter<regina::NTriangulation>(),
        i18n("New Angle Structure Solutions"), i18n("Angle Structures"));
}

void ReginaPart::newCensus() {
    unimplemented();
}

void ReginaPart::newContainer() {
    newPacket(new BasicPacketCreator<regina::NContainer>(), 0,
        i18n("New Container"), i18n("Container"));
}

void ReginaPart::newFilter() {
    newPacket(new NSurfaceFilterCreator(), 0,
        i18n("New Normal Surface Filter"), i18n("Surface Filter"));
}

void ReginaPart::newNormalSurfaces() {
    unimplemented();
}

void ReginaPart::newScript() {
    newPacket(new BasicPacketCreator<regina::NScript>(), 0,
        i18n("New Script"), i18n("Script"));
}

void ReginaPart::newText() {
    newPacket(new BasicPacketCreator<regina::NText>(), 0,
        i18n("New Text Packet"), i18n("Text"));
}

void ReginaPart::newTriangulation() {
    unimplemented();
}

void ReginaPart::newPacket(PacketCreator* creator, PacketFilter* parentFilter,
        const QString& dialogTitle, const QString& suggestedLabel) {
    if (! checkReadWrite())
        return;

    NewPacketDialog dlg(widget(), creator, packetTree,
        treeView->selectedPacket(), parentFilter, dialogTitle, suggestedLabel);
    if (dlg.exec() == QDialog::Accepted) {
        regina::NPacket* newPacket = dlg.createdPacket();
        if (newPacket) {
            QListViewItem* item = treeView->find(newPacket);
            if (item)
                treeView->ensureItemVisible(item);
            packetView(newPacket);

            setModified(true);
        }
    }
}

