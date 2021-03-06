
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
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

#include "dim4/dim4triangulation.h"
#include "packet/ncontainer.h"
#include "packet/nscript.h"
#include "packet/ntext.h"
#include "snappea/nsnappeatriangulation.h"

#include "newpacketdialog.h"
#include "packetcreator.h"
#include "packetfilter.h"
#include "packettreeview.h"
#include "reginamain.h"
#include "packettypes/dim2tricreator.h"
#include "packettypes/dim4tricreator.h"
#include "packettypes/nanglestructurecreator.h"
#include "packettypes/nhypersurfacecreator.h"
#include "packettypes/nnormalsurfacecreator.h"
#include "packettypes/nsurfacefiltercreator.h"
#include "packettypes/ntriangulationcreator.h"
#include "packettypes/nsnappeacreator.h"

void ReginaMain::newAngleStructures() {
    newPacket(new NAngleStructureCreator(),
        new SubclassFilter<regina::NTriangulation>(),
        tr("New Angle Structure Solutions"));
}

void ReginaMain::newContainer() {
    newPacket(new BasicPacketCreator<regina::NContainer>(), 0,
        tr("New Container"));
}

void ReginaMain::newDim2Triangulation() {
    newPacket(new Dim2TriangulationCreator(), 0,
        tr("New 2-Manifold Triangulation"));
}

void ReginaMain::newDim4Triangulation() {
    newPacket(new Dim4TriangulationCreator(), 0,
        tr("New 4-Manifold Triangulation"));
}

void ReginaMain::newFilter() {
    newPacket(new NSurfaceFilterCreator(), 0,
        tr("New Normal Surface Filter"));
}

void ReginaMain::newNormalSurfaces() {
    newPacket(new NNormalSurfaceCreator(),
        new SubclassFilter<regina::NTriangulation>(),
        tr("New Normal Surface List"));
}

void ReginaMain::newNormalHypersurfaces() {
    newPacket(new NHypersurfaceCreator(),
        new SubclassFilter<regina::Dim4Triangulation>(),
        tr("New Normal Hypersurface List"));
}

void ReginaMain::newPDF() {
    // A PDF must be an import; just use the import routine directly.
    importPDF();
}

void ReginaMain::newScript() {
    newPacket(new BasicPacketCreator<regina::NScript>(), 0,
        tr("New Script"));
}

void ReginaMain::newSnapPeaTriangulation() {
    newPacket(new NSnapPeaTriangulationCreator(this), 0,
        tr("New SnapPea Triangulation"));
}

void ReginaMain::newText() {
    newPacket(new BasicPacketCreator<regina::NText>(), 0,
        tr("New Text Packet"));
}

void ReginaMain::newTriangulation() {
    newPacket(new NTriangulationCreator(), 0,
        tr("New 3-Manifold Triangulation"));
}

void ReginaMain::newPacket(PacketCreator* creator, PacketFilter* parentFilter,
        const QString& dialogTitle) {
    NewPacketDialog dlg(this, creator, packetTree,
        treeView->selectedPacket(), parentFilter, dialogTitle);
    if (dlg.validate() && dlg.exec() == QDialog::Accepted) {
        regina::NPacket* newPacket = dlg.createdPacket();
        if (newPacket) {
            // Open a UI for the new packet, and select it in the tree.
            packetView(newPacket, true, true);
        }
    }
}

