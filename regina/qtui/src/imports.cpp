
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Qt User Interface                                                    *
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

#include "packet/npacket.h"

#include "packettreeview.h"
#include "reginamain.h"
#include "foreign/dehydrationhandler.h"
#include "foreign/importdialog.h"
#include "foreign/isosighandler.h"
#include "foreign/orbhandler.h"
#include "foreign/pdfhandler.h"
#include "foreign/pythonhandler.h"
#include "foreign/reginahandler.h"
#include "foreign/snappeahandler.h"
#include "reginafilter.h"

#include <QFileDialog>

void ReginaMain::importDehydration() {
    importFile(DehydrationHandler::instance, 0, tr(FILTER_ALL),
        tr("Import Dehydrated Triangulation List"));
}

void ReginaMain::importIsoSig3() {
    importFile(IsoSigHandler::instance3, 0, tr(FILTER_ALL),
        tr("Import Isomorphism Signature List"));
}

void ReginaMain::importPDF() {
    importFile(PDFHandler::instance, 0, tr(FILTER_PDF),
        tr("Import PDF Document"));
}

void ReginaMain::importPython() {
    importFile(PythonHandler::instance, 0, tr(FILTER_PYTHON_SCRIPTS),
        tr("Import Python Script"));
}

void ReginaMain::importRegina() {
    importFile(ReginaHandler(), 0, tr(FILTER_REGINA),
        tr("Import Regina Data File"));
}

void ReginaMain::importSnapPea() {
    importFile(SnapPeaHandler::instance, 0, tr(FILTER_SNAPPEA),
        tr("Import SnapPea Triangulation"));
}

void ReginaMain::importOrb() {
    importFile(OrbHandler::instance, 0, tr(FILTER_ORB),
        tr("Import Orb or Casson Triangulation"));
}

void ReginaMain::importFile(const PacketImporter& importer,
        PacketFilter* parentFilter, const QString& fileFilter,
        const QString& dialogTitle) {
    QString file = QFileDialog::getOpenFileName(this,
        dialogTitle, QString(), fileFilter);
    if (file.isEmpty())
        return;
    regina::NPacket* newTree = importer.importData(file, this);

    if (newTree) {
        ImportDialog dlg(this, newTree, packetTree,
            treeView->selectedPacket(), parentFilter,
            importer.useImportEncoding(), dialogTitle);
        if (dlg.validate() && dlg.exec() == QDialog::Accepted)
            packetView(newTree, true);
        else
            delete newTree;
    }
}

