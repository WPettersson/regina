
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

#include "foreign/csvsurfacelist.h"
#include "surfaces/nnormalsurfacelist.h"

#include "csvsurfacehandler.h"
#include "reginasupport.h"
#include "../packetfilter.h"

#include <QFile>
#include <QTextDocument>

const CSVSurfaceHandler CSVSurfaceHandler::instance;

PacketFilter* CSVSurfaceHandler::canExport() const {
    return new SingleTypeFilter<regina::NNormalSurfaceList>();
}

bool CSVSurfaceHandler::exportData(regina::NPacket* data,
        const QString& fileName, QWidget* parentWidget) const {
    regina::NNormalSurfaceList* list =
        dynamic_cast<regina::NNormalSurfaceList*>(data);
    if (! regina::writeCSVStandard(
            static_cast<const char*>(QFile::encodeName(fileName)), *list)) {
        ReginaSupport::warn(parentWidget,
            QObject::tr("The export failed."), 
            QObject::tr("<qt>An unknown error occurred, probably related "
            "to file I/O.  Please check that you have permissions to write "
            "to the file <tt>%1</tt>.</qt>").arg(Qt::escape(fileName)));
        return false;
    }
    return true;
}
