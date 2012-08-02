
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

/*! \file packetimporter.h
 *  \brief Provides a basic infrastructure for importing packet trees
 *  from foreign file formats.
 */

#ifndef __PACKETIMPORTER_H
#define __PACKETIMPORTER_H

class QString;
class QTextCodec;
class QWidget;

namespace regina {
    class NPacket;
};

/**
 * An object responsible for importing a packet tree from a foreign file
 * format.  Different foreign file formats should correspond to different
 * subclasses of PacketImporter.
 */
class PacketImporter {
    public:
        /**
         * Default destructor that does nothing.
         */
        virtual ~PacketImporter();

        /**
         * Import a packet tree from the given file.
         * If a text encoding is required, the preferred codec
         * ReginaPrefSet::importExportCodec() should be used.
         *
         * If the import is unsuccessful, this routine should display
         * an appropriate error to the user (using the argument
         * \a parentWidget as a parent for the message box) and return 0.
         * Otherwise the imported packet tree should be returned.
         *
         * Sensible packet labels must be assigned to all packets in the
         * imported tree.
         */
        virtual regina::NPacket* importData(const QString& fileName,
            QWidget* parentWidget) const = 0;

        /**
         * Should the GUI inform the user that their preferred codec
         * will be used?
         *
         * The default implementation returns \c false.
         */
        virtual bool useImportEncoding() const;
};

inline PacketImporter::~PacketImporter() {
}

inline bool PacketImporter::useImportEncoding() const {
    return false;
}

#endif
