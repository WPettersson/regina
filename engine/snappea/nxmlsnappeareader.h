
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

/*! \file snappea/nxmlsnappeareader.h
 *  \brief Deals with parsing XML data for SnapPea triangulations.
 */

#ifndef __NXMLSNAPPEAREADER_H
#ifndef __DOXYGEN
#define __NXMLSNAPPEAREADER_H
#endif

#include "regina-core.h"
#include "packet/nxmlpacketreader.h"
#include "snappea/nsnappeatriangulation.h"

namespace regina {

/**
 * \weakgroup snappea
 * @{
 */

/**
 * An XML packet reader that reads a single SnapPea triangulation.
 *
 * \ifacespython Not present.
 */
class REGINA_API NXMLSnapPeaReader : public NXMLPacketReader {
    private:
        NSnapPeaTriangulation* snappea_;
            /**< The SnapPea triangulation currently being read. */

    public:
        /**
         * Creates a new SnapPea triangulation reader.
         *
         * @param resolver the master resolver that will be used to fix
         * dangling packet references after the entire XML file has been read.
         */
        NXMLSnapPeaReader(NXMLTreeResolver& resolver);

        virtual NPacket* packet() override;
        virtual NXMLElementReader* startContentSubElement(
            const std::string& subTagName,
            const regina::xml::XMLPropertyDict& subTagProps) override;
        virtual void endContentSubElement(const std::string& subTagName,
            NXMLElementReader* subReader) override;
};

/*@}*/

// Inline functions for NXMLSnapPeaReader

inline NXMLSnapPeaReader::NXMLSnapPeaReader(NXMLTreeResolver& resolver) :
        NXMLPacketReader(resolver), snappea_(new NSnapPeaTriangulation()) {
}

inline NPacket* NXMLSnapPeaReader::packet() {
    return snappea_;
}

inline NXMLElementReader* NXMLSnapPeaReader::startContentSubElement(
        const std::string& subTagName, const regina::xml::XMLPropertyDict&) {
    if (subTagName == "snappea")
        return new NXMLCharsReader();
    else
        return new NXMLElementReader();
}

} // namespace regina

#endif

