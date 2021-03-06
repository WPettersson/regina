
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

/*! \file angle/nxmlanglestructreader.h
 *  \brief Deals with parsing XML data for angle structure lists.
 */

#ifndef __NXMLANGLESTRUCTREADER_H
#ifndef __DOXYGEN
#define __NXMLANGLESTRUCTREADER_H
#endif

#include "regina-core.h"
#include "packet/nxmlpacketreader.h"
#include "angle/nanglestructurelist.h"

namespace regina {

/**
 * \weakgroup angle
 * @{
 */

/**
 * An XML element reader that reads a single angle structure.
 *
 * \ifacespython Not present.
 */
class REGINA_API NXMLAngleStructureReader : public NXMLElementReader {
    private:
        NAngleStructure* angles;
            /**< The angle structure currently being read. */
        NTriangulation* tri;
            /**< The triangulation on which this angle structure is placed. */
        long vecLen;
            /**< The length of corresponding angle structure vector. */

    public:
        /**
         * Creates a new angle structure reader.
         *
         * @param newTri the triangulation on which this angle structure lies.
         */
        NXMLAngleStructureReader(NTriangulation* newTri);

        /**
         * Returns the angle structure that has been read.
         *
         * @return the newly allocated angle structure, or 0 if an error
         * occurred.
         */
        NAngleStructure* structure();

        virtual void startElement(const std::string& tagName,
            const regina::xml::XMLPropertyDict& tagProps,
            NXMLElementReader* parentReader);
        virtual void initialChars(const std::string& chars);
        virtual NXMLElementReader* startSubElement(
            const std::string& subTagName,
            const regina::xml::XMLPropertyDict& subTagProps);
};

/**
 * An XML packet reader that reads a single angle structure list.
 *
 * \pre The parent XML element reader is in fact an
 * NXMLTriangulationReader.
 *
 * \ifacespython Not present.
 */
class REGINA_API NXMLAngleStructureListReader : public NXMLPacketReader {
    private:
        NAngleStructureList* list;
            /**< The angle structure list currently being read. */
        NTriangulation* tri;
            /**< The triangulation on which these angle structures
                 are placed. */

    public:
        /**
         * Creates a new angle structure list reader.
         *
         * @param newTri the triangulation on which these angle
         * structures are placed.
         * @param resolver the master resolver that will be used to fix
         * dangling packet references after the entire XML file has been read.
         */
        NXMLAngleStructureListReader(NTriangulation* newTri,
            NXMLTreeResolver& resolver);

        virtual NPacket* packet() override;
        virtual NXMLElementReader* startContentSubElement(
            const std::string& subTagName,
            const regina::xml::XMLPropertyDict& subTagProps) override;
        virtual void endContentSubElement(const std::string& subTagName,
            NXMLElementReader* subReader) override;
};

/*@}*/

// Inline functions for NXMLAngleStructureReader

inline NXMLAngleStructureReader::NXMLAngleStructureReader(
        NTriangulation* newTri) : angles(0), tri(newTri), vecLen(-1) {
}

inline NAngleStructure* NXMLAngleStructureReader::structure() {
    return angles;
}

// Inline functions for NXMLAngleStructureListReader

inline NXMLAngleStructureListReader::NXMLAngleStructureListReader(
        NTriangulation* newTri, NXMLTreeResolver& resolver) :
        NXMLPacketReader(resolver),
        list(new NAngleStructureList(false)), tri(newTri) {
}

inline NPacket* NXMLAngleStructureListReader::packet() {
    return list;
}

} // namespace regina

#endif

