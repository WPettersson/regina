
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

/*! \file hypersurface/nxmlhypersurfacereader.h
 *  \brief Deals with parsing XML data for normal hypersurface lists.
 */

#ifndef __NXMLHYPERSURFACEREADER_H
#ifndef __DOXYGEN
#define __NXMLHYPERSURFACEREADER_H
#endif

#include "regina-core.h"
#include "hypersurface/nnormalhypersurfacelist.h"
#include "packet/nxmlpacketreader.h"

namespace regina {

/**
 * \weakgroup hypersurface
 * @{
 */

/**
 * An XML element reader that reads a single normal hypersurface in a
 * 4-manifold triangulation.
 *
 * \ifacespython Not present.
 */
class REGINA_API NXMLNormalHypersurfaceReader : public NXMLElementReader {
    private:
        NNormalHypersurface* surface_;
            /**< The normal hypersurface currently being read. */
        const Dim4Triangulation* tri_;
            /**< The triangulation in which this hypersurface lives. */
        HyperCoords coords_;
            /**< The coordinate system used by this hypersurface. */
        long vecLen_;
            /**< The length of corresponding normal hypersurface vector. */
        std::string name_;
            /**< The optional name associated with this normal hypersurface. */

    public:
        /**
         * Creates a new normal hypersurface reader.
         *
         * @param tri the triangulation in which this normal hypersurface lives.
         * @param coords the coordinate system used by this normal hypersurface.
         */
        NXMLNormalHypersurfaceReader(const Dim4Triangulation* tri,
            HyperCoords coords);

        /**
         * Returns the normal hypersurface that has been read.
         *
         * @return the newly allocated normal hypersurface, or 0 if an error
         * occurred.
         */
        NNormalHypersurface* hypersurface();

        virtual void startElement(const std::string& tagName,
            const regina::xml::XMLPropertyDict& tagProps,
            NXMLElementReader* parentReader);
        virtual void initialChars(const std::string& chars);
        virtual NXMLElementReader* startSubElement(
            const std::string& subTagName,
            const regina::xml::XMLPropertyDict& subTagProps);
};

/**
 * An XML packet reader that reads a single normal hypersurface list.
 *
 * \pre The parent XML element reader is in fact an
 * NXMLTriangulationReader.
 *
 * \ifacespython Not present.
 */
class REGINA_API NXMLNormalHypersurfaceListReader : public NXMLPacketReader {
    private:
        NNormalHypersurfaceList* list_;
            /**< The normal hypersurface list currently being read. */
        const Dim4Triangulation* tri_;
            /**< The triangulation in which these normal hypersurfaces live. */

    public:
        /**
         * Creates a new normal hypersurface list reader.
         *
         * @param tri the triangulation in which these normal hypersurfaces
         * live.
         * @param resolver the master resolver that will be used to fix
         * dangling packet references after the entire XML file has been read.
         */
        NXMLNormalHypersurfaceListReader(const Dim4Triangulation* tri,
            NXMLTreeResolver& resolver);

        virtual NPacket* packet() override;
        virtual NXMLElementReader* startContentSubElement(
            const std::string& subTagName,
            const regina::xml::XMLPropertyDict& subTagProps) override;
        virtual void endContentSubElement(const std::string& subTagName,
            NXMLElementReader* subReader) override;
};

/*@}*/

// Inline functions for NXMLNormalHypersurfaceReader

inline NXMLNormalHypersurfaceReader::NXMLNormalHypersurfaceReader(
        const Dim4Triangulation* tri, HyperCoords coords) :
        surface_(0), tri_(tri), coords_(coords), vecLen_(-1) {
}

inline NNormalHypersurface* NXMLNormalHypersurfaceReader::hypersurface() {
    return surface_;
}

// Inline functions for NXMLNormalHypersurfaceListReader

inline NXMLNormalHypersurfaceListReader::NXMLNormalHypersurfaceListReader(
        const Dim4Triangulation* tri, NXMLTreeResolver& resolver) :
        NXMLPacketReader(resolver), list_(0), tri_(tri) {
}

inline NPacket* NXMLNormalHypersurfaceListReader::packet() {
    return list_;
}

} // namespace regina

#endif

