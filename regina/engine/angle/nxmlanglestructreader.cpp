
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include <vector>
#include "angle/nxmlanglestructreader.h"
#include "triangulation/ntriangulation.h"
#include "utilities/stringutils.h"

namespace regina {

void NXMLAngleStructureReader::startElement(const std::string&,
        const regina::xml::XMLPropertyDict& props,
        NXMLElementReader*) {
    if (! valueOf(props.lookup("len"), vecLen))
        vecLen = -1;
}

void NXMLAngleStructureReader::initialChars(const std::string& chars) {
    if (vecLen < 0 || tri == 0)
        return;

    std::vector<std::string> tokens;
    if (basicTokenise(back_inserter(tokens), chars) % 2 != 0)
        return;

    // Create a new vector and read all non-zero entries.
    NAngleStructureVector* vec = new NAngleStructureVector(vecLen);

    long pos;
    NLargeInteger value;
    for (unsigned long i = 0; i < tokens.size(); i += 2) {
        if (valueOf(tokens[i], pos))
            if (valueOf(tokens[i + 1], value))
                if (pos >= 0 && pos < vecLen) {
                    // All looks valid.
                    vec->setElement(pos, value);
                    continue;
                }

        // Found something invalid.
        delete vec;
        return;
    }

    angles = new NAngleStructure(tri, vec);
}

NXMLElementReader* NXMLAngleStructureReader::startSubElement(
        const std::string& subTagName,
        const regina::xml::XMLPropertyDict& props) {
    if (! angles)
        return new NXMLElementReader();

    /** Flags in data files are deprecated as of Regina 4.93.
    if (subTagName == "flags") {
        if (! valueOf(props.lookup("value"), angles->flags))
            angles->flags = 0;
    }
    */
    return new NXMLElementReader();
}

NXMLElementReader* NXMLAngleStructureListReader::startContentSubElement(
        const std::string& subTagName,
        const regina::xml::XMLPropertyDict& props) {
    bool b;
    if (subTagName == "angleparams") {
        if (valueOf(props.lookup("tautonly"), b))
            list->tautOnly_ = b;
    } else if (subTagName == "struct") {
        return new NXMLAngleStructureReader(tri);
    } else if (subTagName == "spanstrict") {
        if (valueOf(props.lookup("value"), b))
            list->doesSpanStrict = b;
    } else if (subTagName == "spantaut") {
        if (valueOf(props.lookup("value"), b))
            list->doesSpanTaut = b;
    } else if (subTagName == "allowstrict") {
        if (valueOf(props.lookup("value"), b))
            list->doesSpanStrict = b;
    } else if (subTagName == "allowtaut") {
        if (valueOf(props.lookup("value"), b))
            list->doesSpanTaut = b;
    }
    return new NXMLElementReader();
}

void NXMLAngleStructureListReader::endContentSubElement(
        const std::string& subTagName,
        NXMLElementReader* subReader) {
    if (subTagName == "struct")
        if (NAngleStructure* s =
                dynamic_cast<NXMLAngleStructureReader*>(subReader)->
                getStructure())
            list->structures.push_back(s);
}

NXMLPacketReader* NAngleStructureList::getXMLReader(NPacket* parent) {
    return new NXMLAngleStructureListReader(
        dynamic_cast<NTriangulation*>(parent));
}

} // namespace regina

