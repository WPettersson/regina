
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

#include <cstdio>
#include <cstring>
#include <sstream>
#include "file/nfile.h"
#include "file/nfileinfo.h"
#include "utilities/zstream.h"

namespace regina {

#define STARTS_FALSE 0
#define STARTS_TRUE 1
#define STARTS_COULD_NOT_OPEN 2

const int NFileInfo::TYPE_BINARY = 1;
const int NFileInfo::TYPE_XML = 2;

namespace {
    /**
     * Does the given file begin with the given set of characters?
     *
     * Returns STARTS_FALSE, STARTS_TRUE or STARTS_COULD_NOT_OPEN.
     */
    int fileStartsWith(const char* file, const char* prefix) {
        FILE* f = fopen(file, "rb");
        if (! f)
            return STARTS_COULD_NOT_OPEN;

        int ans = STARTS_FALSE;
        size_t len = strlen(prefix);
        char* buf = new char[len];

        if (fread(buf, sizeof(char), len, f) == len)
            if (strncmp(buf, prefix, len) == 0)
                ans = STARTS_TRUE;

        delete[] buf;
        fclose(f);
        return ans;
    }
}

NFileInfo* NFileInfo::identify(const std::string& idPathname) {
    // Check for an old-style binary file.
    int starts = fileStartsWith(idPathname.c_str(), "Normal");
    if (starts == STARTS_FALSE)
        starts = fileStartsWith(idPathname.c_str(), "Regina");

    if (starts == STARTS_COULD_NOT_OPEN)
        return 0;
    else if (starts == STARTS_TRUE) {
        NFileInfo* ans = new NFileInfo();
        ans->pathname = idPathname;
        ans->type = NFileInfo::TYPE_BINARY;
        ans->typeDescription = "Binary Regina data file (obsolete format)";
        ans->compressed = false;

        regina::NFile file;
        if (! file.open(idPathname.c_str(),
                regina::NRandomAccessResource::READ)) {
            ans->invalid = true;
        } else {
            std::ostringstream out;
            out << file.getMajorVersion() << '.' << file.getMinorVersion();
            ans->engine = out.str();
            ans->invalid = false;
        }
        file.close();
        return ans;
    }

    // Check for an XML file.
    NFileInfo* ans = 0;
    if (fileStartsWith(idPathname.c_str(), "<?xml") == STARTS_TRUE) {
        ans = new NFileInfo();
        ans->compressed = false;
    } else {
        regina::DecompressionStream in(idPathname.c_str());
        if (in) {
            std::string s;
            in >> s;
            if ((! in.eof()) && (s == "<?xml")) {
                ans = new NFileInfo();
                ans->compressed = true;
            }
        }
    }

    if (ans) {
        ans->pathname = idPathname;
        ans->type = NFileInfo::TYPE_XML;
        ans->typeDescription = "XML Regina data file";

        // Make it an invalid file until we know otherwise.
        ans->invalid = true;

        regina::DecompressionStream in(idPathname.c_str());
        if (! in)
            return ans;

        std::string s;

        // Start by slurping in the opening "<?xml".
        if (in.eof())
            return ans;
        in >> s;
        if (s != "<?xml")
            return ans;

        // Hunt for the matching "...?>".
        // Try skipping through several strings in case there are extra
        // arguments in the XML prologue (such as encoding or standalone
        // declarations).
        int i;
        for (i = 0; ; i++) {
            if (in.eof())
                return ans;
            in >> s;
            if (s.length() >= 2 &&
                    s[s.length() - 2] == '?' &&
                    s[s.length() - 1] == '>')
                break;

            // If we can't find it after enough tries, just give up.
            // Ten tries should be more than sufficient, since the current XML
            // spec supports only version, encoding and standalone arguments
            // at present.
            if (i >= 10)
                return ans;
        }

        // The next thing we see should be the <reginadata ...> element.
        if (in.eof())
            return ans;
        in >> s;
        if (s != "<reginadata")
            return ans;

        // Next should be the engine version.
        if (in.eof())
            return ans;
        in >> s;
        if (s.length() < 8)
            return ans;
        if (s.substr(0, 8).compare("engine=\"") != 0)
            return ans;

        // We've found the engine attribute; extract its value.
        std::string::size_type pos = s.find('"', 8);
        if (pos == std::string::npos)
            return ans;
        ans->engine = s.substr(8, pos - 8);

        // That's as far as we need to go; we've extracted everything we want.
        ans->invalid = false;
        return ans;
    }

    // Unknown format.
    return 0;
}

void NFileInfo::writeTextShort(std::ostream& out) const {
    out << "File information: " << typeDescription;
    if (compressed)
        out << " (compressed)";
}

void NFileInfo::writeTextLong(std::ostream& out) const {
    out << "Regina data\n" << typeDescription;
    if (compressed)
        out << " (compressed)";
    out << '\n';

    if (invalid)
        out << "File contains invalid metadata.\n";
    else
        out << "Engine " << engine << '\n';
}

} // namespace regina

