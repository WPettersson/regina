
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

#include <cstdlib>
#include <cctype>
#include "maths/nlargeinteger.h"
#include "utilities/stringutils.h"
#include "utilities/nbooleans.h"

namespace regina {

char* duplicate(const std::string& str) {
    char* ans = new char[str.length() + 1];

    char* pos = ans;
    for (std::string::const_iterator it = str.begin(); it != str.end(); it++)
        *(pos++) = *it;
    *pos = 0;

    return ans;
}

bool startsWith(const std::string& str, const std::string& prefix) {
    if (str.length() < prefix.length())
        return false;
    return (str.substr(0, prefix.length()).compare(prefix) == 0);
}

std::string stripWhitespace(const std::string& str) {
    std::string::size_type start = 0;
    std::string::size_type end = str.length();

    while (start < end && isspace(str[start]))
        start++;
    while (start < end && isspace(str[end - 1]))
        end--;

    return str.substr(start, end - start);
}

bool valueOf(const std::string& str, int& dest) {
    // TODO: Check errno, and check for overflow when casting back to int.
    char* endPtr;
    dest = static_cast<int>(strtol(str.c_str(), &endPtr, 10));
    return ((! str.empty()) && (*endPtr == 0));
}

bool valueOf(const std::string& str, unsigned& dest) {
    // TODO: Check errno, and check for overflow when casting back to unsigned.
    char* endPtr;
    dest = static_cast<unsigned>(strtoul(str.c_str(), &endPtr, 10));
    return ((! str.empty()) && (*endPtr == 0));
}

bool valueOf(const std::string& str, long& dest) {
    char* endPtr;
    dest = strtol(str.c_str(), &endPtr, 10);
    return ((! str.empty()) && (*endPtr == 0));
}

bool valueOf(const std::string& str, unsigned long& dest) {
    char* endPtr;
    dest = strtoul(str.c_str(), &endPtr, 10);
    return ((! str.empty()) && (*endPtr == 0));
}

bool valueOf(const std::string& str, NLargeInteger& dest) {
    bool valid;
    dest = NLargeInteger(str.c_str(), 10, &valid);
    return valid;
}

bool valueOf(const std::string& str, double& dest) {
    char* endPtr;
    dest = strtod(str.c_str(), &endPtr);
    return ((! str.empty()) && (*endPtr == 0));
}

bool valueOf(const std::string& str, bool& dest) {
    if (str.empty()) {
        dest = false;
        return false;
    }
    if (str[0] == 't' || str[0] == 'T') {
        dest = true;
        return true;
    }
    dest = false;
    return (str[0] == 'F' || str[0] == 'f');
}

bool valueOf(const std::string& str, NTriBool& dest) {
    if (str.empty()) {
        dest.setUnknown();
        return false;
    }
    if (str[0] == 't' || str[0] == 'T' || str == "1") {
        dest = true;
        return true;
    }
    if (str[0] == 'f' || str[0] == 'F' || str == "-1") {
        dest = false;
        return true;
    }
    dest.setUnknown();
    return (str[0] == 'u' || str[0] == 'U' || str == "0");
}

bool valueOf(const std::string& str, NBoolSet& dest) {
    if (str.length() != 2) {
        dest = NBoolSet::sNone;
        return false;
    }
    char t = str[0];
    char f = str[1];
    if (t != '-' && t != 'T' && t != 't') {
        dest = NBoolSet::sNone;
        return false;
    }
    if (f != '-' && f != 'F' && f != 'f') {
        dest = NBoolSet::sNone;
        return false;
    }

    dest = NBoolSet(t != '-', f != '-');
    return true;
}

} // namespace regina
