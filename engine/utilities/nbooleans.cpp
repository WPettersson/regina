
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

#include "utilities/nbooleans.h"

#include <iostream>

namespace regina {

const unsigned char NBoolSet::eltTrue = 1;
const unsigned char NBoolSet::eltFalse = 2;
const NBoolSet NBoolSet::sNone;
const NBoolSet NBoolSet::sTrue(true);
const NBoolSet NBoolSet::sFalse(false);
const NBoolSet NBoolSet::sBoth(true, true);

std::ostream& operator << (std::ostream& out, const NBoolSet& set) {
    if (set == NBoolSet::sNone)
        out << "{ }";
    else if (set == NBoolSet::sTrue)
        out << "{ true }";
    else if (set == NBoolSet::sFalse)
        out << "{ false }";
    else
        out << "{ true, false }";
    return out;
}

} // namespace regina

