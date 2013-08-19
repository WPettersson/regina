
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

/*! \file utilities/intutils.h
 *  \brief Miscellaneous utility classes for analysing the built-in
 *  integer types at compile time.
 */

#ifndef __INTUTILS_H
#ifndef __DOXYGEN
#define __INTUTILS_H
#endif

#include "regina-config.h"

namespace regina {

/**
 * Gives access to k-byte integer types, where \a k is a constant that
 * is not known until compile time.
 */
template <int bytes>
struct IntOfSize {
    /**
     * A signed integer type with \a k bytes, where \a k is the template
     * parameter.
     *
     * The default is \c void, which indicates that Regina does not know
     * how to access an integer type of the requested size.
     */
    typedef void type;
};

template <>
struct IntOfSize<1> {
    typedef __int8_t type;
};

template <>
struct IntOfSize<2> {
    typedef __int16_t type;
};

template <>
struct IntOfSize<4> {
    typedef __int32_t type;
};

template <>
struct IntOfSize<8> {
    typedef __int64_t type;
};

#ifdef INT128_FOUND
template <>
struct IntOfSize<16> {
    typedef __int128_t type;
};
#endif

} // namespace regina

#endif
