
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

/*! \file surfaces/ndisctype.h
 *  \brief Deals with normal and almost normal disc types.
 */

#ifndef __NDISCTYPE_H
#ifndef __DOXYGEN
#define __NDISCTYPE_H
#endif

#include <iostream>
#include "regina-core.h"

namespace regina {

/**
 * \weakgroup surfaces
 * @{
 */

/**
 * Identifies a single normal or almost normal disc type within a
 * triangulation.
 *
 * A disc type is identified by a tetrahedron index (the data member
 * \a tetIndex), and a disc type within that tetrahedron (the data
 * member \a type).  The latter could mean any number of things according
 * to the application at hand.  For instance, if we are tracking quad types
 * then \a type might be an integer between 0 and 2 inclusive, or if we
 * are tracking all normal discs in standard coordinates then \a type
 * might be an integer between 0 and 6 inclusive.  Ultimately, the
 * specific meaning of \a type is left to the user.
 *
 * It is however assumed that \a type will always be non-negative for
 * "meaningful" disc types; this is to ensure that the constant NONE
 * does not clash with any meaningful values.
 *
 * Note that this class tracks disc \a types, not discs themselves.
 * To track individual normal discs, see the NDiscSpec class instead.
 */
struct REGINA_API NDiscType {
    static const NDiscType NONE;
        /**< Represents a "null" disc type.  Here the \a type member is
             negative, to distinguish it from "meaningful" disc types
             in which \a type is always zero or positive. */

    unsigned long tetIndex;
        /**< The index within the triangulation of the tetrahedron
             containing this disc type.  This must be between 0 and
             NTriangulation::getNumberOfTetrahedra()-1 inclusive. */
    int type;
        /**< Identifies the disc type within the specified tetrahedron.
             The precise meaning of this member is left up to the user,
             though it must be non-negative for "meaningful" disc types.
             See the NDiscType class notes for details. */

    /**
     * Creates a new disc type initialised to NONE.
     */
    NDiscType();
    /**
     * Creates a new disc type initialised with the given values.
     *
     * @param newTetIndex the index within the triangulation of the
     * tetrahedron containing this disc type.
     * @param newType the specific disc type within the given
     * tetrahedron; see the class notes for the meaning of this field.
     */
    NDiscType(unsigned long newTetIndex, int newType);
    /**
     * Creates a copy of the given disc type.
     *
     * @param cloneMe the disc type to clone.
     */
    NDiscType(const NDiscType& cloneMe);

    /**
     * Sets this to a copy of the given disc type.
     *
     * @param cloneMe the disc type to clone.
     * @return a reference to this disc type.
     */
    NDiscType& operator = (const NDiscType& cloneMe);
    /**
     * Determines if this and the given disc type are identical.
     *
     * Note that NONE is considered identical to NONE, and that NONE will
     * not be equal to any "meaningful" disc type (specifically, a disc type
     * for which \a type is non-negative).
     *
     * @return \c true if this and the given disc type are identical, or
     * \c false if they are different.
     */
    bool operator == (const NDiscType& compare) const;
    /**
     * Determines if this and the given disc type are different.
     *
     * This is the negation of the equality test; see operator == for
     * further details.
     *
     * @return \c true if this and the given disc type are different, or
     * \c false if they are identical.
     */
    bool operator != (const NDiscType& compare) const;
    /**
     * Provides an ordering of disc types.  Types are ordered first by
     * \a tetrahedron and then by \a type.  NONE is considered less than
     * all "meaningful" disc types.
     *
     * @return \c true if this disc type appears before the given disc type
     * in the ordering, or \c false if not.
     */
    bool operator < (const NDiscType& compare) const;
};

/**
 * Writes the given disc type to the given output stream.
 * The disc type will be written as a pair <tt>(tetIndex, type)</tt>.
 *
 * @param out the output stream to which to write.
 * @param type the disc type to write.
 * @return a reference to the given output stream.
 */
REGINA_API std::ostream& operator << (std::ostream& out, const NDiscType& type);

/*@}*/

// Inline functions for NDiscType

inline NDiscType::NDiscType() : tetIndex(0), type(-1) {
}

inline NDiscType::NDiscType(unsigned long newTetIndex, int newType) :
        tetIndex(newTetIndex), type(newType) {
}

inline NDiscType::NDiscType(const NDiscType& cloneMe) :
        tetIndex(cloneMe.tetIndex), type(cloneMe.type) {
}

inline NDiscType& NDiscType::operator = (const NDiscType& cloneMe) {
    tetIndex = cloneMe.tetIndex;
    type = cloneMe.type;
    return *this;
}

inline bool NDiscType::operator == (const NDiscType& compare) const {
    return (tetIndex == compare.tetIndex && type == compare.type);
}

inline bool NDiscType::operator != (const NDiscType& compare) const {
    return (tetIndex != compare.tetIndex || type != compare.type);
}

inline bool NDiscType::operator < (const NDiscType& compare) const {
    return (tetIndex < compare.tetIndex ||
        (tetIndex == compare.tetIndex && type < compare.type));
}

inline std::ostream& operator << (std::ostream& out, const NDiscType& type) {
    return out << '(' << type.tetIndex << ", " << type.type << ')';
}

} // namespace regina

#endif

