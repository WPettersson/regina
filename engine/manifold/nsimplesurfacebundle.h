
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

/*! \file manifold/nsimplesurfacebundle.h
 *  \brief Deals with simple closed surface bundles.
 */

#ifndef __NSIMPLESURFACEBUNDLE_H
#ifndef __DOXYGEN
#define __NSIMPLESURFACEBUNDLE_H
#endif

#include "regina-core.h"
#include "nmanifold.h"

namespace regina {

/**
 * \weakgroup manifold
 * @{
 */

/**
 * Represents a particularly simple closed surface bundle over the circle.
 * Only 2-sphere bundles, twisted 2-sphere bundles and projective plane
 * bundles are considered.
 *
 * All optional NManifold routines are implemented for this class.
 */
class REGINA_API NSimpleSurfaceBundle : public NManifold {
    public:
        /**
         * Represents the orientable 2-sphere bundle over the circle.
         */
        static const int S2xS1;
        /**
         * Represents the non-orientable twisted 2-sphere bundle over the
         * circle.
         */
        static const int S2xS1_TWISTED;
        /**
         * Represents the projective plane bundle over the circle.
         */
        static const int RP2xS1;
    private:
        int type_;
            /**< The specific surface bundle being represented.
                 This must be one of the 3-manifold constants defined in
                 this class. */

    public:
        /**
         * Creates a new surface bundle of the given type.
         *
         * @param newType the specific type of surface bundle to
         * represent.  This must be one of the 3-manifold constants
         * defined in this class.
         */
        NSimpleSurfaceBundle(int newType);
        /**
         * Creates a clone of the given surface bundle.
         *
         * @param cloneMe the surface bundle to clone.
         */
        NSimpleSurfaceBundle(const NSimpleSurfaceBundle& cloneMe);
        /**
         * Returns the specific type of surface bundle being represented.
         *
         * @return the type of surface bundle.  This will be one of the
         * 3-manifold constants defined in this class.
         */
        int type() const;
        /**
         * Deprecated routine that returns the specific type of surface bundle
         * being represented.
         *
         *
         * \deprecated This routine has been renamed to type().
         * See the type() documentation for further details.
         */
        REGINA_DEPRECATED int getType() const;
        /**
         * Determines whether this and the given surface bundle represent
         * the same 3-manifold.
         *
         * @param compare the surface bundle with which this will be compared.
         * @return \c true if and only if this and the given surface bundle
         * are homeomorphic.
         */
        bool operator == (const NSimpleSurfaceBundle& compare) const;
        /**
         * Determines whether this and the given surface bundle represent
         * different 3-manifolds.
         *
         * @param compare the surface bundle with which this will be compared.
         * @return \c true if and only if this and the given surface bundle
         * are non-homeomorphic.
         */
        bool operator != (const NSimpleSurfaceBundle& compare) const;

        virtual NTriangulation* construct() const;
        NAbelianGroup* homology() const;
        bool isHyperbolic() const;
        std::ostream& writeName(std::ostream& out) const;
        std::ostream& writeTeXName(std::ostream& out) const;
};

/*@}*/

// Inline functions for NSimpleSurfaceBundle

inline NSimpleSurfaceBundle::NSimpleSurfaceBundle(
        int newType) : type_(newType) {
}
inline NSimpleSurfaceBundle::NSimpleSurfaceBundle(
        const NSimpleSurfaceBundle& cloneMe) : NManifold(),
        type_(cloneMe.type_) {
}
inline int NSimpleSurfaceBundle::type() const {
    return type_;
}
inline int NSimpleSurfaceBundle::getType() const {
    return type_;
}
inline bool NSimpleSurfaceBundle::operator ==
        (const NSimpleSurfaceBundle& compare) const {
    return (type_ == compare.type_);
}
inline bool NSimpleSurfaceBundle::operator !=
        (const NSimpleSurfaceBundle& compare) const {
    return (type_ != compare.type_);
}

inline bool NSimpleSurfaceBundle::isHyperbolic() const {
    return false;
}

} // namespace regina

#endif

