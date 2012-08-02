
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

/*! \file subcomplex/npillowtwosphere.h
 *  \brief Deals with 2-spheres made from two faces glued along their
 *  three edges.
 */

#ifndef __NPILLOWTWOSPHERE_H
#ifndef __DOXYGEN
#define __NPILLOWTWOSPHERE_H
#endif

#include "regina-core.h"
#include "shareableobject.h"
#include "maths/nperm4.h"

namespace regina {

class NFace;
class NTriangulation;

/**
 * \weakgroup subcomplex
 * @{
 */

/**
 * Represents a 2-sphere made from two faces glued together along their
 * three edges.  The two faces must be distinct and the three edges of
 * each face must also be distinct.  Neither of the faces may be boundary
 * faces.
 * These two faces together form an embedded 2-sphere in the triangulation
 * (with the exception that two or three points of the sphere corresponding
 * to the face vertices may be identified).
 *
 * This 2-sphere can be cut along and the two resulting 2-sphere
 * boundaries filled in with 3-balls, and the resulting triangulation has
 * the same number of tetrahedra as the original.  If the original
 * 2-sphere was separating, the resulting triangulation will contain the
 * two terms of the corresponding connected sum.
 */
class REGINA_API NPillowTwoSphere : public ShareableObject {
    private:
        NFace* face[2];
            /**< The two faces whose edges are joined. */
        NPerm4 faceMapping;
            /**< A mapping from vertices (0,1,2) of the first face to
                 vertices (0,1,2) of the second face describing how the
                 face boundaries are joined. */
    
    public:
        /**
         * Returns a newly created clone of this structure.
         *
         * @return a newly created clone.
         */
        NPillowTwoSphere* clone() const;

        /**
         * Returns one of the two faces whose boundaries are joined.
         *
         * @param index specifies which of the two faces to return;
         * this must be either 0 or 1.
         * @return the corresponding face.
         */
        NFace* getFace(int index) const;
        /**
         * Returns a permutation describing how the boundaries of the two
         * faces are joined.
         *
         * The permutation will map vertices (0,1,2) of
         * <tt>getFace(0)</tt> to vertices (0,1,2) of
         * <tt>getFace(1)</tt>.  The map will represent how the vertices
         * of the faces are identified by the three edge gluings.
         *
         * @return a permutation describing how the face boundaries are
         * joined.
         */
        NPerm4 getFaceMapping() const;

        /**
         * Determines if the two given faces together form a pillow
         * 2-sphere.
         *
         * \pre The two given faces are distinct.
         *
         * @param face1 the first face to examine.
         * @param face2 the second face to examine.
         * @return a newly created structure containing details of the
         * pillow 2-sphere, or \c null if the given faces do not
         * form a pillow 2-sphere.
         */
        static NPillowTwoSphere* formsPillowTwoSphere(NFace* face1,
            NFace* face2);

        void writeTextShort(std::ostream& out) const;

    private:
        /**
         * Creates a new uninitialised structure.
         */
        NPillowTwoSphere();
};

/*@}*/

// Inline functions for NPillowTwoSphere

inline NPillowTwoSphere::NPillowTwoSphere() {
}
inline NFace* NPillowTwoSphere::getFace(int index) const {
    return face[index];
}
inline NPerm4 NPillowTwoSphere::getFaceMapping() const {
    return faceMapping;
}
inline void NPillowTwoSphere::writeTextShort(std::ostream& out) const {
    out << "Pillow 2-sphere";
}

} // namespace regina

#endif

