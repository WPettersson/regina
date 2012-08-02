
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

/*! \file subcomplex/nl31pillow.h
 *  \brief Deals with triangular pillow L(3,1) components of a triangulation.
 */

#ifndef __NL31PILLOW_H
#ifndef __DOXYGEN
#define __NL31PILLOW_H
#endif

#include "regina-core.h"
#include "subcomplex/nstandardtri.h"

namespace regina {

class NTetrahedron;

/**
 * \weakgroup subcomplex
 * @{
 */

/**
 * Represents a triangular pillow L(3,1) component of a triangulation.
 *
 * A triangular pillow L(3,1) is a two-tetrahedron two-vertex
 * triangulation of the lens space L(3,1) formed as follows.
 *
 * A triangular pillow is formed from two tetrahedra with a single
 * degree three vertex in the interior of the pillow.  The two boundary
 * faces of this pillow are then identified with a one-third twist.
 *
 * All optional NStandardTriangulation routines are implemented for this
 * class.
 *
 * \testpart
 */
class REGINA_API NL31Pillow : public NStandardTriangulation {
    private:
        NTetrahedron* tet[2];
            /**< The two tetrahedra in the triangular pillow. */
        unsigned interior[2];
            /**< The vertex of each tetrahedron that corresponds to the
                 interior vertex of the triangular pillow. */

    public:
        /**
         * Destroys this structure.
         */
        virtual ~NL31Pillow();
        /**
         * Returns a newly created clone of this structure.
         *
         * @return a newly created clone.
         */
        NL31Pillow* clone() const;

        /**
         * Returns one of the two tetrahedra involved in this structure.
         *
         * @param whichTet specifies which tetrahedron to return; this
         * must be either 0 or 1.
         * @return the requested tetrahedron.
         */
        NTetrahedron* getTetrahedron(int whichTet) const;
        /**
         * Returns the vertex number of the given tetrahedron
         * corresponding to the degree three vertex in the interior of
         * the triangular pillow.  See the general class notes for
         * further details.
         *
         * The specific tetrahedron to examine is determined by the
         * argument \a whichTet; this will be the tetrahedron
         * <tt>getTetrahedron(whichTet)</tt>.
         *
         * @param whichTet specifies which tetrahedron to examine;
         * this must be either 0 or 1.
         * @return the vertex of tetrahedron \a whichTet corresponding
         * to the vertex in the interior of the triangular pillow; this
         * will be between 0 and 3 inclusive.
         */
        unsigned getInteriorVertex(int whichTet) const;

        /**
         * Determines if the given triangulation component is a
         * triangular pillow L(3,1).
         *
         * @param comp the triangulation component to examine.
         * @return a newly created structure containing details of the
         * triangular pillow L(3,1), or \c null if the given component is
         * not a triangular pillow L(3,1).
         */
        static NL31Pillow* isL31Pillow(const NComponent* comp);

        NManifold* getManifold() const;
        NAbelianGroup* getHomologyH1() const;
        std::ostream& writeName(std::ostream& out) const;
        std::ostream& writeTeXName(std::ostream& out) const;
        void writeTextLong(std::ostream& out) const;

    private:
        /**
         * Creates a new uninitialised structure.
         */
        NL31Pillow();
};

/*@}*/

// Inline functions for NL31Pillow

inline NL31Pillow::NL31Pillow() {
}
inline NL31Pillow::~NL31Pillow() {
}

inline NTetrahedron* NL31Pillow::getTetrahedron(int whichTet) const {
    return tet[whichTet];
}
inline unsigned NL31Pillow::getInteriorVertex(int whichTet) const {
    return interior[whichTet];
}
inline std::ostream& NL31Pillow::writeName(std::ostream& out) const {
    return out << "L'(3,1)";
}
inline std::ostream& NL31Pillow::writeTeXName(std::ostream& out) const {
    return out << "L'_{3,1}";
}
inline void NL31Pillow::writeTextLong(std::ostream& out) const {
    out << "Triangular pillow lens space L(3,1)";
}

} // namespace regina

#endif

