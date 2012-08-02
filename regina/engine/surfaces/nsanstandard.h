
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

/*! \file surfaces/nsanstandard.h
 *  \brief Implements almost normal surface vectors using standard
 *  triangle-quad-oct coordinates.
 */

#ifndef __NSANSTANDARD_H
#ifndef __DOXYGEN
#define __NSANSTANDARD_H
#endif

#include "regina-core.h"
#include "surfaces/nnormalsurface.h"

namespace regina {

class NMatrixInt;

/**
 * \weakgroup surfaces
 * @{
 */

/**
 * An almost normal surface vector using standard triangle-quad-oct
 * coordinates.
 *
 * If there are \a t tetrahedra in the underlying
 * triangulation, there must be precisely 10<i>t</i> coordinates.
 * The first ten coordinates will be for the first tetrahedron, the
 * next ten for the second tetrahedron and so on.  For each
 * tetrahedron, the first four represent the number of
 * triangular discs about vertex 0, 1, 2 and 3, the next
 * three represent the number of quadrilateral discs of type 0,
 * 1 and 2 (see NNormalSurface::getQuadCoord()) and the final three
 * represent the number of octagonal discs of type 0, 1 and 2 (see
 * NNormalSurface::getOctCoord()).
 *
 * \ifacespython Not present.
 */
class REGINA_API NNormalSurfaceVectorANStandard : public NNormalSurfaceVector {
    public:
        /**
         * Creates a new vector all of whose entries are initialised to
         * zero.
         *
         * @param length the number of elements in the new vector.
         */
        NNormalSurfaceVectorANStandard(unsigned length);
        /**
         * Creates a new vector that is a clone of the given vector.
         *
         * @param cloneMe the vector to clone.
         */
        NNormalSurfaceVectorANStandard(const NVector<NLargeInteger>& cloneMe);

        virtual bool allowsAlmostNormal() const;
        virtual bool allowsSpun() const;
        virtual bool allowsOriented() const;

        virtual NLargeInteger getTriangleCoord(unsigned long tetIndex,
            int vertex, NTriangulation* triang) const;
        virtual NLargeInteger getQuadCoord(unsigned long tetIndex,
            int quadType, NTriangulation* triang) const;
        virtual NLargeInteger getOctCoord(unsigned long tetIndex,
            int octType, NTriangulation* triang) const;
        virtual NLargeInteger getEdgeWeight(unsigned long edgeIndex,
            NTriangulation* triang) const;
        virtual NLargeInteger getFaceArcs(unsigned long faceIndex,
            int faceVertex, NTriangulation* triang) const;

        virtual NNormalSurfaceVector* clone() const;

        static NNormalSurfaceVector* makeZeroVector(
            const NTriangulation* triangulation);
        static NMatrixInt* makeMatchingEquations(NTriangulation* triangulation);
        static NEnumConstraintList* makeEmbeddedConstraints(
            NTriangulation* triangulation);
};

/*@}*/

// Inline functions for NNormalSurfaceVectorANStandard

inline NNormalSurfaceVectorANStandard::NNormalSurfaceVectorANStandard(
        unsigned length) : NNormalSurfaceVector(length) {
}
inline NNormalSurfaceVectorANStandard::NNormalSurfaceVectorANStandard(
        const NVector<NLargeInteger>& cloneMe) :
        NNormalSurfaceVector(cloneMe) {
}

inline NLargeInteger NNormalSurfaceVectorANStandard::getTriangleCoord(
        unsigned long tetIndex, int vertex, NTriangulation*) const {
    return (*this)[10 * tetIndex + vertex];
}
inline NLargeInteger NNormalSurfaceVectorANStandard::getQuadCoord(
        unsigned long tetIndex, int quadType, NTriangulation*) const {
    return (*this)[10 * tetIndex + 4 + quadType];
}
inline NLargeInteger NNormalSurfaceVectorANStandard::getOctCoord(
        unsigned long tetIndex, int octType, NTriangulation*) const {
    return (*this)[10 * tetIndex + 7 + octType];
}

} // namespace regina

#endif

