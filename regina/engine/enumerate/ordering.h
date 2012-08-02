
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

/*! \file enumerate/ordering.h
 *  \brief Provides different ways of sorting hyperplanes (or matching
 *  equations) when performing normal surface enumeration.
 */

#ifndef __ORDERING_H
#ifndef __DOXYGEN
#define __ORDERING_H
#endif

#include "regina-core.h"
#include "maths/nmatrixint.h"

namespace regina {

/**
 * \weakgroup enumerate
 * @{
 */

/**
 * A comparison object that sorts hyperplanes by position vectors.
 * This ordering is described in "Optimizing the double description
 * method for normal surface enumeration", B.A. Burton,
 * Mathematics of Computation 79 (2010), 453-484.
 *
 * This comparison object is used to sort hyperplanes into a good
 * order before enumerating vertex or fundamental normal surfaces
 * A hyperplane is described by a row of the \a subspace matrix,
 * as passed to an enumeration routine such as
 * NDoubleDescription::enumerateExtremalRays() or
 * NHilbertDual::enumerateHilbertBasis().
 *
 * The ordering is defined as follows.  For each hyperplane, we
 * create a <i>position vector</i> (h_1, ..., h_f), where h_i is 0 if the
 * hyperplane contains the ith coordinate axis, or 1 if not.
 * We then compare these position vectors lexicographically.
 *
 * \ifacespython Not present.
 */
class NPosOrder {
    private:
        const NMatrixInt& matrix_;
            /**< The \a subspace matrix as passed to the enumeration routine. */

    public:
        /**
         * Creates a new helper object for comparing hyperplanes.
         *
         * @param matrix the \a subspace matrix as passed to
         * the normal surface enumeration routine.
         */
        inline NPosOrder(const NMatrixInt& matrix);

        /**
         * Determines whether the hyperplane described by
         * row \a i of the matrix is smaller
         * than the hyperplane described by row \a j.
         * Here "smaller" is defined by position vectors;
         * see the NPosOrder class notes for details.
         *
         * @param i the first matrix row index; this must be between
         * 0 and matrix.rows()-1 inclusive, where \a matrix is
         * the matrix passed to the class constructor.
         * @param j the second matrix row index; this must also be
         * between 0 and matrix.rows()-1 inclusive.
         * @return \c true if and only if the hyperplane described by
         * row \a i is smaller than the hyperplane described by row \a j.
         */
        inline bool operator () (int i, int j) const;
};

/*@}*/

// Inline functions for NPosOrder

inline NPosOrder::NPosOrder(const NMatrixInt& matrix) :
        matrix_(matrix) {
}

inline bool NPosOrder::operator () (
        int i, int j) const {
    for (unsigned c = 0; c < matrix_.columns(); ++c) {
        if (matrix_.entry(i, c) == 0 && matrix_.entry(j, c) != 0)
            return true;
        if (matrix_.entry(i, c) != 0 && matrix_.entry(j, c) == 0)
            return false;
    }
    return false;
}

} // namespace regina

#endif

