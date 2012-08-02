
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

/*! \file subcomplex/nlayeredsurfacebundle.h
 *  \brief Deals with layered surface bundle triangulations.
 */

#ifndef __NLAYEREDSURFACEBUNDLE_H
#ifndef __DOXYGEN
#define __NLAYEREDSURFACEBUNDLE_H
#endif

#include <memory>
#include "regina-core.h"
#include "maths/nmatrix2.h"
#include "subcomplex/nstandardtri.h"
#include "triangulation/ntriangulation.h"

namespace regina {

class NTxICore;

/**
 * \weakgroup subcomplex
 * @{
 */

/**
 * Describes a layered torus bundle.  This is a triangulation of a
 * torus bundle over the circle formed as follows.
 *
 * We begin with a thin I-bundle over the torus, i.e,. a triangulation
 * of the product <tt>T x I</tt> that is only one tetrahedron thick.
 * This is referred to as the \a core, and is described by an object
 * of type NTxICore.
 *
 * We then identify the upper and lower torus boundaries of this core
 * according to some homeomorphism of the torus.  This may be impossible
 * due to incompatible boundary edges, and so we allow a layering of
 * tetrahedra over one of the boundari tori in order to adjust the
 * boundary edges accordingly.  Layerings are described in more detail
 * in the NLayering class.
 *
 * Given the parameters of the core <tt>T x I</tt> and the specific
 * layering, the monodromy for this torus bundle over the circle can be
 * calculated.  The getManifold() routine returns details of the
 * corresponding 3-manifold.
 *
 * All optional NStandardTriangulation routines are implemented for this
 * class.
 *
 * \testpart
 */
class REGINA_API NLayeredTorusBundle : public NStandardTriangulation {
    private:
        const NTxICore& core_;
            /**< The core <tt>T x I</tt> triangulation whose boundaries
                 are joined (possibly via a layering of tetrahedra). */
        NIsomorphism* coreIso_;
            /**< Describes how the tetrahedra and vertices of the core
                 <tt>T x I</tt> triangulation returned by NTxICore::core()
                 map to the tetrahedra and vertices of the larger layered
                 torus bundle under consideration. */
        NMatrix2 reln_;
            /**< Describes how the layering of tetrahedra maps the
                 lower boundary curves to the upper boundary curves.
                 See layeringReln() for details. */

    public:
        /**
         * Destroys this layered torus bundle and all of its internal
         * components.
         */
        virtual ~NLayeredTorusBundle();

        /**
         * Returns the <tt>T x I</tt> triangulation at the core of this
         * layered surface bundle.  This is the product <tt>T x I</tt>
         * whose boundaries are joined (possibly via some layering of
         * tetrahedra).
         *
         * Note that the triangulation returned by NTxICore::core()
         * (that is, NLayeredSurfaceBundle::core().core()) may
         * well use different tetrahedron and vertex numbers.  That is,
         * an isomorphic copy of it appears within this layered surface
         * bundle but the individual tetrahedra and vertices may have
         * been permuted.  For a precise mapping from the NTxICore::core()
         * triangulation to this triangulation, see the routine coreIso().
         *
         * @return the core <tt>T x I</tt> triangulation.
         */
        const NTxICore& core() const;

        /**
         * Returns the isomorphism describing how the core <tt>T x I</tt>
         * appears as a subcomplex of this layered surface bundle.
         *
         * As described in the core() notes, the core <tt>T x I</tt>
         * triangulation returned by NTxICore::core() appears within this
         * layered surface bundle, but not necessarily with the same
         * tetrahedron or vertex numbers.
         *
         * This routine returns an isomorphism that maps the tetrahedra
         * and vertices of the core <tt>T x I</tt> triangulation (as
         * returned by NLayeredSurfaceBundle::core().core()) to the
         * tetrahedra and vertices of this overall layered surface bundle.
         *
         * The isomorphism that is returned belongs to this object, and
         * should not be modified or destroyed.
         *
         * @return the isomorphism from the core <tt>T x I</tt> to this
         * layered surface bundle.
         */
        const NIsomorphism* coreIso() const;

        /**
         * Returns a 2-by-2 matrix describing how the layering of
         * tetrahedra relates curves on the two torus boundaries of the
         * core <tt>T x I</tt>.
         *
         * The NTxICore class documentation describes generating \a alpha
         * and \a beta curves on the two torus boundaries of the core
         * <tt>T x I</tt> (which are referred to as the \e upper and
         * \e lower boundaries).  The two boundary tori are parallel in
         * two directions: through the core, and through the layering.
         * It is desirable to know the parallel relationship between
         * the two sets of boundary curves in each direction.
         *
         * The relationship through the core is already described by
         * NTxICore::parallelReln().  This routine describes the
         * relationship through the layering.
         *
         * Let \a a_u and \a b_u be the \a alpha and \a beta curves on
         * the upper boundary torus, and let \a a_l and \a b_l be the
         * \a alpha and \a beta curves on the lower boundary torus.
         * Suppose that the upper \a alpha is parallel to
         * \a w.\a a_l + \a x.\a b_l, and that the upper \a beta is
         * parallel to \a y.\a a_l + \a z.\a b_l.  Then the matrix
         * returned will be
         *
         * <pre>
         *     [ w  x ]
         *     [      ] .
         *     [ y  z ]
         * </pre>
         *
         * In other words,
         *
         * <pre>
         *     [ a_u ]                      [ a_l ]
         *     [     ]  =  layeringReln() * [     ] .
         *     [ b_u ]                      [ b_l ]
         * </pre>
         *
         * It can be observed that this matrix expresses the upper
         * boundary curves in terms of the lower, whereas
         * NTxICore::parallelReln() expresses the lower boundary curves
         * in terms of the upper.  This means that the monodromy
         * describing the overall torus bundle over the circle can be
         * calculated as
         * <pre>
         *     M  =  layeringReln() * core().parallelReln()
         * </pre>
         * or alternatively using the similar matrix
         * <pre>
         *     M'  =  core().parallelReln() * layeringReln() .
         * </pre>
         *
         * Note that in the degenerate case where there is no layering at
         * all, this matrix is still perfectly well defined; in this
         * case it describes a direct identification between the upper
         * and lower boundary tori.
         *
         * @return the relationship through the layering between the
         * upper and lower boundary curves of the core <tt>T x I</tt>.
         */
        const NMatrix2& layeringReln() const;

        /**
         * Determines if the given triangulation is a layered surface bundle.
         *
         * @param tri the triangulation to examine.
         * @return a newly created structure containing details of the
         * layered surface bundle, or \c null if the given triangulation
         * is not a layered surface bundle.
         */
        static NLayeredTorusBundle* isLayeredTorusBundle(NTriangulation* tri);

        NManifold* getManifold() const;
        NAbelianGroup* getHomologyH1() const;
        std::ostream& writeName(std::ostream& out) const;
        std::ostream& writeTeXName(std::ostream& out) const;
        void writeTextLong(std::ostream& out) const;

    private:
        /**
         * Creates a new structure based upon the given core <tt>T x I</tt>
         * triangulation.  Aside from this core, the structure will
         * remain uninitialised.
         *
         * Note that only a reference to the core <tt>T x I</tt> is stored.
         * This class does not manage the life span of the core; it is
         * assumed that the core will remain in existence for at least
         * as long as this object does.  (Usually the core is a static or
         * global variable that is not destroyed until the program exits.)
         *
         * @param whichCore a reference to the core <tt>T x I</tt>
         * triangulation upon which this layered surface bundle is based.
         */
        NLayeredTorusBundle(const NTxICore& whichCore);

        /**
         * Contains code common to both writeName() and writeTeXName().
         *
         * @param out the output stream to which to write.
         * @param tex \c true if this routine is called from
         * writeTeXName() or \c false if it is called from writeName().
         * @return a reference to \a out.
         */
        std::ostream& writeCommonName(std::ostream& out, bool tex) const;

        /**
         * Internal to isLayeredTorusBundle().  Determines if the given
         * triangulation is a layered surface bundle with the given core
         * <tt>T x I</tt> triangulation (up to isomorphism).
         *
         * @param tri the triangulation to examine.
         * @param core the core <tt>T x I</tt> to search for.
         * @return a newly created structure containing details of the
         * layered surface bundle, or \c null if the given triangulation is
         * not a layered surface bundle with the given <tt>T x I</tt> core.
         */
        static NLayeredTorusBundle* hunt(NTriangulation* tri,
            const NTxICore& core);
};

/*@}*/

// Inline functions for NLayeredTorusBundle

inline NLayeredTorusBundle::NLayeredTorusBundle(const NTxICore& whichCore) :
        core_(whichCore), coreIso_(0) {
}

inline const NTxICore& NLayeredTorusBundle::core() const {
    return core_;
}

inline const NIsomorphism* NLayeredTorusBundle::coreIso() const {
    return coreIso_;
}

inline const NMatrix2& NLayeredTorusBundle::layeringReln() const {
    return reln_;
}

inline std::ostream& NLayeredTorusBundle::writeName(std::ostream& out) const {
    return writeCommonName(out, false);
}

inline std::ostream& NLayeredTorusBundle::writeTeXName(std::ostream& out)
        const {
    return writeCommonName(out, true);
}

} // namespace regina

#endif

