
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

/*! \file subcomplex/nsatannulus.h
 *  \brief Deals with saturated two-face annuli within a Seifert fibred space.
 */

#ifndef __NSATANNULUS_H
#ifndef __DOXYGEN
#define __NSATANNULUS_H
#endif

#include "regina-core.h"
#include "maths/nperm4.h"

namespace regina {

class NIsomorphism;
class NMatrix2;
class NTetrahedron;
class NTriangulation;

/**
 * \weakgroup subcomplex
 * @{
 */

/**
 * Represents an annulus formed from a pair of faces in a Seifert fibred
 * space.  This annulus is saturated, i.e., a union of fibres.  More than
 * that, the fibres run parallel to the two boundary edges of the annulus.
 *
 * The annulus is described from one side only.  The description
 * includes an array of indices \a tet[] describing which two tetrahedra
 * provide the faces of the annulus, as well as an array of permutations
 * \a roles[] detailing how the annulus matches up with the individual
 * tetrahedron vertices.
 *
 * The annulus can be drawn as follows, with the upper edge identified
 * with the lower:
 *
 * <pre>
 *            *--->---*
 *            |0  2 / |
 *     First  |    / 1|  Second
 *     face   |   /   |   face
 *            |1 /    |
 *            | / 2  0|
 *            *--->---*
 * </pre>
 *
 * Suppose that \a tet[0] and \a tet[1] are the tetrahedra providing the
 * first and second faces respectively.  Then the markings 0..2 on the
 * first face above correspond to vertices \a roles[0][0..2] of tetrahedron
 * \a tet[0], and likewise the markings 0..2 on the second face above
 * correspond to vertices \a roles[1][0..2] of tetrahedron \a tet[1].
 *
 * Note that the diagram above can also be drawn as follows.
 *
 * <pre>
 *            *--->---*
 *            | \ 2  1|
 *     First  |0 \    |  Second
 *     face   |   \   |   face
 *            |    \ 0|
 *            |1  2 \ |
 *            *--->---*
 * </pre>
 *
 * Note also that the labelling of the tetrahedra and their vertices
 * establishes an orientation on the vertical fibres, as well as a
 * left-to-right direction across the annulus.
 *
 * For convenience we refer to edges \a roles[][0-1] as \e vertical,
 * edges \a roles[][0-2] as \e horizontal, and edge \a roles[][1-2] as
 * \e diagonal.  This is illustrated in the following diagrams.
 *
 * <pre>
 *         V  Horizontal       V   Diagonal
 *         e  *--->---*        e  *--->---*
 *         r  |   g / |        r  |H\ 2  1|
 *         t  |  a / 1|        t  | o\    |
 *         i  | i /   |        i  |  r\   |
 *         c  |D /    |        c  |   i\ 0|
 *         a  | / 2  0|        a  |    z\ |
 *         l  *--->---*        l  *--->---*
 * </pre>
 *
 * \ifacespython The member arrays \a tet and \a roles are accessed for
 * reading through functions \a tet() and \a roles() respectively.  For
 * instance, the first face tetrahedron for the saturated annulus \a a can
 * be accessed as <tt>a.tet(0)</tt>.  These same member arrays are
 * accessed for writing through functions \a setTet() and \a setRoles(),
 * so for instance the second face vertex roles for the saturated annulus
 * \a a can be modified by calling <tt>a.setRoles(1, newRoles)</tt>.
 */
struct REGINA_API NSatAnnulus {
    NTetrahedron* tet[2];
        /**< Describes which tetrahedra provide the first and second
             faces.  See the class notes for details. */
    NPerm4 roles[2];
        /**< Describes how the first and second faces match up with
             individual tetrahedron vertices.  See the class notes for
             details. */

    /**
     * Creates a new uninitialised structure.  Both tetrahedra will be
     * set to null pointers.
     */
    NSatAnnulus();
    /**
     * Creates a clone of the given structure.
     *
     * @param cloneMe the structure to clone.
     */
    NSatAnnulus(const NSatAnnulus& cloneMe);
    /**
     * Creates a new structure initialised to the given values.  See the
     * class notes for what the various tetrahedra and permutations mean.
     *
     * @param t0 the tetrahedron to assign to \a tet[0].
     * @param r0 the permutation to assign to \a roles[0].
     * @param t1 the tetrahedron to assign to \a tet[1].
     * @param r1 the permutation to assign to \a roles[1].
     */
    NSatAnnulus(NTetrahedron* t0, NPerm4 r0, NTetrahedron* t1, NPerm4 r1);
    /**
     * Makes this equal to a clone of the given structure.
     *
     * @param cloneMe the structure to clone.
     * @return a reference to this structure.
     */
    NSatAnnulus& operator = (const NSatAnnulus& cloneMe);

    /**
     * Determines whether or not this and the given structure describe
     * the same annulus with the same representation.  This requires
     * both structures to have identical \a tet[] and \a roles[] arrays.
     *
     * @param other the structure to compare with this.
     * @return \c true if the structures describe the same annulus with
     * the same representation, or \c false if they do not.
     */
    bool operator == (const NSatAnnulus& other) const;
    /**
     * Determines whether or not this and the given structure describe
     * the same annulus with the same representation.  This requires
     * both structures to have identical \a tet[] and \a roles[] arrays.
     *
     * @param other the structure to compare with this.
     * @return \c true if the structures do not describe the same annulus
     * with the same representation, or \c false if they do.
     */
    bool operator != (const NSatAnnulus& other) const;

    /**
     * Determines how many faces of this annulus lie on the boundary
     * of the triangulation.
     *
     * Note that this routine can also be used as a boolean function
     * to determine whether any faces of the annulus lie on the
     * triangulation boundary.
     *
     * @return the number of faces of this annulus that lie on the boundary
     * of the triangulation; this will be 0, 1 or 2.
     */
    unsigned meetsBoundary() const;

    /**
     * Converts this into a representation of the same annulus from the
     * other side.  The first and second faces and their 0..2 markings
     * (as described in the class notes) remain unchanged.  However, the
     * two tetrahedra that are used to describe the annulus will be
     * replaced by their counterparts on the other side of the annulus
     * (i.e., the two new tetrahedra that meet the two original tetrahedra
     * along the annulus itself).
     *
     * \pre Neither face of this annulus is a boundary face of the
     * triangulation.
     */
    void switchSides();
    /**
     * Returns a representation of the same annulus from the other side.
     * This structure will not be changed.  See switchSides() for further
     * details.
     *
     * \pre Neither face of this annulus is a boundary face of the
     * triangulation.
     *
     * @return a new representation of this annulus from the other side.
     */
    NSatAnnulus otherSide() const;

    /**
     * Reverses the direction of the vertical fibres in this annulus
     * representation.  The first and second faces (as described in the
     * class notes) will remain unchanged, but the markings 0 and 1 on
     * each face will be switched.
     */
    void reflectVertical();
    /**
     * Returns a representation of this annulus in which the vertical
     * direction of the fibres has been reversed.  This structure will
     * not be changed.  See reflectVertical() for further details.
     *
     * @return a new representation of this annulus in which fibres have
     * been reversed.
     */
    NSatAnnulus verticalReflection() const;

    /**
     * Performs a left-to-right reflection of this annulus
     * representation.  The vertical direction of the fibres will remain
     * unchanged, but the first and second faces will be switched (and
     * the 0..2 markings changed to compensate).
     */
    void reflectHorizontal();
    /**
     * Returns a left-to-right reflected representation of this annulus.
     * This structure will not be changed.  See reflectHorizontal() for
     * further details.
     *
     * @return a new left-to-right reflection of this annulus.
     */
    NSatAnnulus horizontalReflection() const;

    /**
     * Rotates the representation of this annulus by 180 degrees.
     * This has the effect of switching the first and second faces and
     * also reversing the direction of the vertical fibres.
     *
     * Calling this routine is equivalent to calling reflectVertical() and
     * then reflectHorizontal().
     */
    void rotateHalfTurn();
    /**
     * Returns a 180 degree rotated representation of this annulus.
     * This structure will not be changed.  See rotateHalfTurn() for
     * further details.
     *
     * @return a new 180 degree rotation of this annulus.
     */
    NSatAnnulus halfTurnRotation() const;

    /**
     * Determines whether this and the given annulus are adjacent,
     * possibly modulo vertical or horizontal reflections.  That is,
     * this routine determines whether this and the given structure
     * represent opposite sides of the same saturated annulus, where the
     * fibres for both structures are consistent (though possibly reversed).
     * See switchSides() for details on what "opposite sides" means in
     * this context, and see reflectVertical() and reflectHorizontal()
     * for descriptions of the various types of reflection.
     *
     * Information regarding reflections is returned via the two boolean
     * pointers \a refVert and \a refHoriz.  If the two annuli are
     * identically opposite each other as described by switchSides(),
     * both booleans will be set to \c false.  If the two annuli are
     * identically opposite after one undergoes a vertical and/or
     * horizontal reflection, then the booleans \a refVert and/or
     * \a refHoriz will be set to \c true accordingly.
     *
     * The critical difference between this routine and isJoined() is
     * that this routine insists that the fibres on each annulus be
     * consistent.  This routine is thus suitable for examining joins
     * between different sections of the same Seifert fibred space,
     * for example.
     *
     * \ifacespython This routine only takes a single argument (the
     * annulus \a other).  The return value is a tuple of three
     * booleans: the usual return value, the value returned in \a refVert,
     * and the value returned in \a refHoriz.
     *
     * @param other the annulus to compare with this.
     * @param refVert returns information on whether the annuli are
     * adjacent modulo a vertical reflection.  This is set to \c true
     * if a vertical reflection is required and \c false if it is not.
     * If no adjacency was found at all, this boolean is not touched.
     * A null pointer may be passed, in which case this information will
     * not be returned at all.
     * @param refHoriz returns information on whether the annuli are
     * adjacent modulo a horizontal reflection.  This is set to \c true
     * if a horizontal reflection is required and \c false if it is not.
     * If no adjacency was found at all, this boolean is not touched.
     * A null pointer may be passed, in which case this information will
     * not be returned at all.
     * @return \c true if some adjacency was found (either with or
     * without reflections), or \c false if no adjacency was found at all.
     */
    bool isAdjacent(const NSatAnnulus& other, bool* refVert, bool* refHoriz)
        const;

    /**
     * Determines whether this and the given annulus are joined in some
     * form, even if the fibres on each annulus are not consistent.
     *
     * This routine treats each annulus as though its boundaries are
     * identified to form a torus (though it does not actually test
     * whether this is true).  It then examines whether this and the
     * given annulus represent opposite sides of the same torus.
     * More specifically, it tests whether both annuli are formed from
     * the same pair of faces, and whether the mapping of 0/1/2 markings
     * from one annulus to the other is the same for each face.  Note that
     * the faces are allowed to be switched (i.e., the first face of one
     * annulus may be the second face of the other).
     *
     * The critical difference between this routine and isAdjacent() is
     * that this routine allows the fibres on each annulus to be
     * inconsistent.  This routine is thus suitable for examining joins
     * between different Seifert fibred blocks in a graph manifold, for
     * example.
     *
     * If the two annuli are joined, the precise relationship between
     * the curves on each annulus will be returned in the matrix
     * \a matching.  Specifically, let \a x and \a y be the oriented
     * curves running from markings 0-1 and 0-2 respectively on the
     * first face of this annulus.  Likewise, let \a x' and \a y' run
     * from markings 0-1 and 0-2 respectively on the first face of the
     * annulus \a other.  Then the joining between the two annuli can
     * be expressed as follows:
     *
     * <pre>
     *     [x ]                [x']
     *     [  ]  =  matching * [  ].
     *     [y ]                [y']
     * </pre>
     *
     * @param other the annulus to compare with this.
     * @param matching returns details on how the curves on each annulus
     * are related.  If the this and the given annulus are not joined,
     * then this matrix is not touched.
     * @return \c true if this and the given annulus are found to be
     * joined, or \c false if they are not.
     */
    bool isJoined(const NSatAnnulus& other, NMatrix2& matching) const;

    /**
     * Determines whether this annulus has its boundaries identified to
     * form an embedded two-sided torus within the surrounding triangulation.
     *
     * It will be verified that:
     * - the two faces of this annulus are joined along all three pairs
     *   of edges to form a torus;
     * - the three edges of this torus remain distinct (i.e., different edges
     *   of the torus do not become identified within the larger triangulation);
     * - this torus is two-sided within the surrounding triangulation.
     *
     * @return \c true if this annulus forms an embedded two-sided torus as
     * described above, or \c false if it does not.
     */
    bool isTwoSidedTorus() const;

    /**
     * Adjusts this annulus representation according to the given
     * isomorphism between triangulations.
     *
     * The given isomorphism must describe a mapping from \a originalTri
     * to \a newTri, and this annulus must refer to tetrahedra in
     * \a originalTri.  This routine will adjust this annulus according
     * to the given isomorphism, so that it refers to the corresponding
     * tetrahedra in \a newTri (with the \a roles permutations also
     * updated accordingly).
     *
     * \pre This annulus refers to tetrahedra in \a originalTri, and
     * \a iso describes a mapping from \a originalTri to \a newTri.
     *
     * @param originalTri the triangulation currently used by this
     * annulus representation.
     * @param iso the mapping from \a originalTri to \a newTri.
     * @param newTri the triangulation to be used by the updated annulus
     * representation.
     */
    void transform(const NTriangulation* originalTri,
            const NIsomorphism* iso, NTriangulation* newTri);
    /**
     * Returns the image of this annulus representation under the given
     * isomorphism between triangulations.  This annulus representation
     * will not be changed.  See transform() for further details.
     *
     * @param originalTri the triangulation currently used by this
     * annulus representation.
     * @param iso the mapping from \a originalTri to \a newTri.
     * @param newTri the triangulation to be used by the new annulus
     * representation.
     */
    NSatAnnulus image(const NTriangulation* originalTri,
            const NIsomorphism* iso, NTriangulation* newTri) const;

    /**
     * Attaches a layered solid torus to the this saturated annulus.
     *
     * The layered solid torus will be attached so that the
     * given values \a alpha and \a beta describe how the
     * meridinal disc cuts the vertical and horizontal edges of the
     * annulus respectively.
     *
     * The result will effectively insert an (\a alpha, \a beta)
     * exceptional fibre into the Seifert fibred space space, where
     * the vertical edges run parallel to the fibres and the horizontal
     * edges represent the base orbifold.  The sign of the fibre is
     * consistent with the fibre inserted by NSatLST::adjustSFS()
     * (in particular, negating \a beta will negate the fibre).
     *
     * In the case of a (2,1) fibre, the layered solid torus will be
     * degenerate (i.e., the two faces of the annulus will simply be
     * joined together).
     *
     * \pre The given value \a alpha is not zero.
     * \pre The given values \a alpha and \a beta are coprime.
     *
     * @param tri the triangulation into which the new tetrahedra should
     * be inserted.
     * @param alpha describes how the meridinal disc of the torus should
     * cut the vertical edges.  This may be positive or negative.
     * @param beta describes how the meridinal disc of the torus should
     * cut the horizontal edges.  Again this may be positive or negative.
     */
    void attachLST(NTriangulation* tri, long alpha, long beta) const;
};

/*@}*/

// Inline functions for NSatAnnulus

inline NSatAnnulus::NSatAnnulus() {
    tet[0] = tet[1] = 0;
}

inline NSatAnnulus::NSatAnnulus(const NSatAnnulus& cloneMe) {
    tet[0] = cloneMe.tet[0]; tet[1] = cloneMe.tet[1];
    roles[0] = cloneMe.roles[0]; roles[1] = cloneMe.roles[1];
}

inline NSatAnnulus::NSatAnnulus(NTetrahedron* t0, NPerm4 r0,
        NTetrahedron* t1, NPerm4 r1) {
    tet[0] = t0; tet[1] = t1;
    roles[0] = r0; roles[1] = r1;
}

inline NSatAnnulus& NSatAnnulus::operator = (const NSatAnnulus& cloneMe) {
    tet[0] = cloneMe.tet[0]; tet[1] = cloneMe.tet[1];
    roles[0] = cloneMe.roles[0]; roles[1] = cloneMe.roles[1];
    return *this;
}

inline bool NSatAnnulus::operator == (const NSatAnnulus& other) const {
    return (tet[0] == other.tet[0] && tet[1] == other.tet[1] &&
            roles[0] == other.roles[0] && roles[1] == other.roles[1]);
}

inline bool NSatAnnulus::operator != (const NSatAnnulus& other) const {
    return (tet[0] != other.tet[0] || tet[1] != other.tet[1] ||
            roles[0] != other.roles[0] || roles[1] != other.roles[1]);
}

inline NSatAnnulus NSatAnnulus::otherSide() const {
    NSatAnnulus a(*this);
    a.switchSides();
    return a;
}

inline void NSatAnnulus::reflectVertical() {
    roles[0] = roles[0] * NPerm4(0, 1);
    roles[1] = roles[1] * NPerm4(0, 1);
}

inline NSatAnnulus NSatAnnulus::verticalReflection() const {
    return NSatAnnulus(tet[0], roles[0] * NPerm4(0, 1),
                       tet[1], roles[1] * NPerm4(0, 1));
}

inline void NSatAnnulus::reflectHorizontal() {
    NTetrahedron* t = tet[0];
    tet[0] = tet[1];
    tet[1] = t;

    NPerm4 r = roles[0];
    roles[0] = roles[1] * NPerm4(0, 1);
    roles[1] = r * NPerm4(0, 1);
}

inline NSatAnnulus NSatAnnulus::horizontalReflection() const {
    return NSatAnnulus(tet[1], roles[1] * NPerm4(0, 1),
                       tet[0], roles[0] * NPerm4(0, 1));
}

inline void NSatAnnulus::rotateHalfTurn() {
    NTetrahedron* t = tet[0];
    tet[0] = tet[1];
    tet[1] = t;

    NPerm4 r = roles[0];
    roles[0] = roles[1];
    roles[1] = r;
}

inline NSatAnnulus NSatAnnulus::halfTurnRotation() const {
    return NSatAnnulus(tet[1], roles[1], tet[0], roles[0]);
}

inline NSatAnnulus NSatAnnulus::image(const NTriangulation* originalTri,
        const NIsomorphism* iso, NTriangulation* newTri) const {
    NSatAnnulus a(*this);
    a.transform(originalTri, iso, newTri);
    return a;
}

} // namespace regina

#endif

