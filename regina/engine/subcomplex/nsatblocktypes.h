
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

/*! \file subcomplex/nsatblocktypes.h
 *  \brief Describes several types of saturated blocks within Seifert fibred
 *  space triangulations.
 */

#ifndef __NSATBLOCKTYPES_H
#ifndef __DOXYGEN
#define __NSATBLOCKTYPES_H
#endif

#include "regina-core.h"
#include "subcomplex/nsatblock.h"

namespace regina {

class NLayeredSolidTorus;

/**
 * \weakgroup subcomplex
 * @{
 */

/**
 * A degenerate zero-tetrahedron saturated block that corresponds to
 * attaching a Mobius band to a single annulus boundary.
 *
 * This is a degenerate case of the layered solid torus (see the class
 * NSatLST), where instead of joining a solid torus to an annulus
 * boundary we join a Mobius band.  The Mobius band can be thought of as
 * a zero-tetrahedron solid torus with two boundary faces, which in fact
 * are opposite sides of the same face.  By attaching a zero-tetrahedron
 * Mobius band to an annulus boundary, we are effectively joining the
 * two faces of the annulus together.
 *
 * The meridinal disc of this zero-tetrahedron solid torus meets the
 * three edges of the annulus in 1, 1 and 2 places, so it is in fact
 * a degenerate (1,1,2) layered solid torus.  Note that the weight 2 edge
 * is the boundary edge of the Mobius strip.
 */
class REGINA_API NSatMobius : public NSatBlock {
    private:
        int position_;
            /**< Describes how the Mobius band is attached to the
                 boundary annulus.  This can take the value 0, 1 or 2.
                 See the position() documentation for further details. */

    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatMobius(const NSatMobius& cloneMe);

        /**
         * Describes how the Mobius band is attached to the
         * boundary annulus.
         *
         * The class notes discuss the weight two edge of the Mobius band
         * (or equivalently the boundary edge of the Mobius band).  The
         * return value of this routine indicates which edge of the
         * boundary annulus this weight two edge is joined to.
         *
         * In the NSatAnnulus class notes, the three edges of the
         * annulus are denoted vertical, horizontal and boundary, and
         * the vertices of each face are given markings 0, 1 and 2.
         *
         * The return value of this routine takes the value 0, 1 or 2 as
         * follows:
         * - 0 means that the weight two edge is joined to the diagonal
         *   edge of the annulus (markings 1 and 2);
         * - 1 means that the weight two edge is joined to the horizontal
         *   edge of the annulus (markings 0 and 2);
         * - 2 means that the weight two edge is joined to the vertical
         *   edge of the annulus (markings 0 and 1).
         *
         * @return the value 0, 1 or 2 as described above.
         */
        int position() const;

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (Mobius band).  This routine is
         * a specific case of NSatBlock::isBlock(); see that routine for
         * further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatMobius* isBlockMobius(const NSatAnnulus& annulus,
            TetList& avoidTets);

    private:
        /**
         * Constructs a partially initialised block.  The boundary
         * annulus will remain uninitialised, and must be initialised
         * before this block can be used.
         *
         * @param position indicates which edge of the boundary annulus
         * meets the weight two edge of the Mobius strip, as described in
         * the position() member function documentation.  This value
         * must be 0, 1 or 2.
         */
        NSatMobius(int position);
};

/**
 * A saturated block that is a layered solid torus.  See the
 * NLayeredSolidTorus class for details.
 *
 * The three boundary edges of the layered solid torus are attached to
 * the vertical, horizontal and diagonal edges of the boundary annulus;
 * see the NSatAnnulus class notes for details on precisely what
 * vertical, horizontal and diagonal mean.
 */
class REGINA_API NSatLST : public NSatBlock {
    private:
        NLayeredSolidTorus* lst_;
            /**< Contains details of the layered solid torus that this
                 block represents. */
        NPerm4 roles_;
            /**< Describes how the layered solid torus is attached to the
                 boundary annulus.  In particular, edge groups \a roles_[0],
                 \a roles_[1] and \a roles_[2] of the layered solid torus are
                 attached to the vertical, horizontal and diagonal edges of
                 the annulus respectively. */

    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatLST(const NSatLST& cloneMe);
        /**
         * Destroys this structure and its internal components.
         */
        ~NSatLST();

        /**
         * Returns details of the layered solid torus that this block
         * represents.
         *
         * @return details of the layered solid torus.
         */
        const NLayeredSolidTorus* lst() const;
        /**
         * Describes how the layered solid torus is attached to the
         * boundary annulus.
         *
         * The NLayeredSolidTorus class notes describe top-level edge
         * groups 0, 1 and 2 for a layered solid torus.  On the other
         * hand, the NSatAnnulus class notes define vertical, horizontal
         * and diagonal edges on the boundary annulus.
         *
         * Suppose that the permutation returned by this routine is \a r.
         * This indicates that:
         * - edge group \a r[0] is attached to the vertical annulus edges;
         * - edge group \a r[1] is attached to the horizontal annulus edges;
         * - edge group \a r[2] is attached to the diagonal annulus edges.
         *
         * The image \a r[3] will always be 3.
         *
         * @return a description of how the layered solid torus is
         * attached to the boundary annulus.
         */
        NPerm4 roles() const;

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void transform(const NTriangulation* originalTri,
            const NIsomorphism* iso, NTriangulation* newTri);
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (layered solid torus).  This routine is
         * a specific case of NSatBlock::isBlock(); see that routine for
         * further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatLST* isBlockLST(const NSatAnnulus& annulus,
            TetList& avoidTets);

    private:
        /**
         * Constructs a partially initialised block.  The boundary
         * annulus will remain uninitialised, and must be initialised
         * before this block can be used.
         *
         * @param lst details of the layered solid torus.
         * @param roles describes how the layered solid torus is
         * attached to the boundary annulus, as explained in the
         * \a roles_ data member documentation.
         */
        NSatLST(NLayeredSolidTorus* lst, NPerm4 roles);
};

/**
 * A saturated block that is a three-tetrahedron triangular prism.
 *
 * Such a prism may be of major type or of minor type.  In a \e major
 * type prism, the horizontal edges of the boundary annuli are all
 * major (degree three) edges of the prism.  Likewise, in a \e minor
 * type prism, the horizontal boundary edges are all minor (degree two)
 * edges of the prism.  See the NSatAnnulus class notes for a definition
 * of "horizontal" and the NTriSolidTorus class notes for further
 * details regarding "major" and "minor".
 */
class REGINA_API NSatTriPrism : public NSatBlock {
    private:
        bool major_;
            /**< Is this prism of major type or of minor type? */

    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatTriPrism(const NSatTriPrism& cloneMe);

        /**
         * Is this prism of major type or minor type?  See the class
         * notes for further details.
         *
         * Note that this routine cannot be called major(), since on
         * some compilers that name clashes with a macro for isolating
         * major/minor bytes.
         *
         * @return \c true if this prism is of major type, or \c false
         * if it is of minor type.
         */
        bool isMajor() const;

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (triangular prism).  This routine is
         * a specific case of NSatBlock::isBlock(); see that routine for
         * further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatTriPrism* isBlockTriPrism(const NSatAnnulus& annulus,
            TetList& avoidTets);

        /**
         * Inserts a new copy of a triangular prism block into the given
         * triangulation, and returns the corresponding block structure.
         *
         * The given triangulation will not be emptied before the new
         * tetrahedra are inserted.
         *
         * @param tri the triangulation into which the new block should
         * be inserted.
         * @param major \c true if a block of major type should be inserted,
         * or \c false if a block of minor type should be inserted.
         * @return structural details of the newly inserted block.
         */
        static NSatTriPrism* insertBlock(NTriangulation& tri, bool major);

    protected:
        /**
         * Constructs a partially initialised block.  The boundary
         * annuli will remain uninitialised, and must be initialised
         * before this block can be used.
         *
         * @param major \c true if this block is of major type, or
         * \c false if it is of minor type.
         */
        NSatTriPrism(bool major);

    private:
        /**
         * Implements a special case of isBlockTriPrism() to search for
         * a block of major type.  See isBlockTriPrism() for further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatTriPrism* isBlockTriPrismMajor(const NSatAnnulus& annulus,
            TetList& avoidTets);
};

/**
 * A saturated block that is a six-tetrahedron cube.
 *
 * There are several ways of triangulating a cube with six tetrahedra;
 * the specific method used here is illustrated in the diagram below
 * (where the top face of the cube is identified with the bottom).
 *
 * \image html cube.png
 *
 * Note that none of the four tetrahedra that meet the boundary annuli
 * touch each other, and that each of these four boundary tetrahedra
 * meet both central tetrahedra.  Note also that (unlike other
 * triangulations) this cube cannot be split vertically into two
 * triangular prisms.
 */
class REGINA_API NSatCube : public NSatBlock {
    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatCube(const NSatCube& cloneMe);

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (cube).  This routine is a specific case
         * of NSatBlock::isBlock(); see that routine for further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatCube* isBlockCube(const NSatAnnulus& annulus,
            TetList& avoidTets);

        /**
         * Inserts a new copy of a cube block into the given triangulation,
         * and returns the corresponding block structure.
         *
         * The given triangulation will not be emptied before the new
         * tetrahedra are inserted.
         *
         * @param tri the triangulation into which the new block should
         * be inserted.
         * @return structural details of the newly inserted block.
         */
        static NSatCube* insertBlock(NTriangulation& tri);

    protected:
        /**
         * Constructs an uninitialised block.  The boundary annuli
         * must be initialised before this block can be used.
         */
        NSatCube();
};

/**
 * A saturated block that is a reflector strip.
 *
 * A reflector strip is a ring of triangular prisms arranged end-to-end,
 * as illustrated in the diagram below.  The top rectangle of each prism
 * is identified with the bottom in an orientation-reversing fashion
 * (the back edge moves to the front and vice versa), and the prisms
 * are joined in a loop from left to right.  The Seifert fibres run
 * vertically in the diagram, with each saturated boundary annulus shaded
 * at the rear of each prism.
 *
 * \image html reflector.png
 *
 * The effect of a reflector strip is to create a reflector boundary in
 * the base orbifold of the surrounding Seifert fibred space.  Each prism
 * provides a segment of this reflector boundary.
 *
 * A reflector strip may have arbitrary length, and it may also include
 * a twist as the ring of prisms wraps back around to meet itself.  Note
 * that a twisted reflector strip will have a twisted ring of boundary
 * annuli, as described by NSatBlock::twistedBoundary().
 *
 * The \e length of a reflector strip is defined to be the number of
 * prisms that are joined together, or equivalently the number of
 * saturated annuli on the boundary.
 */
class REGINA_API NSatReflectorStrip : public NSatBlock {
    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatReflectorStrip(const NSatReflectorStrip& cloneMe);

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (reflector strip).  This routine is a specific
         * case of NSatBlock::isBlock(); see that routine for further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatReflectorStrip* isBlockReflectorStrip(
            const NSatAnnulus& annulus, TetList& avoidTets);

        /**
         * Inserts a new reflector strip into the given triangulation,
         * and returns the corresponding block structure.
         *
         * The given triangulation will not be emptied before the new
         * tetrahedra are inserted.
         *
         * @param tri the triangulation into which the new block should
         * be inserted.
         * @param length the length of the new reflector strip, i.e.,
         * the number of boundary annuli; this must be strictly positive.
         * @param twisted \c true if the new reflector strip should be twisted
         * (causing its ring of boundary annuli to be twisted also), or
         * \c false if the new strip should not be twisted.
         * @return structural details of the newly inserted block.
         */
        static NSatReflectorStrip* insertBlock(NTriangulation& tri,
            unsigned length, bool twisted);

    protected:
        /**
         * Constructs a partially initialised block of the given length.
         * The boundary annuli will remain uninitialised, and must be
         * initialised before this block can be used.
         *
         * @param length the length of the new reflector strip, i.e.,
         * the number of boundary annuli; this must be strictly positive.
         * @param twisted \c true if the strip should be twisted (giving
         * a twisted ring of boundary annuli), or \c false if not.
         */
        NSatReflectorStrip(unsigned length, bool twisted);
};

/**
 * A degenerate saturated block that is a single tetrahedron wrapped
 * around so that two opposite edges touch.  This forms a degenerate
 * one-tetrahedron solid torus that is pinched along a single meridinal
 * curve.
 *
 * The four faces of this tetrahedron form two boundary annuli, and the
 * tetrahedron is effectively layered onto each boundary annulus.  See
 * the NLayering class notes for more discussion on layerings in general.
 *
 * Although this block is degenerate (the fibres are all pinched
 * together where the opposite edges of the tetrahedron meet), it can be
 * used without problems as long as the entire Seifert fibred space is
 * not formed from degenerate blocks.  In other words, using such blocks
 * is fine as long as they eventually meet a real (non-degenerate) block,
 * which will give room for the fibres to separate so that they are no
 * longer pinched together.
 *
 * The NSatAnnulus class notes describe horizontal and diagonal edges of
 * a saturated annulus.  This block may be one of two types, according
 * to how the tetrahedron is layered onto the boundary annuli.  Either
 * the tetrahedron can be layered over the horizontal edge of each
 * annulus (with the fibres pinched together between the two diagonal
 * edges), or the tetrahedron can be layered over the diagonal edge of
 * each annulus (with the fibres pinched together between the two
 * horizontal edges).
 */
class REGINA_API NSatLayering : public NSatBlock {
    private:
        bool overHorizontal_;
            /**< Do we layer over the horizontal annulus edge, or the
                 diagonal annulus edge? */

    public:
        /**
         * Constructs a clone of the given block structure.
         *
         * @param cloneMe the block structure to clone.
         */
        NSatLayering(const NSatLayering& cloneMe);

        /**
         * Does this describe a layering over the horizontal edge of the
         * boundary annulus, or a layering over the diagonal edge?
         *
         * See the NSatAnnulus class notes for definitions of horizontal
         * and diagonal in this context.
         */
        bool overHorizontal() const;

        virtual NSatBlock* clone() const;
        virtual void adjustSFS(NSFSpace& sfs, bool reflect) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeAbbr(std::ostream& out, bool tex = false) const;

        /**
         * Determines whether the given annulus is a boundary annulus for
         * a block of this type (single layering).  This routine is
         * a specific case of NSatBlock::isBlock(); see that routine for
         * further details.
         *
         * @param annulus the proposed boundary annulus that should form
         * part of the new saturated block.
         * @param avoidTets the list of tetrahedra that should not be
         * considered, and to which any new tetrahedra will be added.
         * @return details of the saturated block if one was found, or
         * \c null if none was found.
         */
        static NSatLayering* isBlockLayering(const NSatAnnulus& annulus,
            TetList& avoidTets);

    protected:
        /**
         * Constructs a partially initialised block.  The boundary
         * annuli will remain uninitialised, and must be initialised
         * before this block can be used.
         *
         * @param overHorizontal \c true if this block describes a
         * layering over the horizontal edge of the boundary annulus, or
         * \c false if it describes a layering over the diagonal edge.
         */
        NSatLayering(bool overHorizontal);
};

/*@}*/

// Inline functions for NSatMobius

inline NSatMobius::NSatMobius(const NSatMobius& cloneMe) : NSatBlock(cloneMe),
        position_(cloneMe.position_) {
}

inline NSatMobius::NSatMobius(int position) : NSatBlock(1),
        position_(position) {
}

inline int NSatMobius::position() const {
    return position_;
}

inline NSatBlock* NSatMobius::clone() const {
    return new NSatMobius(*this);
}

// Inline functions for NSatLST

inline NSatLST::NSatLST(NLayeredSolidTorus* lst, NPerm4 roles) : NSatBlock(1),
        lst_(lst), roles_(roles) {
}

inline const NLayeredSolidTorus* NSatLST::lst() const {
    return lst_;
}

inline NPerm4 NSatLST::roles() const {
    return roles_;
}

inline NSatBlock* NSatLST::clone() const {
    return new NSatLST(*this);
}

// Inline functions for NSatTriPrism

inline NSatTriPrism::NSatTriPrism(const NSatTriPrism& cloneMe) :
        NSatBlock(cloneMe), major_(cloneMe.major_) {
}

inline NSatTriPrism::NSatTriPrism(bool major) : NSatBlock(3), major_(major) {
}

inline bool NSatTriPrism::isMajor() const {
    return major_;
}

inline NSatBlock* NSatTriPrism::clone() const {
    return new NSatTriPrism(*this);
}

inline void NSatTriPrism::writeTextShort(std::ostream& out) const {
    out << "Saturated triangular prism of "
        << (major_ ? "major" : "minor") << " type";
}

inline void NSatTriPrism::writeAbbr(std::ostream& out, bool tex) const {
    if (tex)
        out << "\\triangle";
    else
        out << "Tri";
}

// Inline functions for NSatCube

inline NSatCube::NSatCube(const NSatCube& cloneMe) : NSatBlock(cloneMe) {
}

inline NSatCube::NSatCube() : NSatBlock(4) {
}

inline NSatBlock* NSatCube::clone() const {
    return new NSatCube(*this);
}

inline void NSatCube::writeTextShort(std::ostream& out) const {
    out << "Saturated cube";
}

inline void NSatCube::writeAbbr(std::ostream& out, bool tex) const {
    if (tex)
        out << "\\square";
    else
        out << "Cube";
}

// Inline functions for NSatReflectorStrip

inline NSatReflectorStrip::NSatReflectorStrip(
        const NSatReflectorStrip& cloneMe) : NSatBlock(cloneMe) {
}

inline NSatReflectorStrip::NSatReflectorStrip(unsigned length, bool twisted) :
        NSatBlock(length, twisted) {
}

inline NSatBlock* NSatReflectorStrip::clone() const {
    return new NSatReflectorStrip(*this);
}

inline void NSatReflectorStrip::writeTextShort(std::ostream& out) const {
    out << "Saturated reflector strip of length " << nAnnuli();
    if (twistedBoundary())
        out << " (twisted)";
}

inline void NSatReflectorStrip::writeAbbr(std::ostream& out, bool tex) const {
    if (twistedBoundary()) {
        if (tex)
            out << "\\tilde{\\circledash}_" << nAnnuli();
        else
            out << "Ref~(" << nAnnuli() << ')';
    } else {
        if (tex)
            out << "\\circledash_" << nAnnuli();
        else
            out << "Ref(" << nAnnuli() << ')';
    }
}

// Inline functions for NSatLayering

inline NSatLayering::NSatLayering(const NSatLayering& cloneMe) :
        NSatBlock(cloneMe), overHorizontal_(cloneMe.overHorizontal_) {
}

inline NSatLayering::NSatLayering(bool overHorizontal) :
        NSatBlock(2), overHorizontal_(overHorizontal) {
}

inline bool NSatLayering::overHorizontal() const {
    return overHorizontal_;
}

inline NSatBlock* NSatLayering::clone() const {
    return new NSatLayering(*this);
}

inline void NSatLayering::writeTextShort(std::ostream& out) const {
    out << "Saturated single layering over "
        << (overHorizontal_ ? "horizontal" : "diagonal") << " edge";
}

inline void NSatLayering::writeAbbr(std::ostream& out, bool tex) const {
    if (tex)
        out << "lozenge";
    else
        out << "Layer";
}

} // namespace regina

#endif

