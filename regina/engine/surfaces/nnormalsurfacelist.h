
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

/*! \file surfaces/nnormalsurfacelist.h
 *  \brief Contains a packet representing a collection of normal
 *  surfaces in a 3-manifold.
 */

#ifndef __NNORMALSURFACELIST_H
#ifndef __DOXYGEN
#define __NNORMALSURFACELIST_H
#endif

#include <algorithm>
#include <iterator>
#include <vector>
#include "regina-core.h"
#include "packet/npacket.h"
#include "surfaces/nnormalsurface.h"
#include "surfaces/nsurfaceset.h"
#include "utilities/memutils.h"
#include "utilities/nthread.h"

namespace regina {

class NTriangulation;
class NMatrixInt;
class NProgressManager;
class NProgressMessage;
class NProgressNumber;
class NXMLPacketReader;
class NXMLNormalSurfaceListReader;

/**
 * \weakgroup surfaces
 * @{
 */

/**
 * A packet representing a collection of normal surfaces in a 3-manifold.
 * Such a packet must always be a child packet of the triangulation from
 * which the surfaces were obtained.  If this triangulation changes, the
 * information contained in this packet will become invalid.
 *
 * See the NNormalSurfaceVector class notes for details of what to do
 * when introducing a new flavour of coordinate system.
 *
 * Normal surface lists should be created using the routine enumerate(),
 * which is new as of Regina 3.95.
 *
 * \testpart
 *
 * \todo \feature Allow custom matching equations.
 * \todo \feature Allow enumeration with some coordinates explicitly set
 * to zero.
 * \todo \feature Allow generating only closed surfaces.
 * \todo \feature Generate facets of the solution space representing
 * embedded surfaces.
 */
class REGINA_API NNormalSurfaceList : public NPacket, public NSurfaceSet {
    public:
        static const int packetType;

        /**
         * Represents standard triangle-quadrilateral coordinates for
         * normal surfaces.
         */
        static const int STANDARD;
        /**
         * Represents standard triangle-quadrilateral-octagon coordinates
         * for octagonal almost normal surfaces.
         */
        static const int AN_STANDARD;
        /**
         * Represents quadrilateral coordinates for normal surfaces.
         * For details, see "Normal surface Q-theory", Jeffrey L. Tollefson,
         * Pacific J. Math. 183 (1998), no. 2, 359--374.
         */
        static const int QUAD;
        /**
         * Represents quadrilateral-octagon coordinates for octagonal
         * almost normal surfaces.  For details, see
         * "Quadrilateral-octagon coordinates for almost normal surfaces",
         * Benjamin A. Burton, Experiment. Math. 19 (2010), 285-315.
         */
        static const int AN_QUAD_OCT;

        /**
         * Represents edge weight coordinates for normal surfaces.
         * This flavour is for representation only; surface
         * vectors and lists of this flavour cannot be created.
         */
        static const int EDGE_WEIGHT;
        /**
         * Represents face arc coordinates for normal surfaces.
         * This flavour is for representation only; surface
         * vectors and lists of this flavour cannot be created.
         */
        static const int FACE_ARCS;
        /**
         * Indicates that a list of almost normal surfaces was created
         * using Regina 4.5.1 or earlier, where surfaces with more than
         * one octagon of the same type were stripped out of the final
         * solution set.  As of Regina 4.6 such surfaces are now
         * included in the solution set, since we need them if we
         * wish to enumerate \e all almost normal surfaces (not just
         * the \e vertex almost normal surfaces).
         *
         * This flavour is only used with legacy data files; new vectors
         * and lists of this flavour cannot be created.  The underlying
         * coordinate system is identical to AN_STANDARD.
         */
        static const int AN_LEGACY;
        /**
         * Represents standard triangle-quadrilateral coordinates for
         * transversely oriented normal surfaces.
         */
        static const int ORIENTED;
        /**
         * Represents quadrilateral coordinates for transversely oriented 
         * normal surfaces.
         */
        static const int ORIENTED_QUAD;

        class VectorIterator;

    protected:
        std::vector<NNormalSurface*> surfaces;
            /**< Contains the normal surfaces stored in this packet. */
        int flavour;
            /**< Stores which flavour of coordinate system is being
             *   used by the normal surfaces in this packet. */
        bool embedded;
            /**< Stores whether we are only interested in embedded
             *   normal surfaces. */

    public:
        /**
         * Destroys this list and all the surfaces within.
         */
        virtual ~NNormalSurfaceList();

        /**
         * Enumerates all vertex normal surfaces in the given
         * triangulation using the given flavour of coordinate system.
         * These vertex normal surfaces will be stored in a new normal
         * surface list.  Their representations will
         * use the smallest possible integer coordinates.
         * The option is offered to find only embedded normal surfaces
         * or to also include immersed and singular normal surfaces.
         *
         * The normal surface list that is created will be inserted as the
         * last child of the given triangulation.  This triangulation \b must
         * remain the parent of this normal surface list, and must not
         * change while this normal surface list remains in existence.
         *
         * If a progress manager is passed, the normal surface
         * enumeration will take place in a new thread and this routine
         * will return immediately.  The NProgress object assigned to
         * this progress manager is guaranteed to be of the class
         * NProgressNumber.
         *
         * If no progress manager is passed, the enumeration will run
         * in the current thread and this routine will return only when
         * the enumeration is complete.  Note that this enumeration can
         * be extremely slow for larger triangulations.
         *
         * \todo \feature Allow picking up the first ``interesting'' surface
         * and bailing en route.
         * \todo \featurelong Determine the faces of the normal solution
         * space.
         * \todo \featurelong Allow either subsets of normal surface
         * lists or allow deletion of surfaces from lists.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @param newFlavour the flavour of coordinate system to be used;
         * this must be one of the predefined coordinate system
         * constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced; this defaults to
         * \c true.
         * @param manager a progress manager through which progress will
         * be reported, or 0 if no progress reporting is required.  If
         * non-zero, \a manager must point to a progress manager for
         * which NProgressManager::isStarted() is still \c false.
         * @return the newly created normal surface list.  Note that if
         * a progress manager is passed then this list may not be completely
         * filled when this routine returns.  If a progress manager is
         * passed and a new thread could not be started, this routine
         * returns 0 (and no normal surface list is created).
         */
        static NNormalSurfaceList* enumerate(NTriangulation* owner,
            int newFlavour, bool embeddedOnly = true,
            NProgressManager* manager = 0);

        /**
         * Enumerates all fundamental normal surfaces in the given
         * triangulation using the given flavour of coordinate system,
         * using the primal Hilbert basis algorithm.
         * These fundamental normal surfaces will be stored in a new normal
         * surface list.  The option is offered to find only embedded
         * normal surfaces or to also include immersed and singular
         * normal surfaces.
         *
         * The primal algorithm is the recommended method for
         * enumerating fundamental normal surfaces, although other
         * algorithms are made available in this class also.
         * For full details of the procedure, see "Fundamental normal
         * surfaces and the enumeration of Hilbert bases", Burton,
         * arXiv:1111.7055, Nov 2011.
         *
         * The normal surface list that is created will be inserted as the
         * last child of the given triangulation.  This triangulation \b must
         * remain the parent of this normal surface list, and must not
         * change while this normal surface list remains in existence.
         *
         * The first step of the primal algorithm is to enumerate all
         * vertex normal surfaces.  If you have already done this, you
         * may pass the list of vertex normal surfaces as the (optional)
         * parameter \a vtxSurfaces.
         *
         * If a progress manager is passed, the normal surface
         * enumeration will take place in a new thread and this routine
         * will return immediately.  The NProgress object assigned to
         * this progress manager is guaranteed to be of the class
         * NProgressMessage.
         *
         * If no progress manager is passed, the enumeration will run
         * in the current thread and this routine will return only when
         * the enumeration is complete.  Note that this enumeration can
         * be extremely slow for larger triangulations.
         *
         * \pre If non-zero, the argument \a vtxSurfaces is precisely
         * the set of all vertex normal surfaces in the given triangulation,
         * enumerated using the same coordinate system and embedded-only
         * constraints as given here.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @param newFlavour the flavour of coordinate system to be used;
         * this must be one of the predefined coordinate system
         * constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced; this defaults to
         * \c true.
         * @param vtxSurfaces the set of all \e vertex normal surfaces
         * as enumerated under the same coordinate system and
         * constraints as given here; this may be 0 if unknown.
         * @param manager a progress manager through which progress will
         * be reported, or 0 if no progress reporting is required.  If
         * non-zero, \a manager must point to a progress manager for
         * which NProgressManager::isStarted() is still \c false.
         * @return the newly created normal surface list.  Note that if
         * a progress manager is passed then this list may not be completely
         * filled when this routine returns.  If a progress manager is
         * passed and a new thread could not be started, this routine
         * returns 0 (and no normal surface list is created).
         */
        static NNormalSurfaceList* enumerateFundPrimal(
            NTriangulation* owner, int newFlavour, bool embeddedOnly = true,
            NNormalSurfaceList* vtxSurfaces = 0,
            NProgressManager* manager = 0);

        /**
         * Enumerates all fundamental normal surfaces in the given
         * triangulation using the given flavour of coordinate system,
         * using the dual Hilbert basis algorithm.
         * These fundamental normal surfaces will be stored in a new normal
         * surface list.  The option is offered to find only embedded
         * normal surfaces or to also include immersed and singular
         * normal surfaces.
         *
         * The dual algorithm is fast but its performance is highly variable;
         * for this reason the primal algorithm is recommended instead.
         * For full details of both procedures, see
         * "Fundamental normal surfaces and the enumeration of Hilbert bases",
         * Burton, arXiv:1111.7055, Nov 2011.
         *
         * The normal surface list that is created will be inserted as the
         * last child of the given triangulation.  This triangulation \b must
         * remain the parent of this normal surface list, and must not
         * change while this normal surface list remains in existence.
         *
         * If a progress manager is passed, the normal surface
         * enumeration will take place in a new thread and this routine
         * will return immediately.  The NProgress object assigned to
         * this progress manager is guaranteed to be of the class
         * NProgressNumber.
         *
         * If no progress manager is passed, the enumeration will run
         * in the current thread and this routine will return only when
         * the enumeration is complete.  Note that this enumeration can
         * be extremely slow for larger triangulations.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @param newFlavour the flavour of coordinate system to be used;
         * this must be one of the predefined coordinate system
         * constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced; this defaults to
         * \c true.
         * @param manager a progress manager through which progress will
         * be reported, or 0 if no progress reporting is required.  If
         * non-zero, \a manager must point to a progress manager for
         * which NProgressManager::isStarted() is still \c false.
         * @return the newly created normal surface list.  Note that if
         * a progress manager is passed then this list may not be completely
         * filled when this routine returns.  If a progress manager is
         * passed and a new thread could not be started, this routine
         * returns 0 (and no normal surface list is created).
         */
        static NNormalSurfaceList* enumerateFundDual(
            NTriangulation* owner, int newFlavour, bool embeddedOnly = true,
            NProgressManager* manager = 0);

        /**
         * Uses a slow-but-direct procedure to enumerate all embedded
         * vertex normal surfaces in standard (tri-quad) coordinates
         * within the given triangulation.
         *
         * The standard enumerate() routine will choose the fastest
         * available algorithm for enumerating vertex normal surfaces.
         * In particular, when enumerating embedded vertex normal
         * surfaces in standard (tri-quad) coordinates, it will often take
         * a two-step approach: (i) enumerate vertex normal surfaces in
         * \e quadrilateral space; (ii) convert the quadrilateral space
         * solution set to a standard tri-quad space solution set.  This
         * two-step procedure is typically \e much faster than enumerating
         * solutions in standard coordinates directly.  For details on
         * this procedure see "Converting between quadrilateral and
         * standard solution sets in normal surface theory",
         * Benjamin A. Burton, Algebr. Geom. Topol. 9 (2009), 2121-2174.
         *
         * This routine allows the user to force a direct enumeration in
         * standard space, \e without going via quadrilateral space.
         * The algorithm used is the souped-up double description method
         * in standard coordinates as described in "Optimizing the
         * double description method for normal surface enumeration",
         * Benjamin A. Burton, Math. Comp. 79 (2010), 453-484.
         *
         * Aside from the underlying algorithm, the behaviour of this
         * routine is identical to enumerate().  See enumerate() for
         * details regarding preconditions, postconditions, ownership
         * and so on.
         *
         * Unlike enumerate(), this routine does not support progress
         * management and does not support running in a separate thread.
         *
         * \warning
         * Users will generally not want to call this routine, since it
         * is often much slower than enumerate() and it gives precisely
         * the same results.  This routine is provided mainly for interest's
         * sake, and to allow comparisons between different algorithms.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @return the newly created normal surface list.
         */
        static NNormalSurfaceList* enumerateStandardDirect(
            NTriangulation* owner);

        /**
         * Uses a slow-but-direct procedure to enumerate all embedded
         * vertex almost normal surfaces in standard (tri-quad-oct)
         * coordinates within the given triangulation.
         *
         * This routine is the almost normal analogue to the
         * enumerateStandardDirect() enumeration routine; see the
         * enumerateStandardDirect() documentation for further information.
         *
         * \warning
         * Users will generally not want to call this routine, since it
         * is often much slower than enumerate() and it gives precisely
         * the same results.  This routine is provided mainly for interest's
         * sake, and to allow comparisons between different algorithms.
         *
         * @param owner the triangulation upon which this list of
         * almost normal surfaces will be based.
         * @return the newly created surface list.
         */
        static NNormalSurfaceList* enumerateStandardANDirect(
            NTriangulation* owner);

        /**
         * Uses an extremely slow procedure to enumerate all embedded
         * fundamental surfaces in the given triangulation,
         * by running Normaliz over the full (and typically very large)
         * solution cone, and only enforcing embedded constraints (such as
         * the quadrilateral constraints) afterwards.
         *
         * Aside from the underlying algorithm, the behaviour of this
         * routine is identical to enumerateFundPrimal() and
         * enumerateFundDual().  See those routines for details
         * regarding preconditions, postconditions, ownership and so on.
         *
         * Unlike enumerateFundPrimal() and enumerateFundDual(), this
         * routine does not support progress management, and does not
         * support running in a separate thread.
         *
         * \warning
         * Users will generally not want to call this routine, since it
         * is typically much slower than either enumerateFundPrimal()
         * or enumerateFundDual(), and it gives precisely the same results.
         * This routine is provided mainly for interest's sake, and to allow
         * comparisons between different algorithms.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @param newFlavour the flavour of coordinate system to be used;
         * this must be one of the predefined coordinate system
         * constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced; this defaults to \c true.
         * @return the newly created normal surface list.
         */
        static NNormalSurfaceList* enumerateFundFullCone(
            NTriangulation* owner, int newFlavour, bool embeddedOnly = true);

        /**
         * Uses an extremely slow modified Contejean-Devie procedure to
         * enumerate all embedded fundamental surfaces in the given
         * triangulation.  For details of the modifications, see
         * "Fundamental normal surfaces and the enumeration of Hilbert bases",
         * Burton, arXiv:1111.7055, Nov 2011.
         *
         * Aside from the underlying algorithm, the behaviour of this
         * routine is identical to enumerateFundPrimal() and
         * enumerateFundDual().  See those routines for details
         * regarding preconditions, postconditions, ownership and so on.
         *
         * Unlike enumerateFundPrimal() and enumerateFundDual(), this
         * routine does not support progress management, and does not
         * support running in a separate thread.
         *
         * \warning
         * Users will generally not want to call this routine, since it
         * is typically much slower than either enumerateFundPrimal()
         * or enumerateFundDual(), and it gives precisely the same results.
         * This routine is provided mainly for interest's sake, and to allow
         * comparisons between different algorithms.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces will be based.
         * @param newFlavour the flavour of coordinate system to be used;
         * this must be one of the predefined coordinate system
         * constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced; this defaults to \c true.
         * @return the newly created normal surface list.
         */
        static NNormalSurfaceList* enumerateFundCD(
            NTriangulation* owner, int newFlavour, bool embeddedOnly = true);

        virtual int getFlavour() const;
        virtual bool allowsAlmostNormal() const;
        virtual bool allowsSpun() const;
        virtual bool allowsOriented() const;
        virtual bool isEmbeddedOnly() const;
        virtual NTriangulation* getTriangulation() const;
        virtual unsigned long getNumberOfSurfaces() const;
        virtual const NNormalSurface* getSurface(unsigned long index) const;
        virtual ShareableObject* getShareableObject();

        virtual int getPacketType() const;
        virtual std::string getPacketTypeName() const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeTextLong(std::ostream& out) const;
        static NXMLPacketReader* getXMLReader(NPacket* parent);
        virtual void writePacket(NFile& out) const;
        static NNormalSurfaceList* readPacket(NFile& in, NPacket* parent);
        virtual bool dependsOnParent() const;

        /**
         * Converts the set of all embedded vertex normal surfaces in
         * quadrilateral space to the set of all embedded vertex normal
         * surfaces in standard (tri-quad) space.  The initial list in
         * quadrilateral space is taken to be this normal surface list;
         * the final list in standard space will be inserted as a new
         * child packet of the underlying triangulation (specifically, as
         * the final child).  As a convenience, the final list will also
         * be returned from this routine.
         *
         * This routine can only be used with normal surfaces, not almost
         * normal surfaces.  For almost normal surfaces, see the similar
         * routine quadOctToStandardAN().
         *
         * This procedure is available for any triangulation whose vertex
         * links are all spheres and/or discs, and is \e much faster than
         * enumerating surfaces directly in standard tri-quad coordinates.
         * The underlying algorithm is described in detail in "Converting
         * between quadrilateral and standard solution sets in normal
         * surface theory", Benjamin A. Burton, Algebr. Geom. Topol. 9 (2009),
         * 2121-2174.
         *
         * Typically users do not need to call this routine directly,
         * since the standard enumerate() routine will use it implicitly
         * where possible.  That is, when asked for standard vertex surfaces,
         * enumerate() will first find all \e quadrilateral vertex surfaces
         * and then use this procedure to convert them to standard vertex
         * surfaces; this is generally orders of magnitude faster than
         * enumerating surfaces directly in standard coordinates.
         *
         * Nevertheless, this standalone routine is provided as a convenience
         * for users who already have a set of quadrilateral vertex surfaces,
         * and who simply wish to convert them to a set of standard
         * vertex surfaces without the cost of implicitly enumerating the
         * quadrilateral vertex surfaces again.
         *
         * It should be noted that this routine does \e not simply convert
         * vectors from one form to another; instead it converts a full
         * solution set of vertex surfaces in quadrilateral coordinates to a
         * full solution set of vertex surfaces in standard coordinates.
         * Typically there are many more vertex surfaces in standard
         * coordinates (all of which this routine will find).
         *
         * This routine will run some very basic sanity checks before
         * starting.  Specifically, it will check the validity and vertex
         * links of the underlying triangulation, and will verify that
         * the coordinate flavour and embedded-only flag are set to
         * NNormalSurfaceList::QUAD and \c true respectively.  If any of
         * these checks fails, this routine will do nothing and return 0.
         *
         * \pre The underlying triangulation (the parent packet of this
         * normal surface list) is valid, and the link of every vertex
         * is either a sphere or a disc.
         * \pre This normal surface list is precisely the set of all
         * embedded vertex normal surfaces in quadrilateral space; no more,
         * no less.  Moreover, these vectors are stored using quadrilateral
         * coordinates.  Typically this means that it was obtained through
         * enumerate(), with the flavour set to NNormalSurfaceList::QUAD and
         * with \a embeddedOnly set to \c true.
         *
         * @return a full list of vertex normal surfaces in standard (tri-quad)
         * coordinates, or 0 if any of the basic sanity checks failed.
         */
        NNormalSurfaceList* quadToStandard() const;

        /**
         * Converts the set of all embedded vertex almost normal surfaces in
         * quadrilateral-octagon space to the set of all embedded vertex
         * almost normal surfaces in the standard tri-quad-oct space.
         *
         * This routine is the almost normal analogue to the
         * quadToStandard() conversion routine; see the quadToStandard()
         * documentation for further information.
         *
         * \pre The underlying triangulation (the parent packet of this
         * normal surface list) is valid, and the link of every vertex
         * is either a sphere or a disc.
         * \pre This surface list is precisely the set of all embedded vertex
         * almost normal surfaces in quadrilateral-octagon space; no more,
         * no less.  Moreover, these vectors are stored using
         * quadrilateral-octagon coordinates.  Typically this means that it
         * was obtained through enumerate(), with the flavour set to
         * NNormalSurfaceList::AN_QUAD_OCT and with \a embeddedOnly set to
         * \c true.
         *
         * @return a full list of vertex almost normal surfaces in standard
         * tri-quad-oct coordinates, or 0 if any of the basic sanity checks
         * failed.
         */
        NNormalSurfaceList* quadOctToStandardAN() const;

        /**
         * Converts the set of all embedded vertex normal surfaces in
         * standard (tri-quad) space to the set of all embedded vertex
         * normal surfaces in quadrilateral space.  The initial list in
         * standard space is taken to be this normal surface list;
         * the final list in quadrilateral space will be inserted as a new
         * child packet of the underlying triangulation (specifically, as
         * the final child).  As a convenience, the final list will also
         * be returned from this routine.
         *
         * This routine can only be used with normal surfaces, not almost
         * normal surfaces.  For almost normal surfaces, see the similar
         * routine standardANToQuadOct().
         *
         * This procedure is available for any triangulation whose vertex
         * links are all spheres and/or discs.  The underlying algorithm
         * is described in detail in "Converting between quadrilateral and
         * standard solution sets in normal surface theory",
         * Benjamin A. Burton, Algebr. Geom. Topol. 9 (2009), 2121-2174.
         *
         * It should be noted that this routine does \e not simply convert
         * vectors from one form to another; instead it converts a full
         * solution set of vertex surfaces in standard coordinates to a
         * full solution set of vertex surfaces in quadrilateral coordinates.
         * Typically there are far fewer vertex surfaces in quadrilateral
         * coordinates (all of which this routine will find).
         *
         * This routine will run some very basic sanity checks before
         * starting.  Specifically, it will check the validity and vertex
         * links of the underlying triangulation, and will verify that
         * the coordinate flavour and embedded-only flag are set to
         * NNormalSurfaceList::STANDARD and \c true respectively.  If any of
         * these checks fails, this routine will do nothing and return 0.
         *
         * \pre The underlying triangulation (the parent packet of this
         * normal surface list) is valid, and the link of every vertex
         * is either a sphere or a disc.
         * \pre This normal surface list is precisely the set of all
         * embedded vertex normal surfaces in standard (tri-quad) space;
         * no more, no less.  Moreover, these vectors are stored using
         * standard coordinates.  Typically this means that this list was
         * obtained through enumerate(), with the flavour set to
         * NNormalSurfaceList::STANDARD and with \a embeddedOnly set to
         * \c true.
         *
         * @return a full list of vertex normal surfaces in quadrilateral
         * coordinates, or 0 if any of the basic sanity checks failed.
         */
        NNormalSurfaceList* standardToQuad() const;

        /**
         * Converts the set of all embedded vertex almost normal surfaces in
         * standard tri-quad-oct space to the set of all embedded vertex
         * almost normal surfaces in the smaller quadrilateral-octagon space.
         *
         * This routine is the almost normal analogue to the
         * standardToQuad() conversion routine; see the standardToQuad()
         * documentation for further information.
         *
         * \pre The underlying triangulation (the parent packet of this
         * normal surface list) is valid, and the link of every vertex
         * is either a sphere or a disc.
         * \pre This normal surface list is precisely the set of all
         * embedded vertex almost normal surfaces in standard tri-quad-oct
         * space; no more, no less.  Typically this means that it was obtained
         * through enumerate(), with the flavour set to
         * NNormalSurfaceList::AN_STANDARD and with \a embeddedOnly set
         * to \c true.
         *
         * @return a full list of vertex almost normal surfaces in
         * quadrilateral-octagon coordinates, or 0 if any of the basic
         * sanity checks failed.
         */
        NNormalSurfaceList* standardANToQuadOct() const;

        /**
         * Creates a new list filled with the surfaces from this list
         * that have at least one locally compatible partner.
         * In other words, a surface \a S from this list will be placed
         * in the new list if and only if there is some other surface \a T
         * in this list for which \a S and \a T are locally compatible.
         * See NNormalSurface::locallyCompatible() for further details on
         * compatibility testing.
         *
         * The new list will be inserted as a new child packet of the
         * underlying triangulation (specifically, as the final child).  As a
         * convenience, the new list will also be returned from this routine.
         *
         * This original list is not altered in any way.  Likewise,
         * the surfaces in the new list are deep copies of the originals
         * (so they can be altered without affecting the original surfaces).
         *
         * \pre This list contains only embedded normal surfaces.  More
         * precisely, isEmbeddedOnly() must return \c true.
         *
         * \warning If this list contains a vertex link (plus at least
         * one other surface), then the new list will be identical to
         * the old (i.e., every surface will be copied across).
         *
         * @return the new list, which will also have been inserted as
         * a new child packet of the underlying triangulation.
         */
        NNormalSurfaceList* filterForLocallyCompatiblePairs() const;

        /**
         * Creates a new list filled with the surfaces from this list
         * that have at least one disjoint partner.
         * In other words, a surface \a S from this list will be placed
         * in the new list if and only if there is some other surface \a T
         * in this list for which \a S and \a T can be made to intersect
         * nowhere at all, without changing either normal isotopy class.
         * See NNormalSurface::disjoint() for further details on disjointness
         * testing.
         *
         * This routine cannot deal with empty, disconnected or
         * non-compact surfaces.  Such surfaces will be silently
         * ignored, and will not be used in any disjointness tests (in
         * particular, they will never be considered as a "disjoint partner"
         * for any other surface).
         *
         * The new list will be inserted as a new child packet of the
         * underlying triangulation (specifically, as the final child).  As a
         * convenience, the new list will also be returned from this routine.
         *
         * This original list is not altered in any way.  Likewise,
         * the surfaces in the new list are deep copies of the originals
         * (so they can be altered without affecting the original surfaces).
         *
         * \pre This list contains only embedded normal surfaces.  More
         * precisely, isEmbeddedOnly() must return \c true.
         * \pre All surfaces within this list are stored using the same
         * flavour of coordinate system (i.e., the same subclass of
         * NNormalSurfaceVector).
         *
         * \warning If this list contains a vertex link (plus at least
         * one other surface), then the new list will be identical to
         * the old (i.e., every surface will be copied across).
         *
         * \todo Deal properly with surfaces that are too large to handle.
         *
         * @return the new list, which will also have been inserted as
         * a new child packet of the underlying triangulation.
         */
        NNormalSurfaceList* filterForDisjointPairs() const;

        /**
         * Creates a new list filled with only the surfaces from this list
         * that "might" represent two-sided incompressible surfaces.
         * More precisely, we consider all two-sided surfaces in this list,
         * as well as the two-sided double covers of all one-sided surfaces
         * in this list (see below for details on how one-sided surfaces
         * are handled).  Each of these surfaces is examined using
         * relatively fast heuristic tests for incompressibility.  Any
         * surface that is definitely \e not incompressible is thrown
         * away, and all other surfaces are placed in the new list.
         *
         * Therefore, it is guaranteed that every incompressible surface
         * from the old list will be placed in the new list.  However,
         * it is not known whether any given surface in the new list is
         * indeed incompressible.
         *
         * See NNormalSurface::isIncompressible() for the definition of
         * incompressibility that is used here.  Note in particular that
         * spheres are \e never considered incompressible.
         *
         * As indicated above, this filter works exclusively with two-sided
         * surfaces.  If a surface in this list is one-sided, the heuristic
         * incompressibility tests will be run on its two-sided double cover.
         * Nevertheless, if the tests pass, the original one-sided surface
         * (not the double cover) will be added to the new list.
         *
         * The new list will be inserted as a new child packet of the
         * underlying triangulation (specifically, as the final child).  As a
         * convenience, the new list will also be returned from this routine.
         *
         * This original list is not altered in any way.  Likewise,
         * the surfaces in the new list are deep copies of the originals
         * (so they can be altered without affecting the original surfaces).
         *
         * Currently the heuristic tests include (i) throwing away
         * all vertex links and thin edge links, and then
         * (ii) cutting along the remaining surfaces and running
         * NTriangulation::hasSimpleCompressingDisc() on the resulting
         * bounded triangulations.  For more details on these tests
         * see "The Weber-Seifert dodecahedral space is non-Haken",
         * Benjamin A. Burton, J. Hyam Rubinstein and Stephan Tillmann,
         * Trans. Amer. Math. Soc. 364:2 (2012), pp. 911-932.
         *
         * \pre The underlying 3-manifold triangulation is valid and closed.
         * In particular, it has no ideal vertices.
         * \pre This list contains only embedded normal surfaces.  More
         * precisely, isEmbeddedOnly() must return \c true.
         * \pre This list contains only compact, connected normal surfaces.
         * \pre No surfaces in this list contain any octagonal discs.
         *
         * \warning The behaviour of this routine is subject to change
         * in future versions of Regina, since additional tests may be
         * added to improve the power of this filtering.
         *
         * \todo Add progress tracking.
         *
         * @return the new list, which will also have been inserted as
         * a new child packet of the underlying triangulation.
         */
        NNormalSurfaceList* filterForPotentiallyIncompressible() const;

        /**
         * Returns a newly created matrix containing the matching
         * equations that were used to create this normal surface list.
         * The destruction of this matrix is the responsibility of the
         * caller of this routine.  Multiple calls to this routine will
         * result in the construction of multiple matrices.  This
         * routine in fact merely calls makeMatchingEquations() with the
         * appropriate parameters.
         *
         * The format of the matrix is identical to that returned by
         * makeMatchingEquations().
         * 
         * @return the matching equations used to create this normal
         * surface list.
         */
        NMatrixInt* recreateMatchingEquations() const;

        /**
         * An iterator that gives access to the raw vectors for surfaces in
         * this list, pointing to the beginning of this surface list.
         *
         * \ifacespython Not present.
         *
         * @return an iterator at the beginning of this surface list.
         */
        VectorIterator beginVectors() const;

        /**
         * An iterator that gives access to the raw vectors for surfaces in
         * this list, pointing past the end of this surface list.
         * This iterator is not dereferenceable.
         *
         * \ifacespython Not present.
         *
         * @return an iterator past the end of this surface list.
         */
        VectorIterator endVectors() const;

        /**
         * A bidirectional iterator that runs through the raw vectors for
         * surfaces in this list.
         *
         * \ifacespython Not present.
         */
        class VectorIterator : public std::iterator<
                std::bidirectional_iterator_tag, const NNormalSurfaceVector*> {
            private:
                std::vector<NNormalSurface*>::const_iterator it_;
                    /**< An iterator into the underlying list of surfaces. */

            public:
                /**
                 * Creates a new uninitialised iterator.
                 */
                VectorIterator();

                /**
                 * Creates a copy of the given iterator.
                 *
                 * @param cloneMe the iterator to clone.
                 */
                VectorIterator(const VectorIterator& cloneMe);

                /**
                 * Makes this a copy of the given iterator.
                 *
                 * @param cloneMe the iterator to clone.
                 * @return a reference to this iterator.
                 */
                VectorIterator& operator = (const VectorIterator& cloneMe);

                /**
                 * Compares this with the given operator for equality.
                 *
                 * @param other the iterator to compare this with.
                 * @return \c true if the iterators point to the same
                 * element of the same normal surface list, or \c false
                 * if they do not.
                 */
                bool operator == (const VectorIterator& other) const;

                /**
                 * Compares this with the given operator for inequality.
                 *
                 * @param other the iterator to compare this with.
                 * @return \c false if the iterators point to the same
                 * element of the same normal surface list, or \c true
                 * if they do not.
                 */
                bool operator != (const VectorIterator& other) const;

                /**
                 * Returns the raw vector for the normal surface that this
                 * iterator is currently pointing to.
                 *
                 * \pre This iterator is dereferenceable (in particular,
                 * it is not past-the-end).
                 *
                 * @return the corresponding normal surface vector.
                 */
                const NNormalSurfaceVector* operator *() const;

                /**
                 * The preincrement operator.
                 *
                 * @return a reference to this iterator after the increment.
                 */
                VectorIterator& operator ++();

                /**
                 * The postincrement operator.
                 *
                 * @return a copy of this iterator before the
                 * increment took place.
                 */
                VectorIterator operator ++(int);

                /**
                 * The predecrement operator.
                 *
                 * @return a reference to this iterator after the decrement.
                 */
                VectorIterator& operator --();

                /**
                 * The postdecrement operator.
                 *
                 * @return a copy of this iterator before the
                 * decrement took place.
                 */
                VectorIterator operator --(int);

            private:
                /**
                 * Initialise a new vector iterator using an iterator for
                 * the internal list of normal surfaces.
                 */
                VectorIterator(
                    const std::vector<NNormalSurface*>::const_iterator& i);

            friend class NNormalSurfaceList;
        };

    protected:
        /**
         * Creates a new normal surface list performing no initialisation
         * whatsoever other than property initialisation.
         */
        NNormalSurfaceList();

        virtual NPacket* internalClonePacket(NPacket* parent) const;
        virtual void writeXMLPacketData(std::ostream& out) const;

        /**
         * An output iterator used to insert surfaces into an
         * NNormalSurfaceList.
         *
         * Objects of type <tt>NNormalSurface*</tt> and
         * <tt>NNormalSurfaceVector*</tt> can be assigned to this
         * iterator.  In the latter case, a surrounding NNormalSurface
         * will be automatically created.
         *
         * \warning The behaviour of this class has changed!
         * As of Regina 4.6, this class happily inserts every surface or
         * vector that it is given.  In previous versions it checked
         * almost normal surface vectors for multiple octagonal discs;
         * this check has been removed to support conversions between
         * quad-oct space and standard almost normal space, and to support
         * the enumeration of \e all almost normal surfaces (as opposed to
         * just vertex surfaces).  Such checks are now left to the user
         * interface (and indeed are now optional, at the user's discretion).
         */
        struct SurfaceInserter : public std::iterator<
                std::output_iterator_tag, NNormalSurfaceVector*> {
            NNormalSurfaceList* list;
                /**< The list into which surfaces will be inserted. */
            NTriangulation* owner;
                /**< The triangulation in which the surfaces to be
                 *   inserted are contained. */

            /**
             * Creates a new uninitialised output iterator.
             *
             * \warning This iterator must not be used until its
             * surface list and triangulation have been initialised.
             */
            SurfaceInserter();
            /**
             * Creates a new output iterator.  The member variables of
             * this iterator will be initialised according to the
             * parameters passed to this constructor.
             *
             * @param newList the list into which surfaces will be
             * inserted.
             * @param newOwner the triangulation in which the surfaces
             * to be inserted are contained.
             */
            SurfaceInserter(NNormalSurfaceList& newList,
                NTriangulation* newOwner);
            /**
             * Creates a new output iterator that is a clone of the
             * given iterator.
             *
             * @param cloneMe the output iterator to clone.
             */
            SurfaceInserter(const SurfaceInserter& cloneMe);

            /**
             * Sets this iterator to be a clone of the given output iterator.
             *
             * @param cloneMe the output iterator to clone.
             * @return this output iterator.
             */
            SurfaceInserter& operator =(const SurfaceInserter& cloneMe);

            /**
             * Appends a normal surface to the end of the appropriate
             * surface list.
             *
             * The given surface will be deallocated with the other
             * surfaces in this list when the list is eventually destroyed.
             *
             * @param surface the normal surface to insert.
             * @return this output iterator.
             */
            SurfaceInserter& operator =(NNormalSurface* surface);
            /**
             * Appends the normal surface corresponding to the given
             * vector to the end of the appropriate surface list.
             *
             * The given vector will be owned by the newly created
             * normal surface and will be deallocated with the other
             * surfaces in this list when the list is eventually destroyed.
             *
             * \warning The behaviour of this routine has changed!
             * As of Regina 4.6, this routine no longer checks for
             * multiple octagonal discs.  See the SurfaceInserter
             * class notes for details.
             *
             * @param vector the vector of the normal surface to insert.
             * @return this output iterator.
             */
            SurfaceInserter& operator =(NNormalSurfaceVector* vector);

            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            SurfaceInserter& operator *();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            SurfaceInserter& operator ++();
            /**
             * Returns a reference to this output iterator.
             *
             * @return this output iterator.
             */
            SurfaceInserter& operator ++(int);
        };

    private:
        /**
         * A helper class containing constants, typedefs and operations
         * for working with normal (as opposed to almost normal) surfaces.
         *
         * This class and its partner AlmostNormalSpec can be used to
         * write generic template code that works with both normal
         * \e and almost normal surfaces.
         *
         * The full definition of this class is in the file normalspec.tcc .
         */
        struct NormalSpec;

        /**
         * A helper class containing constants, typedefs and operations
         * for working with almost normal (as opposed to normal) surfaces.
         *
         * This class and its partner NormalSpec can be used to
         * write generic template code that works with both normal
         * \e and almost normal surfaces.
         *
         * The full definition of this class is in the file normalspec.tcc .
         */
        struct AlmostNormalSpec;

        /**
         * Creates an empty list of normal surfaces with the given
         * parameters.
         *
         * @param newFlavour the flavour of coordinate system to be used
         * for filling this list; this must be one of the predefined
         * coordinate system constants in NNormalSurfaceList.
         * @param embeddedOnly \c true if only embedded normal surfaces
         * are to be produced, or \c false if immersed and singular
         * normal surfaces are also to be produced.
         */
        NNormalSurfaceList(int newFlavour, bool embeddedOnly);

        /**
         * Enumerates all embedded vertex surfaces in (standard normal
         * or standard almost normal) coordinates by first enumerating
         * surfaces in (quad or quad-oct) coordinates and then performing
         * a conversion routine.  See quadToStandard() and
         * quadOctToStandardAN() for further information on this procedure.
         *
         * The resulting surfaces in standard coordinates will be
         * inserted directly into this list.
         *
         * This routine is designed to work with surface lists that are
         * still under construction.  As such, it ignores the packet
         * tree completely.  The parent packet is ignored (and not changed);
         * instead the underlying triangulation is passed explicitly as
         * the argument \a owner.
         *
         * This enumeration will always take place in the current thread.
         * A numeric progress watcher may be passed for progress reporting,
         * in which case this routine will poll for cancellation requests
         * accordingly.
         *
         * \pre The template argument \a Variant is either NormalSpec
         * or AlmostNormalSpec, according to whether we are working with
         * normal or almost normal surfaces.
         * \pre The flavour for this surface list is set to
         * NNormalSurfaceList::STANDARD if we are working with normal
         * surfaces, or NNormalSurfaceList::AN_STANDARD if we are working
         * with almost normal surfaces.  Moreover, the embedded-only flag
         * is set to \c true.
         * \pre The given triangulation is valid, and the link of every
         * vertex is either a sphere or a disc.
         *
         * @param owner the triangulation upon which this list of normal
         * surfaces is to be based.
         * @param progress a numeric progress watcher through which
         * progress will be reported, or 0 if no progress reporting is
         * required.  If a progress watcher is passed, its expected
         * total will be increased immediately by some number of steps
         * and the completed total will be increased gradually by this
         * same number.  NProgress::setFinished() will \e not be called.
         */
        template <class Variant>
        void enumerateStandardViaReduced(NTriangulation* owner,
            NProgressNumber* progress = 0);

        /**
         * Converts a set of embedded vertex normal surfaces in
         * (quad or quad-oct) space to a set of embedded vertex normal
         * surfaces in (standard normal or standard almost normal) space.
         * The original (quad or quad-oct) space surfaces are passed in
         * the argument \a quadList, and the resulting (standard normal
         * or standard almost normal) space surfaces will be inserted
         * directly into this list.
         *
         * See quadToStandard() and quadOctToStandardAN() for full details
         * and preconditions for this procedure.
         *
         * This routine is designed to work with surface lists that are
         * still under construction.  As such, it ignores the packet
         * tree completely.  The parent packet is ignored (and not changed);
         * instead the underlying triangulation is passed explicitly as
         * the argument \a owner.
         *
         * Although this routine takes a vector of non-const pointers, it
         * guarantees not to modify (or destroy) any of the contents.
         *
         * Although this is a template function, it is a private
         * template function and so is only defined in the one <tt>.cpp</tt>
         * file that needs it.
         *
         * \pre The template argument \a Variant is either NormalSpec
         * or AlmostNormalSpec, according to whether we are doing the
         * work for quadToStandard() or quadOctToStandardAN() respectively.
         * \pre The flavour for this surface list is set to
         * NNormalSurfaceList::STANDARD or NNormalSurfaceList::AN_STANDARD,
         * according to whether we are doing the work for quadToStandard()
         * or quadOctToStandardAN() respectively, and the embedded-only
         * flag is set to \c true.
         * \pre The given triangulation is valid and non-empty, and the link
         * of every vertex is either a sphere or a disc.
         *
         * @param owner the triangulation upon which this list of
         * surfaces is to be based.
         * @param reducedList a full list of vertex surfaces in
         * (quad or quad-oct) coordinates for the given triangulation.
         */
        template <class Variant>
        void buildStandardFromReduced(NTriangulation* owner,
            const std::vector<NNormalSurfaceVector*>& reducedList);

        /**
         * Implements the one-template-argument version of
         * buildStandardFromReduced() using the specified bitmask type
         * to store zero sets.  See the one-template-argument
         * buildStandardFromReduced() for further information on this routine,
         * including important preconditions.
         *
         * The one-template-argument buildStandardFromReduced() simply
         * chooses an appropriate bitmask type and then calls this
         * routine, which does the real work.
         *
         * \pre The template argument \a BitmaskType can support
         * bitmasks of size 7 \a n (if we are using normal surfaces) or size
         * 10 \a n (if we are using almost normal surfaces), where \a n is
         * the number of tetrahedra in the given triangulation.
         */
        template <class Variant, class BitmaskType>
        void buildStandardFromReducedUsing(NTriangulation* owner,
            const std::vector<NNormalSurfaceVector*>& reducedList);

        /**
         * Converts a set of embedded vertex surfaces in (quad or quad-oct)
         * space to a set of embedded vertex surfaces in (standard normal or
         * standard almost normal) space.
         *
         * This is a generic implementation that performs the real work
         * for both quadToStandard() and quadOctToStandardAN().  See each
         * of those routines for further details as well as relevant
         * preconditions and postconditions.
         *
         * \pre The template argument \a Variant is either NormalSpec
         * or AlmostNormalSpec, according to whether we are implementing
         * quadToStandard() or quadOctToStandardAN() accordingly.
         */
        template <class Variant>
        NNormalSurfaceList* internalReducedToStandard() const;

        /**
         * Converts a set of embedded vertex surfaces in
         * (standard normal or standard almost normal) space to a set of
         * embedded vertex surfaces in (quad or quad-oct) space.
         *
         * This is a generic implementation that performs the real work
         * for both standardToQuad() and standardANToQuadOct().  See each
         * of those routines for further details as well as relevant
         * preconditions and postconditions.
         *
         * \pre The template argument \a Variant is either NormalSpec
         * or AlmostNormalSpec, according to whether we are implementing
         * standardToQuad() or standardANToQuadOct() accordingly.
         */
        template <class Variant>
        NNormalSurfaceList* internalStandardToReduced() const;

        /**
         * A thread class that actually performs the vertex normal surface
         * enumeration.
         */
        class VertexEnumerator : public NThread {
            private:
                NNormalSurfaceList* list;
                    /**< The normal surface list to be filled. */
                NTriangulation* triang;
                    /**< The triangulation upon which this normal
                         surface list will be based. */
                NProgressManager* manager;
                    /**< The progress manager through which progress is
                         reported, or 0 if no progress manager is in use. */

            public:
                /**
                 * Creates a new enumerator thread with the given
                 * parameters.
                 *
                 * @param newList the normal surface list to be filled.
                 * @param useTriang the triangulation upon which this
                 * normal surface list will be based.
                 * @param useManager the progress manager to use for
                 * progress reporting, or 0 if progress reporting is not
                 * required.
                 */
                VertexEnumerator(NNormalSurfaceList* newList,
                    NTriangulation* useTriang, NProgressManager* useManager);

                void* run(void*);
        };

        /**
         * A thread class that performs fundamental normal surface enumeration
         * using the primal Hilbert basis algorithm.
         */
        class FundPrimalEnumerator : public NThread {
            private:
                NNormalSurfaceList* list;
                    /**< The normal surface list to be filled. */
                NTriangulation* triang;
                    /**< The triangulation upon which this normal
                         surface list will be based. */
                NNormalSurfaceList* vtxSurfaces;
                    /**< The list of all vertex normal surfaces in \a triang,
                         or 0 if this information is not known. */
                NProgressManager* manager;
                    /**< The progress manager through which progress is
                         reported, or 0 if no progress manager is in use. */

            public:
                /**
                 * Creates a new enumerator thread with the given
                 * parameters.
                 *
                 * @param newList the normal surface list to be filled.
                 * @param useTriang the triangulation upon which this
                 * normal surface list will be based.
                 * @param useVtxSurafces the vertex normal surfaces in
                 * \a useTriang, or 0 if this information is not known.
                 * @param useManager the progress manager to use for
                 * progress reporting, or 0 if progress reporting is not
                 * required.
                 */
                FundPrimalEnumerator(NNormalSurfaceList* newList,
                    NTriangulation* useTriang,
                    NNormalSurfaceList* useVtxSurfaces,
                    NProgressManager* useManager);

                void* run(void*);
        };

        /**
         * A thread class that performs fundamental normal surface enumeration
         * using the dual Hilbert basis algorithm.
         */
        class FundDualEnumerator : public NThread {
            private:
                NNormalSurfaceList* list;
                    /**< The normal surface list to be filled. */
                NTriangulation* triang;
                    /**< The triangulation upon which this normal
                         surface list will be based. */
                NProgressManager* manager;
                    /**< The progress manager through which progress is
                         reported, or 0 if no progress manager is in use. */

            public:
                /**
                 * Creates a new enumerator thread with the given
                 * parameters.
                 *
                 * @param newList the normal surface list to be filled.
                 * @param useTriang the triangulation upon which this
                 * normal surface list will be based.
                 * @param useManager the progress manager to use for
                 * progress reporting, or 0 if progress reporting is not
                 * required.
                 */
                FundDualEnumerator(NNormalSurfaceList* newList,
                    NTriangulation* useTriang, NProgressManager* useManager);

                void* run(void*);
        };

    friend class regina::NXMLNormalSurfaceListReader;
};

/**
 * Returns a new normal surface vector of the appropriate length for the
 * given triangulation and the given flavour of coordinate system.
 * All elements of this vector will be initialised to zero.
 *
 * The new vector will be of the subclass of NNormalSurfaceVector
 * corresponding to the given flavour of coordinate system.  The caller
 * of this routine is responsible for destroying the new vector.
 *
 * \ifacespython Not present.
 *
 * @param triangulation the triangulation upon which the underlying
 * coordinate system is based.
 * @param flavour the flavour of coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalSurfaceList.
 * @return a new zero vector of the correct class and length.
 */
REGINA_API NNormalSurfaceVector* makeZeroVector(
        const NTriangulation* triangulation, int flavour);
/**
 * Creates a new set of normal surface matching equations for the
 * given triangulation using the given flavour of coordinate system.
 * The returned matrix will be newly allocated and its destruction will
 * be the responsibility of the caller of this routine.
 *
 * Each equation will be represented as a row of the matrix.
 * Each column of the matrix represents a coordinate in the given
 * flavour of coordinate system.
 *
 * @param triangulation the triangulation upon which these matching equations
 * will be based.
 * @param flavour the flavour of coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalSurfaceList.
 * @return a newly allocated set of matching equations.
 */
REGINA_API NMatrixInt* makeMatchingEquations(NTriangulation* triangulation,
    int flavour);
/**
 * Creates a new set of validity constraints representing the condition that
 * normal surfaces be embedded.  The validity constraints will be expressed
 * relative to the given flavour of coordinate system.
 *
 * \ifacespython Not present.
 *
 * @param triangulation the triangulation upon which these validity constraints
 * will be based.
 * @param flavour the flavour of coordinate system to be used;
 * this must be one of the predefined coordinate system
 * constants in NNormalSurfaceList.
 * @return a newly allocated set of constraints.
 */
REGINA_API NEnumConstraintList* makeEmbeddedConstraints(
    NTriangulation* triangulation, int flavour);

/*@}*/

// Inline functions for NNormalSurfaceList

inline NNormalSurfaceList::NNormalSurfaceList() {
}

inline NNormalSurfaceList::~NNormalSurfaceList() {
    for_each(surfaces.begin(), surfaces.end(), FuncDelete<NNormalSurface>());
}

inline int NNormalSurfaceList::getFlavour() const {
    return flavour;
}

inline bool NNormalSurfaceList::isEmbeddedOnly() const {
    return embedded;
}

inline unsigned long NNormalSurfaceList::getNumberOfSurfaces() const {
    return surfaces.size();
}

inline const NNormalSurface* NNormalSurfaceList::getSurface(
        unsigned long index) const {
    return surfaces[index];
}

inline ShareableObject* NNormalSurfaceList::getShareableObject() {
    return this;
}

inline bool NNormalSurfaceList::dependsOnParent() const {
    return true;
}

inline NMatrixInt* NNormalSurfaceList::recreateMatchingEquations() const {
    return makeMatchingEquations(getTriangulation(), flavour);
}

inline NNormalSurfaceList::VectorIterator::VectorIterator() {
}

inline NNormalSurfaceList::VectorIterator::VectorIterator(
        const NNormalSurfaceList::VectorIterator& cloneMe) :
        it_(cloneMe.it_) {
}

inline NNormalSurfaceList::VectorIterator& NNormalSurfaceList::VectorIterator::
        operator =(const NNormalSurfaceList::VectorIterator& cloneMe) {
    it_ = cloneMe.it_;
    return *this;
}

inline bool NNormalSurfaceList::VectorIterator::operator ==(
        const NNormalSurfaceList::VectorIterator& other) const {
    return (it_ == other.it_);
}

inline bool NNormalSurfaceList::VectorIterator::operator !=(
        const NNormalSurfaceList::VectorIterator& other) const {
    return (it_ != other.it_);
}

inline const NNormalSurfaceVector* NNormalSurfaceList::VectorIterator::
        operator *() const {
    return (*it_)->rawVector();
}

inline NNormalSurfaceList::VectorIterator& NNormalSurfaceList::VectorIterator::
        operator ++() {
    ++it_;
    return *this;
}

inline NNormalSurfaceList::VectorIterator NNormalSurfaceList::VectorIterator::
        operator ++(int) {
    return NNormalSurfaceList::VectorIterator(it_++);
}

inline NNormalSurfaceList::VectorIterator& NNormalSurfaceList::VectorIterator::
        operator --() {
    --it_;
    return *this;
}

inline NNormalSurfaceList::VectorIterator NNormalSurfaceList::VectorIterator::
        operator --(int) {
    return NNormalSurfaceList::VectorIterator(it_--);
}

inline NNormalSurfaceList::VectorIterator::VectorIterator(
        const std::vector<NNormalSurface*>::const_iterator& i) : it_(i) {
}

inline NNormalSurfaceList::VectorIterator NNormalSurfaceList::beginVectors()
        const {
    return VectorIterator(surfaces.begin());
}

inline NNormalSurfaceList::VectorIterator NNormalSurfaceList::endVectors()
        const {
    return VectorIterator(surfaces.end());
}

inline NNormalSurfaceList::SurfaceInserter::SurfaceInserter() : list(0),
        owner(0) {
}

inline NNormalSurfaceList::SurfaceInserter::SurfaceInserter(
        NNormalSurfaceList& newList, NTriangulation* newOwner) :
        list(&newList), owner(newOwner) {
}

inline NNormalSurfaceList::SurfaceInserter::SurfaceInserter(
        const SurfaceInserter& cloneMe) : list(cloneMe.list),
        owner(cloneMe.owner) {
}


inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator =(
        const SurfaceInserter& cloneMe) {
    list = cloneMe.list;
    owner = cloneMe.owner;
    return *this;
}

inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator =(
        NNormalSurface* surface) {
    list->surfaces.push_back(surface);
    return *this;
}

inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator =(
        NNormalSurfaceVector* vector) {
    list->surfaces.push_back(new NNormalSurface(owner, vector));
    return *this;
}

inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator *() {
    return *this;
}

inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator ++() {
    return *this;
}

inline NNormalSurfaceList::SurfaceInserter&
        NNormalSurfaceList::SurfaceInserter::operator ++(int) {
    return *this;
}

inline NNormalSurfaceList::NNormalSurfaceList(int newFlavour,
        bool embeddedOnly) : flavour(newFlavour), embedded(embeddedOnly) {
}

inline NNormalSurfaceList::VertexEnumerator::VertexEnumerator(
        NNormalSurfaceList* newList,
        NTriangulation* useTriang, NProgressManager* useManager) :
        list(newList), triang(useTriang), manager(useManager) {
}

inline NNormalSurfaceList::FundPrimalEnumerator::FundPrimalEnumerator(
        NNormalSurfaceList* newList, NTriangulation* useTriang,
        NNormalSurfaceList* newVtxSurfaces, NProgressManager* useManager) :
        list(newList), triang(useTriang), vtxSurfaces(newVtxSurfaces),
        manager(useManager) {
}

inline NNormalSurfaceList::FundDualEnumerator::FundDualEnumerator(
        NNormalSurfaceList* newList,
        NTriangulation* useTriang, NProgressManager* useManager) :
        list(newList), triang(useTriang), manager(useManager) {
}

} // namespace regina

#endif

