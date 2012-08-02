
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

/*! \file triangulation/ntriangulation.h
 *  \brief Deals with triangulations.
 */

#ifndef __NTRIANGULATION_H
#ifndef __DOXYGEN
#define __NTRIANGULATION_H
#endif

#include <map>
#include <memory>
#include <vector>

#include "regina-core.h"
#include "algebra/nabeliangroup.h"
#include "algebra/ngrouppresentation.h"
#include "file/nfilepropertyreader.h"
#include "packet/npacket.h"
#include "utilities/hashutils.h"
#include "utilities/nbooleans.h"
#include "utilities/nmarkedvector.h"
#include "utilities/nproperty.h"

// The following headers are necessary so that std::auto_ptr can invoke
// destructors where necessary.
#include "triangulation/nisomorphism.h"

// NOTE: More #includes follow after the class declarations.

namespace regina {

class NBoundaryComponent;
class NComponent;
class NEdge;
class NFace;
class NTetrahedron;
class NVertex;
class NGroupPresentation;
class NIsomorphism;
class NXMLPacketReader;
class NXMLTriangulationReader;

/**
 * \addtogroup triangulation Triangulations
 * Triangulations of 3-manifolds.
 * @{
 */

/**
 * Stores the triangulation of a 3-manifold along with its
 * various cellular structures and other information.
 *
 * When the triangulation is deleted, the corresponding
 * tetrahedra, the cellular structure and all other properties
 * will be deallocated.
 *
 * Faces, edges, vertices and components are always temporary;
 * whenever a change
 * occurs with the triangulation, these will be deleted and a new
 * skeletal structure will be calculated.  The same is true of various
 * other triangulation properties.
 *
 * The management of tetrahedra within a triangulation has become simpler and
 * safer as of Regina 4.90.  In older versions (Regina 4.6 and earlier),
 * users were required to create tetrahedra, individually add them to
 * triangulations, and manually notify a triangulation whenever its tetrahedron
 * gluings changed.  As of Regina 4.90, new tetrahedra are created using
 * NTriangulation::newTetrahedron() which automatically places them within
 * a triangulation, and all gluing changes are likewise communicated to the
 * triangulation automatically.  These are part of a larger suite of changes
 * (all designed to help the user avoid inconsistent states and accidental
 * crashes); see the NTetrahedron class notes for further details.
 *
 * \testpart
 *
 * \todo \feature Is the boundary incompressible?
 * \todo \featurelong Am I obviously a handlebody?  (Simplify and see
 * if there is nothing left).  Am I obviously not a handlebody?
 * (Compare homology with boundary homology).
 * \todo \featurelong Is the triangulation Haken?
 * \todo \featurelong What is the Heegaard genus?
 * \todo \featurelong Have a subcomplex as a child packet of a
 * triangulation.  Include routines to crush a subcomplex or to expand a
 * subcomplex to a normal surface.
 * \todo \featurelong Implement writeTextLong() for skeletal objects.
 */
class REGINA_API NTriangulation : public NPacket, public NFilePropertyReader {
    public:
        static const int packetType;

        typedef std::vector<NTetrahedron*>::const_iterator TetrahedronIterator;
            /**< Used to iterate through tetrahedra. */
        typedef std::vector<NFace*>::const_iterator FaceIterator;
            /**< Used to iterate through faces. */
        typedef std::vector<NEdge*>::const_iterator EdgeIterator;
            /**< Used to iterate through edges. */
        typedef std::vector<NVertex*>::const_iterator VertexIterator;
            /**< Used to iterate through vertices. */
        typedef std::vector<NComponent*>::const_iterator ComponentIterator;
            /**< Used to iterate through components. */
        typedef std::vector<NBoundaryComponent*>::const_iterator
                BoundaryComponentIterator;
            /**< Used to iterate through boundary components. */

        typedef std::map<std::pair<unsigned long, unsigned long>, double>
                TuraevViroSet;
            /**< A map from (r, whichRoot) pairs to Turaev-Viro invariants. */

    private:
        mutable bool calculatedSkeleton;
            /**< Has the skeleton been calculated? */

        NMarkedVector<NTetrahedron> tetrahedra;
            /**< The tetrahedra that form the triangulation. */
        mutable NMarkedVector<NFace> faces;
            /**< The faces in the triangulation skeleton. */
        mutable NMarkedVector<NEdge> edges;
            /**< The edges in the triangulation skeleton. */
        mutable NMarkedVector<NVertex> vertices;
            /**< The vertices in the triangulation skeleton. */
        mutable NMarkedVector<NComponent> components;
            /**< The components that form the triangulation. */
        mutable NMarkedVector<NBoundaryComponent> boundaryComponents;
            /**< The components that form the boundary of the
                 triangulation. */

        mutable bool valid;
            /**< Is the triangulation valid? */
        mutable bool ideal;
            /**< Is the triangulation ideal? */
        mutable bool standard;
            /**< Is the triangulation standard? */
        mutable bool orientable;
            /**< Is the triangulation orientable? */

        mutable NProperty<NGroupPresentation, StoreManagedPtr> fundamentalGroup;
            /**< Fundamental group of the triangulation. */
        mutable NProperty<NAbelianGroup, StoreManagedPtr> H1;
            /**< First homology group of the triangulation. */
        mutable NProperty<NAbelianGroup, StoreManagedPtr> H1Rel;
            /**< Relative first homology group of the triangulation
             *   with respect to the boundary. */
        mutable NProperty<NAbelianGroup, StoreManagedPtr> H1Bdry;
            /**< First homology group of the boundary. */
        mutable NProperty<NAbelianGroup, StoreManagedPtr> H2;
            /**< Second homology group of the triangulation. */

        mutable NProperty<bool> twoSphereBoundaryComponents;
            /**< Does the triangulation contain any 2-sphere boundary
                 components? */
        mutable NProperty<bool> negativeIdealBoundaryComponents;
            /**< Does the triangulation contain any boundary components
                 that are ideal and have negative Euler characteristic? */

        mutable NProperty<bool> zeroEfficient;
            /**< Is the triangulation zero-efficient? */
        mutable NProperty<bool> splittingSurface;
            /**< Does the triangulation have a normal splitting surface? */

        mutable NProperty<bool> threeSphere;
            /**< Is this a triangulation of a 3-sphere? */
        mutable NProperty<bool> threeBall;
            /**< Is this a triangulation of a 3-dimensional ball? */
        mutable NProperty<bool> solidTorus;
            /**< Is this a triangulation of the solid torus? */

        mutable TuraevViroSet turaevViroCache;
            /**< The set of Turaev-Viro invariants that have already
                 been calculated. */

    public:
        /**
         * \name Constructors and Destructors
         */
        /*@{*/

        /**
         * Default constructor.
         *
         * Creates an empty triangulation.
         */
        NTriangulation();
        /**
         * Copy constructor.
         *
         * Creates a new triangulation identical to the given triangulation.
         * The packet tree structure and packet label are \e not copied.
         *
         * @param cloneMe the triangulation to clone.
         */
        NTriangulation(const NTriangulation& cloneMe);
        /**
         * Destroys this triangulation.
         *
         * The constituent tetrahedra, the cellular structure and all other
         * properties will also be deallocated.
         */
        virtual ~NTriangulation();

        /*@}*/
        /**
         * (end: Constructors and Destructors)
         */

        /**
         * \name Packet Administration
         */
        /*@{*/

        virtual int getPacketType() const;
        virtual std::string getPacketTypeName() const;

        virtual void writePacket(NFile& out) const;
        virtual void writeTextShort(std::ostream& out) const;
        virtual void writeTextLong(std::ostream& out) const;
        virtual bool dependsOnParent() const;
        virtual void readIndividualProperty(NFile& infile, unsigned propType);

        /*@}*/
        /**
         * (end: Packet Administration)
         */

        /**
         * \name Tetrahedra
         */
        /*@{*/

        /**
         * Returns the number of tetrahedra in the triangulation.
         *
         * @return the number of tetrahedra.
         */
        unsigned long getNumberOfTetrahedra() const;
        /**
         * A dimension-agnostic alias for getNumberOfTetrahedra().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See getNumberOfTetrahedra() for further information.
         */
        unsigned long getNumberOfSimplices() const;
        /**
         * Returns all tetrahedra in the triangulation.
         *
         * The reference returned will remain valid
         * for as long as the triangulation exists,
         * always reflecting the tetrahedra currently in the
         * triangulation.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all tetrahedra.
         */
        const std::vector<NTetrahedron*>& getTetrahedra() const;
        /**
         * A dimension-agnostic alias for getTetrahedra().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See getTetrahedra() for further information.
         */
        const std::vector<NTetrahedron*>& getSimplices() const;
        /**
         * Returns the tetrahedron with the given index number in the
         * triangulation.
         * Note that tetrahedron indexing may change when a tetrahedron
         * is added or removed from the triangulation.
         *
         * @param index specifies which tetrahedron to return; this
         * value should be between 0 and getNumberOfTetrahedra()-1
         * inclusive.
         * @return the <tt>index</tt>th tetrahedron in the
         * triangulation.
         */
        NTetrahedron* getTetrahedron(unsigned long index);
        /**
         * A dimension-agnostic alias for getTetrahedron().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See getTetrahedron() for further information.
         */
        NTetrahedron* getSimplex(unsigned long index);
        /**
         * Returns the tetrahedron with the given index number in the
         * triangulation.
         * Note that tetrahedron indexing may change when a tetrahedron
         * is added or removed from the triangulation.
         *
         * @param index specifies which tetrahedron to return; this
         * value should be between 0 and getNumberOfTetrahedra()-1
         * inclusive.
         * @return the <tt>index</tt>th tetrahedron in the
         * triangulation.
         */
        const NTetrahedron* getTetrahedron(unsigned long index) const;
        /**
         * A dimension-agnostic alias for getTetrahedron().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See getTetrahedron() for further information.
         */
        const NTetrahedron* getSimplex(unsigned long index) const;
        /**
         * Returns the index of the given tetrahedron in the
         * triangulation.
         *
         * Note that tetrahedron indexing may change when a tetrahedron
         * is added or removed from the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getTetrahedronIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * tetrahedron belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given tetrahedron is contained in this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.  If you are passing the result of some other
         * routine that \e might return null (such as
         * NTetrahedron::adjacentTetrahedron), it might be worth explicitly
         * testing for null beforehand.
         *
         * @param tet specifies which tetrahedron to find in the
         * triangulation.
         * @return the index of the specified tetrahedron, where 0 is
         * the first tetrahedron, 1 is the second and so on.
         */
        long tetrahedronIndex(const NTetrahedron* tet) const;
        /**
         * A dimension-agnostic alias for tetrahedronIndex().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See tetrahedronIndex() for further information.
         */
        long simplexIndex(const NTetrahedron* tet) const;
        /**
         * Returns the index of the given tetrahedron in the
         * triangulation.
         *
         * Note that tetrahedron indexing may change when a tetrahedron
         * is added or removed from the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative tetrahedronIndex() is
         * fast, but insists that \a tet belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param tet specifies which tetrahedron to find in the
         * triangulation.
         * @return the index of the specified tetrahedron, where 0 is
         * the first tetrahedron, 1 is the second and so on.  If the
         * tetrahedron is not contained in the triangulation, a negative
         * number is returned.
         */
        long getTetrahedronIndex(const NTetrahedron* tet) const;
        /**
         * Creates a new tetrahedron and adds it to this triangulation.
         * The new tetrahedron will have an empty description.
         * All four faces of the new tetrahedron will be boundary faces.
         *
         * @return the new tetrahedron.
         */
        NTetrahedron* newTetrahedron();
        /**
         * A dimension-agnostic alias for newTetrahedron().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See newTetrahedron() for further information.
         */
        NTetrahedron* newSimplex();
        /**
         * Creates a new tetrahedron with the given description and adds
         * it to this triangulation.
         * All four faces of the new tetrahedron will be boundary faces.
         *
         * @param desc the description to assign to the new tetrahedron.
         * @return the new tetrahedron.
         */
        NTetrahedron* newTetrahedron(const std::string& desc);
        /**
         * A dimension-agnostic alias for newTetrahedron().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See newTetrahedron() for further information.
         */
        NTetrahedron* newSimplex(const std::string& desc);
        /**
         * Inserts the given tetrahedron into the triangulation.
         * No face gluings anywhere will be examined or altered.
         *
         * The new tetrahedron will be assigned a higher index in the
         * triangulation than all tetrahedra already present.
         *
         * \pre The given tetrahedron does not already belong to a
         * different triangulation (though already belonging to \e this
         * triangulation is perfectly fine).
         *
         * \deprecated Users should create tetrahedra by calling
         * newTetrahedron() or newTetrahedron(const std::string&), which
         * will add the tetrahedron to the triangulation automatically.
         *
         * \warning As of Regina 4.90, this routine will also add any
         * neighbouring tetrahedra that do not yet belong to a
         * triangulation; moreover, this addition is recursive.  This is done
         * to ensure that, whenever one tetrahedron belongs to a
         * triangulation, everything that it is joined to (directly or
         * indirectly) also belongs to that same triangulation.
         * See the NTetrahedron class notes for further details on how
         * tetrahedron management has changed in Regina 4.90 and above.
         *
         * \ifacespython Since this triangulation takes ownership
         * of the given tetrahedron, the python object containing the
         * given tetrahedron becomes a null object and should no longer
         * be used.
         *
         * @param tet the tetrahedron to insert.
         */
        void addTetrahedron(NTetrahedron* tet);
        /**
         * Removes the given tetrahedron from the triangulation.
         * All faces glued to this tetrahedron will be unglued.
         * The tetrahedron will be deallocated.
         *
         * \pre The given tetrahedron exists in the triangulation.
         *
         * \warning This routine has changed behaviour as of Regina 4.90.
         * In older versions of Regina, the tetrahedron was returned to
         * the user.  As of Regina 4.90, the tetrahedron is now destroyed
         * immediately.
         *
         * @param tet the tetrahedron to remove.
         */
        void removeTetrahedron(NTetrahedron* tet);
        /**
         * A dimension-agnostic alias for removeTetrahedron().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See removeTetrahedron() for further information.
         */
        void removeSimplex(NTetrahedron* tet);
        /**
         * Removes the tetrahedron with the given index number
         * from the triangulation.  Note that tetrahedron indexing may
         * change when a tetrahedron is added or removed from the
         * triangulation.
         *
         * All faces glued to this tetrahedron will be unglued.
         * The tetrahedron will be deallocated.
         *
         * \warning This routine has changed behaviour as of Regina 4.90.
         * In older versions of Regina, the tetrahedron was returned to
         * the user.  As of Regina 4.90, the tetrahedron is now destroyed
         * immediately.
         *
         * @param index specifies which tetrahedron to remove; this
         * should be between 0 and getNumberOfTetrahedra()-1 inclusive.
         */
        void removeTetrahedronAt(unsigned long index);
        /**
         * A dimension-agnostic alias for removeTetrahedronAt().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See removeTetrahedronAt() for further information.
         */
        void removeSimplexAt(unsigned long index);
        /**
         * Removes all tetrahedra from the triangulation.
         * All tetrahedra will be deallocated.
         */
        void removeAllTetrahedra();
        /**
         * A dimension-agnostic alias for removeAllTetrahedra().
         * This is to assist with writing dimension-agnostic code that
         * can be reused to work in different dimensions.
         * 
         * Here "simplex" refers to a top-dimensional simplex (which for
         * 3-manifold triangulations means a tetrahedron).
         * 
         * See removeAllTetrahedra() for further information.
         */
        void removeAllSimplices();
        /**
         * Swaps the contents of this and the given triangulation.
         * That is, all tetrahedra that belong to this triangulation
         * will be moved to \a other, and all tetrahedra that belong to
         * \a other will be moved to this triangulation.
         *
         * All NTetrahedron pointers or references will remain valid.
         *
         * @param other the triangulation whose contents should be
         * swapped with this.
         */
        void swapContents(NTriangulation& other);
        /**
         * Moves the contents of this triangulation into the given
         * destination triangulation, without destroying any pre-existing
         * contents.  That is, all tetrahedra that currently belong to
         * \a dest will remain there, and all tetrahedra that belong to this
         * triangulation will be moved across as additional
         * tetrahedra in \a dest.
         *
         * All NTetrahedron pointers or references will remain valid.
         * After this operation, this triangulation will be empty.
         *
         * @param dest the triangulation to which tetrahedra should be
         * moved.
         */
        void moveContentsTo(NTriangulation& dest);
        /**
         * This routine now does nothing, and should not be used.
         *
         * \deprecated In Regina versions 4.6 and earlier, this routine
         * was used to manually notify the triangulation that the gluings
         * of tetrahedra had changed.  In Regina 4.90 and later this
         * notification is automatic.  This routine now does nothing at
         * all, and can safely be removed from any existing code.
         */
        void gluingsHaveChanged();

        /*@}*/
        /**
         * (end: Tetrahedra)
         */

        /**
         * \name Skeletal Queries
         */
        /*@{*/

        /**
         * Returns the number of boundary components in this
         * triangulation.  Note that each ideal vertex forms its own
         * boundary component.
         *
         * @return the number of boundary components.
         */
        unsigned long getNumberOfBoundaryComponents() const;
        /**
         * Returns the number of components in this triangulation.
         *
         * @return the number of components.
         */
        unsigned long getNumberOfComponents() const;
        /**
         * Returns the number of vertices in this triangulation.
         *
         * @return the number of vertices.
         */
        unsigned long getNumberOfVertices() const;
        /**
         * Returns the number of edges in this triangulation.
         *
         * @return the number of edges.
         */
        unsigned long getNumberOfEdges() const;
        /**
         * Returns the number of faces in this triangulation.
         *
         * @return the number of faces.
         */
        unsigned long getNumberOfFaces() const;

        /**
         * Returns all components of this triangulation.
         *
         * Bear in mind that each time the triangulation changes, the
         * components will be deleted and replaced with new
         * ones.  Thus the objects contained in this list should be
         * considered temporary only.
         *
         * This reference to the list however will remain valid and
         * up-to-date for as long as the triangulation exists.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all components.
         */
        const std::vector<NComponent*>& getComponents() const;
        /**
         * Returns all boundary components of this triangulation.
         * Note that each ideal vertex forms its own boundary component.
         *
         * Bear in mind that each time the triangulation changes, the
         * boundary components will be deleted and replaced with new
         * ones.  Thus the objects contained in this list should be
         * considered temporary only.
         *
         * This reference to the list however will remain valid and
         * up-to-date for as long as the triangulation exists.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all boundary components.
         */
        const std::vector<NBoundaryComponent*>& getBoundaryComponents() const;
        /**
         * Returns all vertices of this triangulation.
         *
         * Bear in mind that each time the triangulation changes, the
         * vertices will be deleted and replaced with new
         * ones.  Thus the objects contained in this list should be
         * considered temporary only.
         *
         * This reference to the list however will remain valid and
         * up-to-date for as long as the triangulation exists.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all vertices.
         */
        const std::vector<NVertex*>& getVertices() const;
        /**
         * Returns all edges of this triangulation.
         *
         * Bear in mind that each time the triangulation changes, the
         * edges will be deleted and replaced with new
         * ones.  Thus the objects contained in this list should be
         * considered temporary only.
         *
         * This reference to the list however will remain valid and
         * up-to-date for as long as the triangulation exists.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all edges.
         */
        const std::vector<NEdge*>& getEdges() const;
        /**
         * Returns all faces of this triangulation.
         *
         * Bear in mind that each time the triangulation changes, the
         * faces will be deleted and replaced with new
         * ones.  Thus the objects contained in this list should be
         * considered temporary only.
         *
         * This reference to the list however will remain valid and
         * up-to-date for as long as the triangulation exists.
         *
         * \ifacespython This routine returns a python list.
         *
         * @return the list of all faces.
         */
        const std::vector<NFace*>& getFaces() const;
        /**
         * Returns the requested triangulation component.
         *
         * Bear in mind that each time the triangulation changes, the
         * components will be deleted and replaced with new
         * ones.  Thus this object should be considered temporary only.
         *
         * @param index the index of the desired component, ranging from 0
         * to getNumberOfComponents()-1 inclusive.
         * @return the requested component.
         */
        NComponent* getComponent(unsigned long index) const;
        /**
         * Returns the requested triangulation boundary component.
         *
         * Bear in mind that each time the triangulation changes, the
         * boundary components will be deleted and replaced with new
         * ones.  Thus this object should be considered temporary only.
         *
         * @param index the index of the desired boundary
         * component, ranging from 0
         * to getNumberOfBoundaryComponents()-1 inclusive.
         * @return the requested boundary component.
         */
        NBoundaryComponent* getBoundaryComponent(unsigned long index) const;
        /**
         * Returns the requested triangulation vertex.
         *
         * Bear in mind that each time the triangulation changes, the
         * vertices will be deleted and replaced with new
         * ones.  Thus this object should be considered temporary only.
         *
         * @param index the index of the desired vertex, ranging from 0
         * to getNumberOfVertices()-1 inclusive.
         * @return the requested vertex.
         */
        NVertex* getVertex(unsigned long index) const;
        /**
         * Returns the requested triangulation edge.
         *
         * Bear in mind that each time the triangulation changes, the
         * edges will be deleted and replaced with new
         * ones.  Thus this object should be considered temporary only.
         *
         * @param index the index of the desired edge, ranging from 0
         * to getNumberOfEdges()-1 inclusive.
         * @return the requested edge.
         */
        NEdge* getEdge(unsigned long index) const;
        /**
         * Returns the requested triangulation face.
         *
         * Bear in mind that each time the triangulation changes, the
         * faces will be deleted and replaced with new
         * ones.  Thus this object should be considered temporary only.
         *
         * @param index the index of the desired face, ranging from 0
         * to getNumberOfFaces()-1 inclusive.
         * @return the requested face.
         */
        NFace* getFace(unsigned long index) const;
        /**
         * Returns the index of the given component in the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getComponentIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * component belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given component belongs to this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.
         *
         * @param component specifies which component to find in the
         * triangulation.
         * @return the index of the specified component, where 0 is the first
         * component, 1 is the second and so on.
         */
        long componentIndex(const NComponent* component) const;
        /**
         * Returns the index of the given component in the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative componentIndex() is
         * fast, but insists that \a component belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param component specifies which component to find in the
         * triangulation.
         * @return the index of the specified component, where 0 is the first
         * component, 1 is the second and so on.  If the given component
         * is not part of this triangulation, a negative number is returned.
         */
        long getComponentIndex(const NComponent* component) const;
        /**
         * Returns the index of the given boundary component
         * in the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getBoundaryComponentIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * boundary component belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given boundary component belongs to this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.
         *
         * @param bc specifies which boundary component to find in the
         * triangulation.
         * @return the index of the specified boundary component,
         * where 0 is the first boundary component, 1 is the second and so on. 
         */
        long boundaryComponentIndex(const NBoundaryComponent* bc) const;
        /**
         * Returns the index of the given boundary component
         * in the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative boundaryComponentIndex() is
         * fast, but insists that \a bc belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param bc specifies which boundary component to find in the
         * triangulation.
         * @return the index of the specified boundary component,
         * where 0 is the first boundary component,
         * 1 is the second and so on.  If the given boundary component
         * is not part of this triangulation, a negative number is returned.
         */
        long getBoundaryComponentIndex(const NBoundaryComponent* bc) const;
        /**
         * Returns the index of the given vertex in the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getVertexIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * vertex belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given vertex belongs to this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.
         *
         * @param vertex specifies which vertex to find in the
         * triangulation.
         * @return the index of the specified vertex, where 0 is the first
         * vertex, 1 is the second and so on.
         */
        long vertexIndex(const NVertex* vertex) const;
        /**
         * Returns the index of the given vertex in the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative vertexIndex() is
         * fast, but insists that \a vertex belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param vertex specifies which vertex to find in the
         * triangulation.
         * @return the index of the specified vertex, where 0 is the first
         * vertex, 1 is the second and so on.  If the given vertex
         * is not part of this triangulation, a negative number is returned.
         */
        long getVertexIndex(const NVertex* vertex) const;
        /**
         * Returns the index of the given edge in the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getEdgeIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * edge belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given edge belongs to this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.
         *
         * @param edge specifies which edge to find in the
         * triangulation.
         * @return the index of the specified edge, where 0 is the first
         * edge, 1 is the second and so on.
         */
        long edgeIndex(const NEdge* edge) const;
        /**
         * Returns the index of the given edge in the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative edgeIndex() is
         * fast, but insists that \a edge belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param edge specifies which edge to find in the
         * triangulation.
         * @return the index of the specified edge, where 0 is the first
         * edge, 1 is the second and so on.  If the given edge
         * is not part of this triangulation, a negative number is returned.
         */
        long getEdgeIndex(const NEdge* edge) const;
        /**
         * Returns the index of the given face in the triangulation.
         *
         * This routine was introduced in Regina 4.5, and replaces the
         * old getFaceIndex().  The name has been changed
         * because, unlike the old routine, it requires that the given
         * face belongs to the triangulation (a consequence of
         * some significant memory optimisations).
         *
         * \pre The given face belongs to this triangulation.
         *
         * \warning Passing a null pointer to this routine will probably
         * crash your program.
         *
         * @param face specifies which face to find in the
         * triangulation.
         * @return the index of the specified face, where 0 is the first
         * face, 1 is the second and so on.
         */
        long faceIndex(const NFace* face) const;
        /**
         * Returns the index of the given face in the triangulation.
         *
         * \deprecated As of Regina 4.5, this is now a slow linear
         * time routine (a consequence of some significant memory
         * optimisations).  The alternative faceIndex() is
         * fast, but insists that \a face belongs to the triangulation.
         * This slower routine will be removed at some point in the future.
         *
         * @param face specifies which face to find in the
         * triangulation.
         * @return the index of the specified face, where 0 is the first
         * face, 1 is the second and so on.  If the given face
         * is not part of this triangulation, a negative number is returned.
         */
        long getFaceIndex(const NFace* face) const;

        /**
         * Determines if this triangulation contains any two-sphere
         * boundary components.
         *
         * @return \c true if and only if there is at least one
         * two-sphere boundary component.
         */
        bool hasTwoSphereBoundaryComponents() const;
        /**
         * Determines if this triangulation contains any ideal boundary
         * components with negative Euler characteristic.
         *
         * @return \c true if and only if there is at least one such
         * boundary component.
         */
        bool hasNegativeIdealBoundaryComponents() const;

        /*@}*/
        /**
         * (end: Skeletal Queries)
         */

        /**
         * \name Isomorphism Testing
         */
        /*@{*/

        /**
         * Determines if this triangulation is combinatorially
         * isomorphic to the given triangulation.
         *
         * Specifically, this routine determines if there is a
         * one-to-one and onto boundary complete combinatorial
         * isomorphism from this triangulation to \a other.  Boundary
         * complete isomorphisms are described in detail in the
         * NIsomorphism class notes.
         *
         * In particular, note that this triangulation and \a other must
         * contain the same number of tetrahedra for such an isomorphism
         * to exist.
         *
         * \todo \opt Improve the complexity by choosing a tetrahedron
         * mapping from each component and following gluings to
         * determine the others.
         *
         * If a boundary complete isomorphism is found, the details of
         * this isomorphism are returned.  The isomorphism is newly
         * constructed, and so to assist with memory management is
         * returned as a std::auto_ptr.  Thus, to test whether an
         * isomorphism exists without having to explicitly deal with the
         * isomorphism itself, you can call
         * <tt>if (isIsomorphicTo(other).get())</tt> and the newly
         * created isomorphism (if it exists) will be automatically
         * destroyed.
         *
         * @param other the triangulation to compare with this one.
         * @return details of the isomorphism if the two triangulations
         * are combinatorially isomorphic, or a null pointer otherwise.
         */
        std::auto_ptr<NIsomorphism> isIsomorphicTo(const NTriangulation& other)
            const;

        /**
         * Determines if an isomorphic copy of this triangulation is
         * contained within the given triangulation, possibly as a
         * subcomplex of some larger component (or components).
         *
         * Specifically, this routine determines if there is a boundary
         * incomplete combinatorial isomorphism from this triangulation
         * to \a other.  Boundary incomplete isomorphisms are described
         * in detail in the NIsomorphism class notes.
         *
         * In particular, note that boundary faces of this triangulation
         * need not correspond to boundary faces of \a other, and that
         * \a other can contain more tetrahedra than this triangulation.
         *
         * If a boundary incomplete isomorphism is found, the details of
         * this isomorphism are returned.  The isomorphism is newly
         * constructed, and so to assist with memory management is
         * returned as a std::auto_ptr.  Thus, to test whether an
         * isomorphism exists without having to explicitly deal with the
         * isomorphism itself, you can call
         * <tt>if (isContainedIn(other).get())</tt> and the newly
         * created isomorphism (if it exists) will be automatically
         * destroyed.
         *
         * If more than one such isomorphism exists, only one will be
         * returned.  For a routine that returns all such isomorphisms,
         * see findAllSubcomplexesIn().
         *
         * @param other the triangulation in which to search for an
         * isomorphic copy of this triangulation.
         * @return details of the isomorphism if such a copy is found,
         * or a null pointer otherwise.
         */
        std::auto_ptr<NIsomorphism> isContainedIn(const NTriangulation& other)
            const;

        /**
         * Finds all ways in which an isomorphic copy of this triangulation
         * is contained within the given triangulation, possibly as a
         * subcomplex of some larger component (or components).
         *
         * This routine behaves identically to isContainedIn(), except
         * that instead of returning just one isomorphism (which may be
         * boundary incomplete and need not be onto), all such isomorphisms
         * are returned.
         *
         * See the isContainedIn() notes for additional information.
         *
         * The isomorphisms that are found will be inserted into the
         * given list.  These isomorphisms will be newly created, and
         * the caller of this routine is responsible for destroying
         * them.  The given list will not be emptied before the new
         * isomorphisms are inserted.
         *
         * \ifacespython Not present.
         *
         * @param other the triangulation in which to search for
         * isomorphic copies of this triangulation.
         * @param results the list in which any isomorphisms found will
         * be stored.
         * @return the number of isomorphisms that were found.
         */
        unsigned long findAllSubcomplexesIn(const NTriangulation& other,
                std::list<NIsomorphism*>& results) const;

        /*@}*/
        /**
         * (end: Isomorphism Testing)
         */

        /**
         * \name Basic Properties
         */
        /*@{*/

        /**
         * Returns the Euler characteristic of this triangulation.
         * This will be evaluated strictly as \a V-E+F-T.
         *
         * Note that this routine handles cusps in a non-standard way.
         * Since it computes the Euler characteristic of the
         * triangulation (and not the underlying manifold), this routine
         * will treat each cusp as a single vertex, and \e not as
         * a surface boundary component.
         *
         * For a routine that handles cusps properly (i.e., treats them
         * as surface boundary components when computing the Euler
         * characteristic), see getEulerCharManifold() instead.
         *
         * This routine was previously called getEulerCharacteristic() in
         * Regina 4.3.1 and earlier.  It was renamed in Regina 4.4 to
         * clarify the non-standard handling of cusps.
         *
         * @return the Euler characteristic of this triangulation.
         */
        long getEulerCharTri() const;

        /**
         * Returns the Euler characteristic of the corresponding compact
         * 3-manifold.
         *
         * Instead of simply calculating \a V-E+F-T, this routine also:
         *
         * - treats ideal vertices as surface boundary components
         *   (i.e., effectively truncates them);
         * - truncates invalid boundary vertices (i.e., boundary vertices
         *   whose links are not discs);
         * - truncates the projective plane cusps at the midpoints of invalid
         *   edges (edges identified with themselves in reverse).
         *
         * For ideal triangulations, this routine therefore computes
         * the proper Euler characteristic of the manifold (unlike
         * getEulerCharTri(), which does not).
         *
         * For triangulations whose vertex links are all spheres or discs,
         * this routine and getEulerCharTri() give identical results.
         *
         * @return the Euler characteristic of the corresponding compact
         * manifold.
         */
        long getEulerCharManifold() const;

        /**
         * A deprecated alias for getEulerCharTri().
         *
         * This routine calculates the Euler characteristic of this
         * triangulation.  Since it treats cusps in a non-standard way,
         * it was renamed to getEulerCharTri() in Regina 4.4 to clarify
         * that this might differ from the Euler characteristic of the
         * corresponding compact manifold.
         *
         * See getEulerCharTri() for further details.
         *
         * \deprecated This routine will be removed in a future version
         * of Regina.  Please use getEulerCharTri() instead.
         *
         * @return the Euler characteristic of this triangulation.
         */
        long getEulerCharacteristic() const;

        /**
         * Determines if this triangulation is valid.
         * A triangulation is valid unless there is some vertex whose
         * link has boundary but is not a disc (i.e., a vertex for which
         * NVertex::getLink() returns NVertex::NON_STANDARD_BDRY),
         * or unless there is some edge glued to itself in reverse
         * (i.e., an edge for which NEdge::isValid() returns \c false).
         *
         * @return \c true if and only if this triangulation is valid.
         */
        bool isValid() const;
        /**
         * Determines if this triangulation is ideal.
         * This is the case if and only if one of the vertex links
         * is closed and not a 2-sphere.
         * Note that the triangulation is not required to be valid.
         *
         * @return \c true if and only if this triangulation is ideal.
         */
        bool isIdeal() const;
        /**
         * Determines if this triangulation is standard.
         * This is the case if and only if every vertex is standard.
         * See NVertex::isStandard() for further details.
         *
         * @return \c true if and only if this triangulation is
         * standard.
         */
        bool isStandard() const;
        /**
         * Determines if this triangulation has any boundary faces.
         *
         * @return \c true if and only if there are boundary faces.
         */
        bool hasBoundaryFaces() const;
        /**
         * Determines if this triangulation is closed.
         * This is the case if and only if it has no boundary.
         * Note that ideal triangulations are not closed.
         *
         * @return \c true if and only if this triangulation is closed.
         */
        bool isClosed() const;
        /**
         * Determines if this triangulation is orientable.
         *
         * @return \c true if and only if this triangulation is
         * orientable.
         */
        bool isOrientable() const;

        /**
         * Determines if this triangulation is oriented; that is, if
         * tetrahedron vertices are labelled in a way that preserves
         * orientation across adjacent faces.  Specifically, this
         * routine returns \c true if and only if every gluing permutation
         * has negative sign.
         *
         * Note that \e orientable triangulations are not always \e oriented
         * by default.  You can call orient() if you need the tetrahedra
         * to be oriented consistently as described above.
         *
         * A non-orientable triangulation can never be oriented.
         *
         * @return \c true if and only if all tetrahedra are oriented
         * consistently.
         *
         * @author Matthias Goerner
         */
        bool isOriented() const;

        /**
         * Determines if this triangulation is ordered; that is, if
         * tetrahedron vertices are labelled so that all gluing
         * permutations are order-preserving on the faces.
         * Equivalently, this tests whether the edges of the triangulation
         * can all be oriented such that they induce a consistent ordering
         * on the vertices of each tetrahedron.
         *
         * Triangulations are not ordered by default, and indeed some
         * cannot be ordered at all.  The routine order() will attempt
         * to relabel tetrahedron vertices to give an ordered triangulation.
         *
         * @return \c true if and only if all gluing permutations are
         * order preserving on the faces.
         *
         * @author Matthias Goerner
         */
        bool isOrdered() const;

        /**
         * Determines if this triangulation is connected.
         *
         * @return \c true if and only if this triangulation is
         * connected.
         */
        bool isConnected() const;

        /*@}*/
        /**
         * (end: Basic Properties)
         */

        /**
         * \name Algebraic Properties
         */
        /*@{*/

        /**
         * Returns the fundamental group of this triangulation.
         * If this triangulation contains any ideal or non-standard
         * vertices, the fundamental group will be
         * calculated as if each such vertex had been truncated.
         *
         * If this triangulation contains any invalid edges, the
         * calculations will be performed <b>without</b> any truncation
         * of the corresponding projective plane cusp.  Thus if a
         * barycentric subdivision is performed on the triangulation, the
         * result of getFundamentalGroup() will change.
         *
         * Bear in mind that each time the triangulation changes, the
         * fundamental group will be deleted.  Thus the reference that is
         * returned from this routine should not be kept for later use.
         * Instead, getFundamentalGroup() should be called again; this will
         * be instantaneous if the group has already been calculated.
         *
         * Note that this triangulation is not required to be valid
         * (see isValid()).
         *
         * \pre This triangulation has at most one component.
         *
         * @return the fundamental group.
         */
        const NGroupPresentation& getFundamentalGroup() const;
        /**
         * Notifies the triangulation that you have simplified the
         * presentation of its fundamental group.  The old group
         * presentation will be destroyed, and this triangulation will
         * take ownership of the new (hopefully simpler) group that is
         * passed.
         *
         * This routine is useful for situations in which some external
         * body (such as GAP) has simplified the group presentation
         * better than Regina can.
         *
         * Regina does <i>not</i> verify that the new group presentation
         * is equivalent to the old, since this is - well, hard.
         *
         * If the fundamental group has not yet been calculated for this
         * triangulation, this routine will nevertheless take ownership
         * of the new group, under the assumption that you have worked
         * out the group through some other clever means without ever
         * having needed to call getFundamentalGroup() at all.
         *
         * Note that this routine will not fire a packet change event.
         *
         * @param newGroup a new (and hopefully simpler) presentation
         * of the fundamental group of this triangulation.
         */
        void simplifiedFundamentalGroup(NGroupPresentation* newGroup);
        /**
         * Returns the first homology group for this triangulation.
         * If this triangulation contains any ideal or non-standard
         * vertices, the homology group will be
         * calculated as if each such vertex had been truncated.
         *
         * If this triangulation contains any invalid edges, the
         * calculations will be performed <b>without</b> any truncation
         * of the corresponding projective plane cusp.  Thus if a
         * barycentric subdivision is performed on the triangulation, the
         * result of getHomologyH1() will change.
         *
         * Bear in mind that each time the triangulation changes, the
         * homology groups will be deleted.  Thus the reference that is
         * returned from this routine should not be kept for later use.
         * Instead, getHomologyH1() should be called again; this will be
         * instantaneous if the group has already been calculated.
         *
         * Note that this triangulation is not required to be valid
         * (see isValid()).
         *
         * @return the first homology group.
         */
        const NAbelianGroup& getHomologyH1() const;
        /**
         * Returns the relative first homology group with
         * respect to the boundary for this triangulation.
         * Note that ideal vertices are considered part of the boundary.
         *
         * Bear in mind that each time the triangulation changes, the
         * homology groups will be deleted.  Thus the reference that is
         * returned from this routine should not be kept for later use.
         * Instead, getHomologyH1Rel() should be called again; this will be
         * instantaneous if the group has already been calculated.
         *
         * \pre This triangulation is valid.
         *
         * @return the relative first homology group with respect to the
         * boundary.
         */
        const NAbelianGroup& getHomologyH1Rel() const;
        /**
         * Returns the first homology group of the
         * boundary for this triangulation.
         * Note that ideal vertices are considered part of the boundary.
         *
         * Bear in mind that each time the triangulation changes, the
         * homology groups will be deleted.  Thus the reference that is
         * returned from this routine should not be kept for later use.
         * Instead, getHomologyH1Bdry() should be called again; this will be
         * instantaneous if the group has already been calculated.
         *
         * This routine is fairly fast, since it deduces the homology of
         * each boundary component through knowing what kind of surface
         * it is.
         *
         * \pre This triangulation is valid.
         *
         * @return the first homology group of the boundary.
         */
        const NAbelianGroup& getHomologyH1Bdry() const;
        /**
         * Returns the second homology group for this triangulation.
         * If this triangulation contains any ideal vertices,
         * the homology group will be
         * calculated as if each such vertex had been truncated.
         * The algorithm used calculates various first homology groups
         * and uses homology and cohomology theorems to deduce the
         * second homology group.
         *
         * Bear in mind that each time the triangulation changes, the
         * homology groups will be deleted.  Thus the reference that is
         * returned from this routine should not be kept for later use.
         * Instead, getHomologyH2() should be called again; this will be
         * instantaneous if the group has already been calculated.
         *
         * \pre This triangulation is valid.
         *
         * @return the second homology group.
         */
        const NAbelianGroup& getHomologyH2() const;
        /**
         * Returns the second homology group with coefficients in Z_2
         * for this triangulation.
         * If this triangulation contains any ideal vertices,
         * the homology group will be
         * calculated as if each such vertex had been truncated.
         * The algorithm used calculates the relative first homology group
         * with respect to the boundary and uses homology and cohomology
         * theorems to deduce the second homology group.
         *
         * This group will simply be the direct sum of several copies of
         * Z_2, so the number of Z_2 terms is returned.
         *
         * \pre This triangulation is valid.
         *
         * @return the number of Z_2 terms in the second homology group
         * with coefficients in Z_2.
         */
        unsigned long getHomologyH2Z2() const;
        /**
         * Computes the Turaev-Viro state sum invariant of this
         * 3-manifold based upon the given initial data.
         *
         * The initial data is as described in the paper of Turaev and
         * Viro, "State sum invariants of 3-manifolds and quantum
         * 6j-symbols", Topology, vol. 31, no. 4, 1992, pp 865-902.
         *
         * In particular, Section 7 describes the initial data as
         * determined by an integer r >=3 and a root of unity q0 of
         * degree 2r for which q0^2 is a primitive root of unity of
         * degree r.
         *
         * These invariants, although computed in the complex field,
         * should all be reals.  Thus the return type is an ordinary
         * double.
         *
         * \pre This triangulation is valid, closed and non-empty.
         *
         * @param r the integer r as described above; this must be at
         * least 3.
         * @param whichRoot determines q0 to be the root of unity
         * e^(2i * Pi * whichRoot / 2r); this argument must be strictly
         * between 0 and 2r and must have no common factors with r.
         * @return the requested Turaev-Viro invariant.
         *
         * @see allCalculatedTuraevViro
         */
        double turaevViro(unsigned long r, unsigned long whichRoot) const;
        /**
         * Returns the set of all Turaev-Viro state sum invariants that
         * have already been calculated for this 3-manifold.
         *
         * Turaev-Viro invariants are described by an (r, whichRoot)
         * pair as described in the turaevViro() notes.  The set
         * returned by this routine maps (r, whichRoot) pairs to the
         * corresponding invariant values.
         *
         * Each time turaevViro() is called, the result will be stored
         * in this set (as well as being returned to the user).  This
         * set will be emptied whenever the triangulation is modified.
         *
         * \ifacespython Not present.
         *
         * @return the set of all Turaev-Viro invariants that have
         * already been calculated.
         *
         * @see turaevViro
         */
        const TuraevViroSet& allCalculatedTuraevViro() const;

        /*@}*/
        /**
         * (end: Algebraic Properties)
         */

        /**
         * \name Normal Surface Properties
         */
        /*@{*/

        /**
         * Determines if this triangulation is 0-efficient.
         * A triangulation is 0-efficient if its only normal spheres and
         * discs are vertex linking, and if it has no 2-sphere boundary
         * components.
         *
         * @return \c true if and only if this triangulation is
         * 0-efficient.
         */
        bool isZeroEfficient();
        /**
         * Is it already known whether or not this triangulation is
         * 0-efficient?
         * See isZeroEfficient() for further details.
         *
         * If this property is already known, future calls to
         * isZeroEfficient() will be very fast (simply returning the
         * precalculated value).
         *
         * \warning This routine does not actually tell you \e whether
         * this triangulation is 0-efficient; it merely tells you whether
         * the answer has already been computed.
         *
         * @return \c true if and only if this property is already known.
         */
        bool knowsZeroEfficient() const;
        /**
         * Determines whether this triangulation has a normal splitting
         * surface.  See NNormalSurface::isSplitting() for details
         * regarding normal splitting surfaces.
         *
         * \pre This triangulation is connected.  If the triangulation
         * is not connected, this routine will still return a result but
         * that result will be unreliable.
         *
         * @return \c true if and only if this triangulation has a
         * normal splitting surface.
         */
        bool hasSplittingSurface();
        /**
         * Is it already known whether or not this triangulation has a
         * splitting surface?
         * See hasSplittingSurface() for further details.
         *
         * If this property is already known, future calls to
         * hasSplittingSurface() will be very fast (simply returning the
         * precalculated value).
         *
         * \warning This routine does not actually tell you \e whether
         * this triangulation has a splitting surface; it merely tells you
         * whether the answer has already been computed.
         *
         * @return \c true if and only if this property is already known.
         */
        bool knowsSplittingSurface() const;

        /*@}*/
        /**
         * (end: Normal Surface Properties)
         */

        /**
         * \name Skeletal Transformations
         */
        /*@{*/

        /**
         * Produces a maximal forest in the 1-skeleton of the
         * triangulation boundary.
         * Both given sets will be emptied and the edges and vertices of
         * the maximal forest will be placed into them.
         * A vertex that forms its own boundary component (such as an
         * ideal vertex) will still be placed in \c vertexSet.
         *
         * Note that the edge and vertex pointers returned will become
         * invalid once the triangulation has changed.
         *
         * \deprecated In Regina 5.0, this routine will be changed to
         * use std::set instead of stdhash::hash_set.  This is part of
         * the broader plan to remove any need for the non-standard
         * hash_set and hash_map.
         *
         * \ifacespython Not present.
         *
         * @param edgeSet the set to be emptied and into which the edges
         * of the maximal forest will be placed.
         * @param vertexSet the set to be emptied and into which the
         * vertices of the maximal forest will be placed.
         */
        void maximalForestInBoundary(
                stdhash::hash_set<NEdge*, HashPointer>& edgeSet,
                stdhash::hash_set<NVertex*, HashPointer>& vertexSet) const;
        /**
         * Produces a maximal forest in the triangulation's 1-skeleton.
         * The given set will be emptied and will have the edges of the
         * maximal forest placed into it.  It can be specified whether
         * or not different boundary components may be joined by the
         * maximal forest.
         *
         * An edge leading to an ideal vertex is still a
         * candidate for inclusion in the maximal forest.  For the
         * purposes of this algorithm, any ideal vertex will be treated
         * as any other vertex (and will still be considered part of its
         * own boundary component).
         *
         * Note that the edge pointers returned will become
         * invalid once the triangulation has changed.
         *
         * \deprecated In Regina 5.0, this routine will be changed to
         * use std::set instead of stdhash::hash_set.  This is part of
         * the broader plan to remove any need for the non-standard
         * hash_set and hash_map.
         *
         * \ifacespython Not present.
         *
         * @param edgeSet the set to be emptied and into which the edges
         * of the maximal forest will be placed.
         * @param canJoinBoundaries \c true if and only if different
         * boundary components are allowed to be joined by the maximal
         * forest.
         */
        void maximalForestInSkeleton(
                stdhash::hash_set<NEdge*, HashPointer>& edgeSet,
                bool canJoinBoundaries = true) const;
        /**
         * Produces a maximal forest in the triangulation's dual
         * 1-skeleton.  The given set will be emptied and will have the
         * faces corresponding to the edges of the maximal forest in the
         * dual 1-skeleton placed into it.
         *
         * Note that the face pointers returned will become invalid once
         * the triangulation has changed.
         *
         * \deprecated In Regina 5.0, this routine will be changed to
         * use std::set instead of stdhash::hash_set.  This is part of
         * the broader plan to remove any need for the non-standard
         * hash_set and hash_map.
         *
         * \ifacespython Not present.
         *
         * @param faceSet the set to be emptied and into which the faces
         * representing the maximal forest will be placed.
         */
        void maximalForestInDualSkeleton(
                stdhash::hash_set<NFace*, HashPointer>& faceSet) const;

        /**
         * Attempts to simplify the triangulation as intelligently as
         * possible without further input.  This routine will attempt to
         * reduce both the number of tetrahedra and the number of boundary
         * faces (with the number of tetrahedra as its priority).
         *
         * Currently this routine uses simplifyToLocalMinimum() in
         * combination with random 4-4 moves, book opening moves and
         * book closing moves.
         *
         * \warning The specific behaviour of this routine may well
         * change between releases.
         *
         * \todo \opt Include random 2-3 moves to get out of wells.
         *
         * @return \c true if and only if the triangulation was changed.
         */
        bool intelligentSimplify();
        /**
         * Uses all known simplification moves to reduce the triangulation
         * monotonically to some local minimum number of tetrahedra.
         * Note that this will probably not give a globally minimal
         * triangulation; see intelligentSimplify() for further
         * assistance in achieving this goal.
         *
         * The moves used include 3-2, 2-0 (edge and vertex),
         * 2-1 and boundary shelling moves.
         *
         * Note that moves that do not reduce the number of tetrahedra
         * (such as 4-4 moves or book opening moves) are not used in this
         * routine.  Such moves do however feature in intelligentSimplify().
         *
         * \warning The specific behaviour of this routine is
         * very likely to change between releases.
         *
         * @param perform \c true if we are to perform the
         * simplifications, or \c false if we are only to investigate
         * whether simplifications are possible (defaults to \c true).
         * @return if \a perform is \c true, this routine returns
         * \c true if and only if the triangulation was changed to
         * reduce the number of tetrahedra; if \a perform is \c false,
         * this routine returns \c true if and only if it determines
         * that it is capable of performing such a change.
         */
        bool simplifyToLocalMinimum(bool perform = true);

        /**
         * Checks the eligibility of and/or performs a 3-2 move
         * about the given edge.
         * This involves replacing the three tetrahedra joined at that
         * edge with two tetrahedra joined by a face.
         * This can be done iff (i) the edge is valid and non-boundary,
         * and (ii) the three tetrahedra are distinct.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a e)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool threeTwoMove(NEdge* e, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a 2-3 move
         * about the given face.
         * This involves replacing the two tetrahedra joined at that
         * face with three tetrahedra joined by an edge.
         * This can be done iff the two tetrahedra are distinct.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a f)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given face is a face of this triangulation.
         *
         * @param f the face about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool twoThreeMove(NFace* f, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a 4-4 move
         * about the given edge.
         * This involves replacing the four tetrahedra joined at that
         * edge with four tetrahedra joined along a different edge.
         * Consider the octahedron made up of the four original
         * tetrahedra; this has three internal axes.  The initial four
         * tetrahedra meet along the given edge which forms one of these
         * axes; the new tetrahedra will meet along a different axis.
         * This move can be done iff (i) the edge is valid and non-boundary,
         * and (ii) the four tetrahedra are distinct.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a e)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge about which to perform the move.
         * @param newAxis Specifies which axis of the octahedron the new
         * tetrahedra should meet along; this should be 0 or 1.
         * Consider the four original tetrahedra in the order described
         * by NEdge::getEmbeddings(); call these tetrahedra 0, 1, 2 and
         * 3.  If \a newAxis is 0, the new axis will separate tetrahedra
         * 0 and 1 from 2 and 3.  If \a newAxis is 1, the new axis will
         * separate tetrahedra 1 and 2 from 3 and 0.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool fourFourMove(NEdge* e, int newAxis, bool check = true,
                bool perform = true);
        /**
         * Checks the eligibility of and/or performs a 2-0 move
         * about the given edge of degree 2.
         * This involves taking the two tetrahedra joined at that edge
         * and squashing them flat.  This can be done if:
         *
         * - the edge is valid and non-boundary;
         *
         * - the two tetrahedra are distinct;
         *
         * - the edges opposite \c e in each tetrahedron are distinct and
         *   not both boundary;
         *
         * - if faces \a f1 and \a f2 of one tetrahedron are to be flattened
         *   onto faces \a g1 and \a g2 of the other respectively, then
         *   (a) \a f1 and \a g1 are distinct,
         *   (b) \a f2 and \a g2 are distinct,
         *   (c) we do not have both \a f1 = \a g2 and \a g1 = \a f2,
         *   (d) we do not have both \a f1 = \a f2 and \a g1 = \a g2, and
         *   (e) we do not have two of the faces boundary and the other
         *   two identified.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a e)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool twoZeroMove(NEdge* e, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a 2-0 move
         * about the given vertex of degree 2.
         * This involves taking the two tetrahedra joined at that vertex
         * and squashing them flat.
         * This can be done if:
         *
         * - the vertex is non-boundary and has a 2-sphere vertex link;
         *
         * - the two tetrahedra are distinct;
         *
         * - the faces opposite \c v in each tetrahedron are distinct and not
         *   both boundary;
         *
         * - the two tetrahedra meet each other on all three faces touching
         *   the vertex (as opposed to meeting each other on one face and
         *   being glued to themselves along the other two).
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a v)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given vertex is a vertex of this triangulation.
         *
         * @param v the vertex about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool twoZeroMove(NVertex* v, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a 2-1 move
         * about the given edge.
         * This involves taking an edge meeting only one tetrahedron
         * just once and merging that tetrahedron with one of the
         * tetrahedra joining it.
         *
         * This can be done assuming the following conditions:
         *
         * - The edge must be valid and non-boundary.
         *
         * - The two remaining faces of the tetrahedron are not joined, and
         *   the tetrahedron face opposite the given endpoint of the edge is
         *   not boundary.
         *
         * - Consider the second tetrahedron to be merged (the one
         *   joined along the face opposite the given endpoint of the edge).
         *   Moreover, consider the two edges of this second tetrahedron
         *   that run from the (identical) vertices of the original
         *   tetrahedron not touching \c e to the vertex of the second
         *   tetrahedron not touching the original tetrahedron.  These edges
         *   must be distinct and may not both be in the boundary.
         *
         * There are additional "distinct and not both boundary" conditions
         * on faces of the second tetrahedron, but those follow automatically
         * from the final condition above.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a e)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge about which to perform the move.
         * @param edgeEnd the end of the edge \e opposite that at
         * which the second tetrahedron (to be merged) is joined.
         * The end is 0 or 1, corresponding to the labelling (0,1) of
         * the vertices of the edge as described in
         * NEdgeEmbedding::getVertices().
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool twoOneMove(NEdge* e, int edgeEnd,
                bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a book opening move
         * about the given face.
         * This involves taking a face meeting the boundary along two
         * edges, and ungluing it to create two new boundary faces
         * (thus exposing the tetrahedra it initially joined).
         * This move is the inverse of the closeBook() move, and is
         * used to open the way for new shellBoundary() moves.
         *
         * This move can be done if:
         *
         * - the face meets the boundary in precisely two edges (and thus
         *   also joins two tetrahedra);
         *
         * - the vertex between these two edges is a standard boundary
         *   vertex (its link is a disc);
         *
         * - the remaining edge of the face (which is internal to the
         *   triangulation) is valid.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a f)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given face is a face of this triangulation.
         *
         * @param f the face about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool openBook(NFace* f, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a book closing move
         * about the given boundary edge.
         * This involves taking a boundary edge of the triangulation and
         * folding together the two boundary faces on either side.  This
         * move is the inverse of the openBook() move, and is used to
         * simplify the boundary of the triangulation.
         * This move can be done if:
         *
         * - the edge \a e is a boundary edge;
         *
         * - the two boundary faces that it joins are distinct;
         *
         * - the two vertices opposite \a e in each of these boundary faces
         *   are valid and distinct;
         *
         * - if edges \a e1 and \a e2 of one boundary face are to be
         *   folded onto edges \a f1 and \a f2 of the other boundary
         *   face respectively, then we do not have both \a e1 = \a e2
         *   and \a f1 = \a f2.
         *
         * There are in fact several other "distinctness" conditions on
         * the edges \a e1, \a e2, \a f1 and \a f2, but they follow
         * automatically from the "distinct vertices" condition above.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a f)
         * can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge about which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool closeBook(NEdge* e, bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a boundary shelling
         * move on the given tetrahedron.
         * This involves simply popping off a tetrahedron that touches
         * the boundary.
         * This can be done if:
         *
         * - all edges of the tetrahedron are valid;
         *
         * - precisely one, two or three faces of the tetrahedron lie in
         *   the boundary;
         *
         * - if one face lies in the boundary, then the opposite vertex
         *   does not lie in the boundary, and no two of the remaining three
         *   edges are identified;
         *
         * - if two faces lie in the boundary, then the remaining edge does
         *   not lie in the boundary, and the remaining two faces of the
         *   tetrahedron are not identified.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects can no longer be used.
         *
         * \pre If the move is being performed and no
         * check is being run, it must be known in advance that the move
         * is legal.
         * \pre The given tetrahedron is a tetrahedron of this triangulation.
         *
         * @param t the tetrahedron upon which to perform the move.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the requested move may be performed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool shellBoundary(NTetrahedron* t,
                bool check = true, bool perform = true);
        /**
         * Checks the eligibility of and/or performs a collapse of
         * an edge in such a way that the topology of the manifold
         * does not change and the number of vertices of the triangulation
         * decreases by one.
         *
         * If the routine is asked to both check and perform, the move
         * will only be performed if the check shows it is legal.
         *
         * Note that after performing this move, all skeletal objects
         * (faces, components, etc.) will be reconstructed, which means
         * any pointers to old skeletal objects (such as the argument \a e)
         * can no longer be used.
         *
         * The eligibility requirements for this move are somewhat
         * involved, and are discussed in detail in the collapseEdge()
         * source code for those who are interested.
         *
         * \pre If the move is being performed and no check is being run,
         * it must be known in advance that the move is legal.
         * \pre The given edge is an edge of this triangulation.
         *
         * @param e the edge to collapse.
         * @param check \c true if we are to check whether the move is
         * allowed (defaults to \c true).
         * @param perform \c true if we are to perform the move
         * (defaults to \c true).
         * @return If \a check is \c true, the function returns \c true
         * if and only if the given edge may be collapsed
         * without changing the topology of the manifold.  If \a check
         * is \c false, the function simply returns \c true.
         */
        bool collapseEdge(NEdge* e, bool check = true, bool perform = true);

        /**
         * Reorders the tetrahedra of this triangulation using a
         * breadth-first search, so that small-numbered tetrahedra are
         * adjacent to other small-numbered tetrahedra.
         *
         * Specifically, the reordering will operate as follows.
         * Tetrahedron 0 will remain tetrahedron 0.  Its immediate
         * neighbours will be numbered 1, 2, 3 and 4 (though if these
         * neighbours are not distinct then of course fewer labels will
         * be required).  Their immediate neighbours will in turn be
         * numbered 5, 6, and so on, ultimately following a breadth-first
         * search throughout the entire triangulation.
         *
         * If the optional argument \a reverse is \c true, then tetrahedron
         * numbers will be assigned in reverse order.  That is, tetrahedron 0
         * will become tetrahedron \a n-1, its neighbours will become
         * tetrahedra \a n-2 down to \a n-5, and so on.
         *
         * @param reverse \c true if the new tetrahedron numbers should
         * be assigned in reverse order, as described above.
         */
        void reorderTetrahedraBFS(bool reverse = false);

        /**
         * Relabels tetrahedron vertices in this triangulation so that
         * all tetrahedra are oriented consistently, if possible.
         *
         * This routine works by flipping vertices 2 and 3 of each
         * tetrahedron with negative orientation.  The result will be a
         * triangulation where the tetrahedron vertices are labelled in
         * a way that preserves orientation across adjacent faces.
         * In particular, every gluing permutation will have negative sign.
         *
         * If this triangulation includes both orientable and
         * non-orientable components, the orientable components will be
         * oriented as described above and the non-orientable
         * components will be left untouched.
         *
         * @author Matthias Goerner
         */
        void orient();

        /**
         * Relabels tetrahedron vertices in this triangulation to give
         * an ordered triangulation, if possible.
         *
         * To be an ordered triangulation, all face gluings (when restricted
         * to the face) must be order preserving. In other words, it must
         * be possible to orient all edges of the triangulation in such a
         * fashion that they are consistent with the ordering of the vertices
         * in each tetrahedron.
         *
         * If it is possible to order this triangulation, the vertices
         * of each tetrahedron will be relabelled accordingly and this
         * routine will return \c true.  Otherwise, this routine will
         * return \c false and the triangulation will not be changed.
         *
         * \warning This routine may be slow, since it backtracks
         * through all possible edge orientations until a consistent one
         * has been found.
         *
         * @param forceOriented \c true if the triangulation must be
         * both ordered and \e oriented, in which case this routine will
         * return \c false if the triangulation cannot be oriented and
         * ordered at the same time.  See orient() for further details.
         * @return \c true if the triangulation has been successfully ordered
         * as described above, or \c false if not.
         *
         * @author Matthias Goerner
         */
        bool order(bool forceOriented = false);

        /*@}*/
        /**
         * (end: Skeletal Transformations)
         */

        /**
         * \name Decompositions
         */
        /*@{*/

        /**
         * Splits a disconnected triangulation into many smaller
         * triangulations, one for each component.  The new component
         * triangulations will be inserted as children of the given
         * parent packet.  The original triangulation will be left
         * unchanged.
         *
         * If the given parent packet is 0, the new component
         * triangulations will be inserted as children of this
         * triangulation.
         *
         * This routine can optionally assign unique (and sensible)
         * packet labels to each of the new component triangulations.
         * Note however that uniqueness testing may be slow, so this
         * assignment of labels should be disabled if the component
         * triangulations are only temporary objects used as part
         * of a larger routine.
         *
         * @param componentParent the packet beneath which the new
         * component triangulations will be inserted, or 0 if they
         * should be inserted directly beneath this triangulation.
         * @param setLabels \c true if the new component triangulations
         * should be assigned unique packet labels, or \c false if
         * they should be left without labels at all.
         * @return the number of new component triangulations
         * constructed.
         */
        unsigned long splitIntoComponents(NPacket* componentParent = 0,
            bool setLabels = true);
        /**
         * Splits this triangulation into its connected sum
         * decomposition.  The individual prime 3-manifold triangulations
         * that make up this decomposition will be inserted as children
         * of the given parent packet.  The original triangulation will
         * be left unchanged.
         *
         * Note that this routine is currently only available for
         * closed orientable triangulations; see the list of
         * preconditions for full details.  The 0-efficiency prime
         * decomposition algorithm of Jaco and Rubinstein is used.
         *
         * If the given parent packet is 0, the new prime summand
         * triangulations will be inserted as children of this
         * triangulation.
         *
         * This routine can optionally assign unique (and sensible)
         * packet labels to each of the new prime summand triangulations.
         * Note however that uniqueness testing may be slow, so this
         * assignment of labels should be disabled if the summand
         * triangulations are only temporary objects used as part
         * of a larger routine.
         *
         * If this is a triangulation of a 3-sphere, no prime summand
         * triangulations will be created at all.
         *
         * \warning The algorithms used in this routine rely on normal
         * surface theory and so can be very slow for larger triangulations.
         * For 3-sphere testing, see the routine isThreeSphere() which
         * uses faster methods where possible.
         *
         * \pre This triangulation is valid, closed, orientable and
         * connected.
         *
         * @param primeParent the packet beneath which the new prime
         * summand triangulations will be inserted, or 0 if they
         * should be inserted directly beneath this triangulation.
         * @param setLabels \c true if the new prime summand triangulations
         * should be assigned unique packet labels, or \c false if
         * they should be left without labels at all.
         * @return the number of prime summands created, 0 if this
         * triangulation is a 3-sphere or 0 if this triangulation does
         * not meet the preconditions described above.
         */
        unsigned long connectedSumDecomposition(NPacket* primeParent = 0,
            bool setLabels = true);
        /**
         * Determines whether this is a triangulation of a 3-sphere.
         *
         * This routine relies upon a combination of Rubinstein's 3-sphere
         * recognition algorithm and Jaco and Rubinstein's 0-efficiency
         * prime decomposition algorithm.
         *
         * \warning The algorithms used in this routine rely on normal
         * surface theory and so can be very slow for larger
         * triangulations (although faster tests are used where possible).
         * The routine knowsThreeSphere() can be called to see if this
         * property is already known or if it happens to be very fast to
         * calculate for this triangulation.
         *
         * @return \c true if and only if this is a 3-sphere triangulation.
         */
        bool isThreeSphere() const;
        /**
         * Is it already known (or trivial to determine) whether or not this
         * is a triangulation of a 3-sphere?  See isThreeSphere() for
         * further details.
         *
         * If this property is indeed already known, future calls to
         * isThreeSphere() will be very fast (simply returning the
         * precalculated value).
         *
         * If this property is not already known, this routine will
         * nevertheless run some very fast preliminary tests to see if the
         * answer is obviously no.  If so, it will store \c false as the
         * precalculated value for isThreeSphere() and this routine will
         * return \c true.
         *
         * Otherwise a call to isThreeSphere() may potentially require more
         * significant work, and so this routine will return \c false.
         *
         * \warning This routine does not actually tell you \e whether
         * this triangulation forms a 3-sphere; it merely tells you whether
         * the answer has already been computed (or is very easily computed).
         *
         * @return \c true if and only if this property is already known
         * or trivial to calculate.
         */
        bool knowsThreeSphere() const;
        /**
         * Determines whether this is a triangulation of a 3-dimensional ball.
         *
         * This routine is based on isThreeSphere(), which in turn combines
         * Rubinstein's 3-sphere recognition algorithm with Jaco and
         * Rubinstein's 0-efficiency prime decomposition algorithm.
         *
         * \warning The algorithms used in this routine rely on normal
         * surface theory and so can be very slow for larger
         * triangulations (although faster tests are used where possible).
         * The routine knowsBall() can be called to see if this
         * property is already known or if it happens to be very fast to
         * calculate for this triangulation.
         *
         * @return \c true if and only if this is a triangulation of a
         * 3-dimensional ball.
         */
        bool isBall() const;
        /**
         * Is it already known (or trivial to determine) whether or not this
         * is a triangulation of a 3-dimensional ball?  See isBall() for
         * further details.
         *
         * If this property is indeed already known, future calls to isBall()
         * will be very fast (simply returning the precalculated value).
         *
         * If this property is not already known, this routine will
         * nevertheless run some very fast preliminary tests to see if the
         * answer is obviously no.  If so, it will store \c false as the
         * precalculated value for isBall() and this routine will
         * return \c true.
         *
         * Otherwise a call to isBall() may potentially require more
         * significant work, and so this routine will return \c false.
         *
         * \warning This routine does not actually tell you \e whether
         * this triangulation forms a ball; it merely tells you whether
         * the answer has already been computed (or is very easily computed).
         *
         * @return \c true if and only if this property is already known
         * or trivial to calculate.
         */
        bool knowsBall() const;
        /**
         * Converts this into a 0-efficient triangulation of the same
         * underlying 3-manifold.  A triangulation is 0-efficient if its
         * only normal spheres and discs are vertex linking, and if it has
         * no 2-sphere boundary components.
         *
         * Note that this routine is currently only available for
         * closed orientable triangulations; see the list of
         * preconditions for details.  The 0-efficiency algorithm of
         * Jaco and Rubinstein is used.
         *
         * If the underlying 3-manifold is prime, it can always be made
         * 0-efficient (with the exception of the special cases RP3 and
         * S2xS1 as noted below).  In this case the original triangulation
         * will be modified directly and 0 will be returned.
         *
         * If the underyling 3-manifold is RP3 or S2xS1, it cannot
         * be made 0-efficient; in this case the original triangulation
         * will be reduced to a two-tetrahedron minimal triangulation
         * and 0 will again be returned.
         *
         * If the underlying 3-manifold is not prime, it cannot be made
         * 0-efficient.  In this case the original triangulation will
         * remain unchanged and a new connected sum decomposition will
         * be returned.  This will be presented as a newly allocated
         * container packet with one child triangulation for each prime
         * summand.
         *
         * \warning The algorithms used in this routine rely on normal
         * surface theory and so can be very slow for larger triangulations.
         *
         * \pre This triangulation is valid, closed, orientable and
         * connected.
         *
         * @return 0 if the underlying 3-manifold is prime (in which
         * case the original triangulation was modified directly), or
         * a newly allocated connected sum decomposition if the
         * underlying 3-manifold is composite (in which case the
         * original triangulation was not changed).
         */
        NPacket* makeZeroEfficient();
        /**
         * Determines whether this is a triangulation of the solid
         * torus; that is, the unknot complement.  This routine can be
         * used on a triangulation with real boundary faces, or on an
         * ideal triangulation (in which case all ideal vertices will
         * be assumed to be truncated).
         *
         * \warning The algorithms used in this routine rely on normal
         * surface theory and so might be very slow for larger
         * triangulations (although faster tests are used where possible).
         * The routine knowsSolidTorus() can be called to see if this
         * property is already known or if it happens to be very fast to
         * calculate for this triangulation.
         *
         * @return \c true if and only if this is either a real (compact)
         * or ideal (non-compact) triangulation of the solid torus.
         */
        bool isSolidTorus() const;
        /**
         * Is it already known (or trivial to determine) whether or not this
         * is a triangulation of a solid torus (that is, the unknot
         * complement)?  See isSolidTorus() for further details.
         *
         * If this property is indeed already known, future calls to
         * isSolidTorus() will be very fast (simply returning the
         * precalculated value).
         *
         * If this property is not already known, this routine will
         * nevertheless run some very fast preliminary tests to see if the
         * answer is obviously no.  If so, it will store \c false as the
         * precalculated value for isSolidTorus() and this routine will
         * return \c true.
         *
         * Otherwise a call to isSolidTorus() may potentially require more
         * significant work, and so this routine will return \c false.
         *
         * \warning This routine does not actually tell you \e whether
         * this triangulation forms a solid torus; it merely tells you whether
         * the answer has already been computed (or is very easily computed).
         *
         * @return \c true if and only if this property is already known
         * or trivial to calculate.
         */
        bool knowsSolidTorus() const;

        /**
         * Searches for a compressing disc within the underlying
         * 3-manifold.
         *
         * Let \a M be the underlying 3-manifold and let \a B be its
         * boundary.  By a <i>compressing disc</i>, we mean a disc \a D
         * properly embedded in \a M, where the boundary of \a D
         * lies in \a B but does not bound a disc in \a B.
         *
         * This routine will first call the heuristic routine
         * hasSimpleCompressingDisc() in the hope of obtaining a fast
         * answer.  If this fails, it will run a full enumeration of
         * vertex normal surfaces, which could be extremely slow.
         * For further details on the underlying algorithms, see
         * "The Weber-Seifert dodecahedral space is non-Haken",
         * Benjamin A. Burton, J. Hyam Rubinstein and Stephan Tillmann,
         * Trans. Amer. Math. Soc. 364:2 (2012), pp. 911-932; and
         * "Algorithms for the complete decomposition of a closed 3-manifold",
         * William Jaco and Jeffrey L. Tollefson, Illinois J. Math. 39 (1995),
         * 358-406.
         *
         * This routine will work on a copy of this triangulation, not
         * the original.  In particular, the copy will be simplified, which
         * means that there is no harm in calling this routine on an
         * unsimplified triangulation.
         *
         * If this triangulation has no boundary components, this
         * routine will simply return \c false.
         *
         * \pre This triangulation is valid and is not ideal.
         * \pre The underlying 3-manifold is irreducible.
         *
         * \warning This routine can be infeasibly slow for large
         * triangulations, since it may need to perform a full enumeration
         * of vertex normal surfaces, and since it might perform "large"
         * operations on these surfaces such as cutting along them.
         * See hasSimpleCompressingDisc() for a "heuristic shortcut"
         * that is faster but might not give a definitive answer.
         *
         * @return \c true if the underlying 3-manifold contains a
         * compressing disc, or \c false if it does not.
         */
        bool hasCompressingDisc() const;

        /**
         * Searches for a "simple" compressing disc inside this
         * triangulation.
         *
         * Let \a M be the underlying 3-manifold and let \a B be its
         * boundary.  By a <i>compressing disc</i>, we mean a disc \a D
         * properly embedded in \a M, where the boundary of \a D
         * lies in \a B but does not bound a disc in \a B.
         *
         * By a \a simple compressing disc, we mean a compressing disc
         * that has a very simple combinatorial structure (here "simple"
         * is subject to change; see the warning below).  Examples
         * include the compressing disc inside a 1-tetrahedron solid
         * torus, or a compressing disc formed from a single internal face
         * surrounded by three boundary edges.
         *
         * The purpose of this routine is to avoid enumerating normal
         * surfaces within a large triangulation where possible.  This
         * routine is relatively fast, and if it returns \c true then this
         * 3-manifold definitely contains a compressing disc.  If this
         * routine returns \c false then there might or might not be a
         * compressing disc; the user will need to perform a full normal
         * surface enumeration using hasCompressingDisc() to be sure.
         *
         * This routine will work on a copy of this triangulation, not
         * the original.  In particular, the copy will be simplified, which
         * means that there is no harm in calling this routine on an
         * unsimplified triangulation.
         *
         * If this triangulation has no boundary components, this
         * routine will simply return \c false.
         *
         * For further information on this test, see "The Weber-Seifert
         * dodecahedral space is non-Haken", Benjamin A. Burton,
         * J. Hyam Rubinstein and Stephan Tillmann,
         * Trans. Amer. Math. Soc. 364:2 (2012), pp. 911-932.
         *
         * \warning The definition of "simple" is subject to change in
         * future releases of Regina.  That is, this routine may be
         * expanded over time to identify more types of compressing discs
         * (thus making it more useful as a "heuristic shortcut").
         *
         * \pre This triangulation is valid and is not ideal.
         *
         * @return \c true if a simple compressing disc was found,
         * or \c false if not.  Note that even with a return value of
         * \c false, there might still be a compressing disc (just not
         * one with a simple combinatorial structure).
         */
        bool hasSimpleCompressingDisc() const;

        /*@}*/
        /**
         * (end: Decompositions)
         */

        /**
         * \name Subdivisions, Extensions and Covers
         */
        /*@{*/

        /**
         * Converts this triangulation into its double cover.
         * Each orientable component will be duplicated, and each
         * non-orientable component will be converted into its
         * orientable double cover.
         */
        void makeDoubleCover();

        /**
         * Converts an ideal triangulation into a finite triangulation.
         * All ideal or non-standard vertices are truncated and thus
         * converted into real boundary components made from unglued
         * faces of tetrahedra.
         *
         * Note that this operation is a loose converse of finiteToIdeal().
         *
         * \warning Currently, this routine subdivides all tetrahedra as
         * if <i>all</i> vertices (not just some) were ideal.
         * This may lead to more tetrahedra than are necessary.
         *
         * \warning Currently, the presence of an invalid edge will force
         * the triangulation to be subdivided regardless of the value of
         * parameter \a forceDivision.  The final triangulation will
         * still have the projective plane cusp caused by the invalid
         * edge.
         *
         * \todo \optlong Have this routine only use as many tetrahedra
         * as are necessary, leaving finite vertices alone.
         *
         * @param forceDivision specifies what to do if the triangulation
         * has no ideal or non-standard vertices.
         * If \c true, the triangulation will be
         * subdivided anyway, as if all vertices were ideal.  If
         * \c false (the default), the triangulation will be left alone.
         *
         * @return \c true if and only if the triangulation was changed.
         * @author David Letscher
         */
        bool idealToFinite(bool forceDivision = false);

        /**
         * Converts each real boundary component into a cusp (i.e., an ideal
         * vertex).  Only boundary components formed from real
         * tetrahedron faces will be affected; ideal boundary components
         * are already cusps and so will not be changed.
         *
         * One side-effect of this operation is that all spherical
         * boundary components will be filled in with balls.
         *
         * This operation is performed by attaching a new tetrahedron to
         * each boundary face and then gluing these new tetrahedra
         * together in a way that mirrors the adjacencies of the
         * underlying boundary faces.  Each boundary component will
         * thereby be pushed up through the new tetrahedra and converted
         * into a cusp formed using vertices of these new tetrahedra.
         *
         * Note that this operation is a loose converse of idealToFinite().
         *
         * \warning If a real boundary component contains vertices whose
         * links are not discs, this operation may have unexpected results.
         *
         * @return \c true if changes were made, or \c false if the
         * original triangulation contained no real boundary components.
         */
        bool finiteToIdeal();

        /**
         * Does a barycentric subdivision of the triangulation.
         * Each tetrahedron is divided into 24 tetrahedra by placing
         * an extra vertex at the centroid of each tetrahedron, the
         * centroid of each face and the midpoint of each edge.
         *
         * @author David Letscher
         */
        void barycentricSubdivision();

        /**
         * Drills out a regular neighbourhood of the given edge of the
         * triangulation.
         *
         * This is done by (i) performing two barycentric subdivisions,
         * (ii) removing all tetrahedra that touch the original edge,
         * and (iii) simplifying the resulting triangulation.
         *
         * \warning The second barycentric subdivision will multiply the
         * number of tetrahedra by 576; as a result this routine might
         * be slow, and the number of tetrahedra at the end might be
         * large (even taking the simplification into account).
         *
         * @param e the edge to drill out.
         */
        void drillEdge(NEdge* e);

        /*@}*/
        /**
         * (end: Subdivisions and Covers)
         */

        /**
         * \name Building Triangulations
         */
        /*@{*/

        /**
         * Performs a layering upon the given boundary edge of the
         * triangulation.  See the NLayering class notes for further
         * details on what a layering entails.
         *
         * \pre The given edge is a boundary edge of this triangulation,
         * and the two boundary faces on either side of it are distinct.
         *
         * @param edge the boundary edge upon which to layer.
         * @return the new tetrahedron provided by the layering.
         */
        NTetrahedron* layerOn(NEdge* edge);

        /**
         * Inserts a new layered solid torus into the triangulation.
         * The meridinal disc of the layered solid torus will intersect
         * the three edges of the boundary torus in \a cuts0, \a cuts1
         * and (\a cuts0 + \a cuts1) points respectively.
         *
         * The boundary torus will always consist of faces 012 and 013 of the
         * tetrahedron containing this boundary torus (this tetrahedron will be
         * returned).  In face 012, edges 12, 02 and 01 will meet the meridinal
         * disc \a cuts0, \a cuts1 and (\a cuts0 + \a cuts1) times respectively.
         * The only exceptions are if these three intersection numbers are
         * (1,1,2) or (0,1,1), in which case edges 12, 02 and 01 will meet the
         * meridinal disc (1, 2 and 1) or (1, 1 and 0) times respectively.
         *
         * The new tetrahedra will be inserted at the end of the list of
         * tetrahedra in the triangulation.
         *
         * \pre 0 \<= \a cuts0 \<= \a cuts1;
         * \pre \a cuts1 is non-zero;
         * \pre gcd(\a cuts0, \a cuts1) = 1.
         *
         * @param cuts0 the smallest of the three desired intersection numbers.
         * @param cuts1 the second smallest of the three desired intersection
         * numbers.
         * @return the tetrahedron containing the boundary torus.
         *
         * @see NLayeredSolidTorus
         */
        NTetrahedron* insertLayeredSolidTorus(unsigned long cuts0,
            unsigned long cuts1);
        /**
         * Inserts a new layered lens space L(p,q) into the triangulation.
         * The lens space will be created by gluing together two layered
         * solid tori in a way that uses the fewest possible tetrahedra.
         *
         * The new tetrahedra will be inserted at the end of the list of
         * tetrahedra in the triangulation.
         *
         * \pre \a p \> \a q \>= 0 unless (<i>p</i>,<i>q</i>) = (0,1);
         * \pre gcd(\a p, \a q) = 1.
         *
         * @param p a parameter of the desired lens space.
         * @param q a parameter of the desired lens space.
         *
         * @see NLayeredLensSpace
         */
        void insertLayeredLensSpace(unsigned long p, unsigned long q);
        /**
         * Inserts a layered loop of the given length into this triangulation.
         * Layered loops are described in more detail in the NLayeredLoop
         * class notes.
         *
         * The new tetrahedra will be inserted at the end of the list of
         * tetrahedra in the triangulation.
         *
         * @param length the length of the new layered loop; this must
         * be strictly positive.
         * @param twisted \c true if the new layered loop should be twisted,
         * or \c false if it should be untwisted.
         *
         * @see NLayeredLoop
         */
        void insertLayeredLoop(unsigned long length, bool twisted);
        /**
         * Inserts an augmented triangular solid torus with the given
         * parameters into this triangulation.  Almost all augmented
         * triangular solid tori represent Seifert fibred spaces with three
         * or fewer exceptional fibres.  Augmented triangular solid tori
         * are described in more detail in the NAugTriSolidTorus class notes.
         *
         * The resulting Seifert fibred space will be
         * SFS((<i>a1</i>,<i>b1</i>) (<i>a2</i>,<i>b2</i>)
         * (<i>a3</i>,<i>b3</i>) (1,1)), where the parameters
         * <i>a1</i>, ..., <i>b3</i> are passed as arguments to this
         * routine.  The three layered solid tori that are attached to
         * the central triangular solid torus will be
         * LST(|<i>a1</i>|, |<i>b1</i>|, |-<i>a1</i>-<i>b1</i>|), ...,
         * LST(|<i>a3</i>|, |<i>b3</i>|, |-<i>a3</i>-<i>b3</i>|).
         *
         * The new tetrahedra will be inserted at the end of the list of
         * tetrahedra in the triangulation.
         *
         * \pre gcd(\a a1, \a b1) = 1.
         * \pre gcd(\a a2, \a b2) = 1.
         * \pre gcd(\a a3, \a b3) = 1.
         *
         * @param a1 a parameter describing the first layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         * @param b1 a parameter describing the first layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         * @param a2 a parameter describing the second layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         * @param b2 a parameter describing the second layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         * @param a3 a parameter describing the third layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         * @param b3 a parameter describing the third layered solid
         * torus in the augmented triangular solid torus; this may be
         * either positive or negative.
         */
        void insertAugTriSolidTorus(long a1, long b1, long a2, long b2,
            long a3, long b3);
        /**
         * Inserts an orientable Seifert fibred space with at most three
         * exceptional fibres over the 2-sphere into this triangulation.
         *
         * The inserted Seifert fibred space will be
         * SFS((<i>a1</i>,<i>b1</i>) (<i>a2</i>,<i>b2</i>)
         * (<i>a3</i>,<i>b3</i>) (1,1)), where the parameters
         * <i>a1</i>, ..., <i>b3</i> are passed as arguments to this
         * routine.
         *
         * The three pairs of parameters (<i>a</i>,<i>b</i>) do not need
         * to be normalised, i.e., the parameters can be positive or
         * negative and <i>b</i> may lie outside the range [0..<i>a</i>).
         * There is no separate twisting parameter; each additional
         * twist can be incorporated into the existing parameters
         * by replacing some pair (<i>a</i>,<i>b</i>) with the pair
         * (<i>a</i>,<i>a</i>+<i>b</i>).  For Seifert fibred
         * spaces with less than three exceptional fibres, some or all
         * of the parameter pairs may be (1,<i>k</i>) or even (1,0).
         *
         * The new tetrahedra will be inserted at the end of the list of
         * tetrahedra in the triangulation.
         *
         * \pre None of \a a1, \a a2 or \a a3 are 0.
         * \pre gcd(\a a1, \a b1) = 1.
         * \pre gcd(\a a2, \a b2) = 1.
         * \pre gcd(\a a3, \a b3) = 1.
         *
         * @param a1 a parameter describing the first exceptional fibre.
         * @param b1 a parameter describing the first exceptional fibre.
         * @param a2 a parameter describing the second exceptional fibre.
         * @param b2 a parameter describing the second exceptional fibre.
         * @param a3 a parameter describing the third exceptional fibre.
         * @param b3 a parameter describing the third exceptional fibre.
         */
        void insertSFSOverSphere(long a1 = 1, long b1 = 0,
            long a2 = 1, long b2 = 0, long a3 = 1, long b3 = 0);
        /**
         * Inserts a copy of the given triangulation into this
         * triangulation.
         *
         * The new tetrahedra will be inserted into this triangulation
         * in the order in which they appear in the given triangulation,
         * and the numbering of their vertices (0-3) will not change.
         * They will be given the same descriptions as appear in the
         * given triangulation.
         *
         * @param source the triangulation whose copy will be inserted.
         */
        void insertTriangulation(const NTriangulation& source);
        /**
         * Inserts the rehydration of the given string into this triangulation.
         * If you simply wish to convert a dehydration string into a
         * new triangulation, use the static routine rehydrate() instead.
         * See dehydrate() for more information on dehydration strings.
         *
         * This routine will first rehydrate the given string into a proper
         * triangulation.  The tetrahedra from the rehydrated triangulation
         * will then be inserted into this triangulation in the same order in
         * which they appear in the rehydrated triangulation, and the
         * numbering of their vertices (0-3) will not change.
         *
         * The routine dehydrate() can be used to extract a dehydration
         * string from an existing triangulation.  Dehydration followed
         * by rehydration might not produce a triangulation identical to
         * the original, but it is guaranteed to produce an isomorphic
         * copy.  See dehydrate() for the reasons behind this.
         *
         * For a full description of the dehydrated triangulation
         * format, see <i>A Census of Cusped Hyperbolic 3-Manifolds</i>,
         * Callahan, Hildebrand and Weeks, Mathematics of Computation 68/225,
         * 1999.
         *
         * @param dehydration a dehydrated representation of the
         * triangulation to insert.  Case is irrelevant; all letters
         * will be treated as if they were lower case.
         * @return \c true if the insertion was successful, or
         * \c false if the given string could not be rehydrated.
         *
         * @see dehydrate
         * @see rehydrate
         */
        bool insertRehydration(const std::string& dehydration);
        /**
         * Rehydrates the given alphabetical string into a new triangulation.
         * See dehydrate() for more information on dehydration strings.
         *
         * This routine will rehydrate the given string into a new
         * triangulation, and return this new triangulation.
         *
         * The converse routine dehydrate() can be used to extract a
         * dehydration string from an existing triangulation.  Dehydration
         * followed by rehydration might not produce a triangulation identical
         * to the original, but it is guaranteed to produce an isomorphic
         * copy.  See dehydrate() for the reasons behind this.
         *
         * For a full description of the dehydrated triangulation
         * format, see <i>A Census of Cusped Hyperbolic 3-Manifolds</i>,
         * Callahan, Hildebrand and Weeks, Mathematics of Computation 68/225,
         * 1999.
         *
         * @param dehydration a dehydrated representation of the
         * triangulation to construct.  Case is irrelevant; all letters
         * will be treated as if they were lower case.
         * @return a newly allocated triangulation if the rehydration was
         * successful, or null if the given string could not be rehydrated.
         *
         * @see dehydrate
         * @see insertRehydration
         */
        static NTriangulation* rehydrate(const std::string& dehydration);
        /**
         * Dehydrates this triangulation into an alphabetical string.
         *
         * A <i>dehydration string</i> is a compact text representation
         * of a triangulation, introduced by Callahan, Hildebrand and Weeks
         * for their cusped hyperbolic census (see below).  The dehydration
         * string of an <i>n</i>-tetrahedron triangulation consists of
         * approximately (but not precisely) 5<i>n</i>/2 lower-case letters.
         *
         * Dehydration strings come with some restrictions:
         * - They rely on the triangulation being "canonical" in some
         *   combinatorial sense.  This is not enforced here; instead
         *   a combinatorial isomorphism is applied to make the
         *   triangulation canonical, and this isomorphic triangulation
         *   is dehydrated instead.  Note that the original triangulation
         *   is not changed.
         * - They require the triangulation to be connected.
         * - They require the triangulation to have no boundary faces
         *   (though ideal triangulations are fine).
         * - They can only support triangulations with at most 25 tetrahedra.
         *
         * The routine rehydrate() can be used to recover a
         * triangulation from a dehydration string.  Note that the
         * triangulation recovered <b>might not be identical</b> to the
         * original, but it is guaranteed to be an isomorphic copy.
         *
         * For a full description of the dehydrated triangulation
         * format, see <i>A Census of Cusped Hyperbolic 3-Manifolds</i>,
         * Callahan, Hildebrand and Weeks, Mathematics of Computation 68/225,
         * 1999.
         *
         * @return a dehydrated representation of this triangulation
         * (or an isomorphic variant of this triangulation), or the
         * empty string if dehydration is not possible because the
         * triangulation is disconnected, has boundary faces or contains
         * too many tetrahedra.
         *
         * @see rehydrate
         * @see insertRehydration
         */
        std::string dehydrate() const;
        /**
         * Constructs the isomorphism signature for this triangulation.
         *
         * An <i>isomorphism signature</i> is a compact text representation
         * of a triangulation.  Unlike dehydrations, an isomorphism signature
         * uniquely determines a triangulation up to combinatorial isomorphism.
         * That is, two triangulations are combinatorially isomorphic if
         * and only if their isomorphism signatures are the same.
         *
         * The isomorphism signature is constructed entirely of
         * printable characters, and has length proportional to
         * <tt>n log n</tt>, where \a n is the number of tetrahedra.
         *
         * Isomorphism signatures are more general than dehydrations:
         * they can be used with any triangulation (including closed, ideal,
         * bounded, invalid and/or disconnected triangulations, as well
         * as triangulations with large numbers of tetrahedra).
         *
         * The time required to construct the isomorphism signature of a
         * triangulation is <tt>O(n^2 log^2 n)</tt>.
         *
         * The routine fromIsoSig() can be used to recover a
         * triangulation from an isomorphism signature.  The triangulation
         * recovered might not be identical to the original, but it will be
         * combinatorially isomorphic.
         *
         * For a full and precise description of the isomorphism signature
         * format, see <i>Simplification paths in the Pachner graphs of
         * closed orientable 3-manifold triangulations</i>, Burton, 2011,
         * <tt>arXiv:1110.6080</tt>.
         *
         * @return the isomorphism signature of this triangulation.
         */
        std::string isoSig() const;
        /**
         * Recovers a full triangulation from an isomorphism signature.
         * See isoSig() for more information on isomorphism signatures.
         *
         * The triangulation that is returned will be newly created.
         *
         * Calling isoSig() followed by fromIsoSig() is not guaranteed to
         * produce an identical triangulation to the original, but it
         * \e is guaranteed to produce a combinatorially isomorphic
         * triangulation.
         *
         * For a full and precise description of the isomorphism signature
         * format, see <i>Simplification paths in the Pachner graphs of
         * closed orientable 3-manifold triangulations</i>, Burton, 2011,
         * <tt>arXiv:1110.6080</tt>.
         *
         * @param signature the isomorphism signature of the
         * triangulation to construct.  Note that, unlike dehydration
         * strings, case is important for isomorphism signatures.
         * @return a newly allocated triangulation if the reconstruction was
         * successful, or null if the given string was not a valid
         * isomorphism signature.
         */
        static NTriangulation* fromIsoSig(const std::string& signature);
        /**
         * Inserts into this triangulation a set of tetrahedra and their
         * gluings as described by the given integer arrays.
         *
         * This routine is provided to make it easy to hard-code a
         * medium-sized triangulation in a C++ source file.  All of the
         * pertinent data can be hard-coded into a pair of integer arrays at
         * the beginning of the source file, avoiding an otherwise tedious
         * sequence of many joinTo() calls.
         *
         * An additional \a nTetrahedra tetrahedra will be inserted into
         * this triangulation.  The relationships between these tetrahedra
         * should be stored in the two arrays as follows.  Note that the
         * new tetrahedra are numbered from 0 to (\a nTetrahedra - 1), and
         * individual tetrahedron faces are numbered from 0 to 3.
         *
         * The \a adjacencies array describes which tetrahedron faces are
         * joined to which others.  Specifically, <tt>adjacencies[t][f]</tt>
         * should contain the number of the tetrahedron joined to face \a f
         * of tetrahedron \a t.  If this face is to be left as a
         * boundary face, <tt>adjacencies[t][f]</tt> should be -1.
         *
         * The \a gluings array describes the particular gluing permutations
         * used when joining these tetrahedron faces together.  Specifically,
         * <tt>gluings[t][f][0..3]</tt> should describe the permutation
         * used to join face \a f of tetrahedron \a t to its adjacent
         * tetrahedron.  These four integers should be 0, 1, 2 and 3 in some
         * order, so that <tt>gluings[t][f][i]</tt> contains the image of
         * \a i under this permutation.  If face \a f of tetrahedron \a t
         * is to be left as a boundary face, <tt>gluings[t][f][0..3]</tt>
         * may contain anything (and will be duly ignored).
         *
         * It is the responsibility of the caller of this routine to
         * ensure that the given arrays are correct and consistent.
         * No error checking will be performed by this routine.
         *
         * Note that, for an existing triangulation, dumpConstruction()
         * will output a pair of C++ arrays that can be copied into a
         * source file and used to reconstruct the triangulation via
         * this routine.
         *
         * \ifacespython Not present.
         *
         * @param nTetrahedra the number of additional tetrahedra to insert.
         * @param adjacencies describes which of the new tetrahedron
         * faces are to be identified.  This array must have initial
         * dimension at least \a nTetrahedra.
         * @param gluings describes the specific gluing permutations by
         * which these new tetrahedron faces should be identified.  This
         * array must also have initial dimension at least \a nTetrahedra.
         */
        void insertConstruction(unsigned long nTetrahedra,
            const int adjacencies[][4], const int gluings[][4][4]);
        /**
         * Returns C++ code that can be used with insertConstruction()
         * to reconstruct this triangulation.
         *
         * The code produced will consist of the following:
         *
         * - the declaration and initialisation of two integer arrays,
         *   describing the tetrahedron gluings in this trianguation;
         * - two additional lines that declare a new NTriangulation and
         *   call insertConstruction() to rebuild this triangulation.
         *
         * The main purpose of this routine is to generate the two integer
         * arrays, which can be tedious and error-prone to code up by hand.
         *
         * Note that the number of lines of code produced grows linearly
         * with the number of tetrahedra.  If this triangulation is very
         * large, the returned string will be very large as well.
         *
         * @return the C++ code that was generated.
         */
        std::string dumpConstruction() const;

        /*@}*/
        /**
         * (end: Building Triangulations)
         */

        /**
         * Allows the user to interactively enter a triangulation in
         * plain text.  Prompts will be sent to the given output stream
         * and information will be read from the given input stream.
         *
         * \ifacespython This routine is a member of class Engine.
         * It takes no parameters; \a in and \a out are always assumed
         * to be standard input and standard output respectively.
         *
         * @param in the input stream from which text will be read.
         * @param out the output stream to which prompts will be
         * written.
         * @return the triangulation entered in by the user.
         */
        static NTriangulation* enterTextTriangulation(std::istream& in,
                std::ostream& out);

        static NXMLPacketReader* getXMLReader(NPacket* parent);
        static NTriangulation* readPacket(NFile& in, NPacket* parent);

    protected:
        virtual NPacket* internalClonePacket(NPacket* parent) const;
        virtual void writeXMLPacketData(std::ostream& out) const;

        /**
         * Turns this triangulation into a clone of the given
         * triangulation.
         * The tree structure and label of this triangulation are not
         * touched.
         *
         * @param from the triangulation from which this triangulation
         * will be cloned.
         */
        void cloneFrom(const NTriangulation& from);

    private:
        void deleteTetrahedra();
            /**< Deallocates all tetrahedra and empties the list. */
        void deleteSkeleton();
            /**< Deallocates all skeletal objects and empties all
                 corresponding lists. */

        /**
         * Clears any calculated properties and declares them all
         * unknown.  All dynamic memory used for storing known
         * properties is deallocated.
         *
         * In most cases this routine is followed immediately by firing
         * a packet change event.
         */
        virtual void clearAllProperties();

        /**
         * Checks that the permutations on face gluings are valid and
         * that adjacent face gluings match.  That is, if face \a A of
         * tetrahedron \a X is glued to face \a B of tetrahedron \a Y,
         * then face \a B of tetrahedron \a Y is also glued to face \a A
         * of tetrahedron \a X using the inverse permutation.
         *
         * If any of these tests fail, this routine writes a message to
         * standard error and marks this triangulations as invalid.
         * Any such failure indicates a likely programming error, since
         * adjacent face gluings should always match if the NTetrahedron
         * gluing routines have been used correctly.
         *
         * @author Matthias Goerner
         */
        void checkPermutations() const;
        /**
         * Recalculates vertices, edges, faces, components and
         * boundary components, as well as various other skeletal
         * properties such as validity and vertex links.
         * All appropriate lists are filled.
         *
         * \pre All skeletal lists are empty.
         */
        void calculateSkeleton() const;
        /**
         * Calculates the triangulation components and associated
         * properties.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateComponents() const;
        void labelComponent(NTetrahedron*, NComponent*) const;
            /**< Internal to calculateComponents(). */
        /**
         * Calculates the triangulation vertices and associated
         * properties.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateVertices() const;
        void labelVertex(NTetrahedron*, int, NVertex*) const;
            /**< Internal to calculateVertices(). */
        /**
         * Calculates the triangulation edges and associated
         * properties.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateEdges() const;
        void labelEdge(NTetrahedron*, int, NEdge*) const;
            /**< Internal to calculateEdges(). */
        /**
         * Calculates the triangulation faces and associated
         * properties.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateFaces() const;
        /**
         * Calculates the triangulation boundary components and
         * properties of these boundary components.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateBoundary() const;
        void labelBoundaryFace(NFace*, NBoundaryComponent*) const;
            /**< Internal to calculateBoundary(). */
        /**
         * Calculates the triangulation vertex links and associated
         * properties.
         *
         * \warning This should only be called from within
         * calculateSkeleton().
         */
        void calculateVertexLinks() const;

        /**
         * Calculates all properties of the triangulation relating to
         * its boundary components.
         */
        void calculateBoundaryProperties() const;

        /**
         * Determines if an isomorphic copy of this triangulation is
         * contained within the given triangulation.
         *
         * If the argument \a completeIsomorphism is \c true, the
         * isomorphism must be onto and boundary complete.
         * That is, this triangulation must be combinatorially
         * isomorphic to the given triangulation.
         *
         * If the argument \a completeIsomorphism is \c false, the
         * isomorphism may be boundary incomplete and may or may not be
         * onto.  That is, this triangulation must appear as a
         * subcomplex of the given triangulation, possibly with some
         * original boundary faces joined to new tetrahedra.
         *
         * See the NIsomorphism class notes for further details
         * regarding boundary complete and boundary incomplete
         * isomorphisms.
         *
         * The isomorphisms found, if any, will be appended to the
         * list \a results.  This list will not be emptied before
         * calculations begin.  All isomorphisms will be newly created,
         * and the caller of this routine is responsible for destroying
         * them.
         *
         * If \a firstOnly is passed as \c true, only the first
         * isomorphism found (if any) will be returned, after which the
         * routine will return immediately.  Otherwise all isomorphisms
         * will be returned.
         *
         * @param other the triangulation in which to search for an
         * isomorphic copy of this triangulation.
         * @param results the list in which any isomorphisms found will
         * be stored.
         * @param completeIsomorphism \c true if isomorphisms must be
         * onto and boundary complete, or \c false if neither of these
         * restrictions should be imposed.
         * @param firstOnly \c true if only one isomorphism should be
         * returned (if any), or \c false if all isomorphisms should be
         * returned.
         * @return the total number of isomorphisms found.
         */
        unsigned long findIsomorphisms(const NTriangulation& other,
                std::list<NIsomorphism*>& results,
                bool completeIsomorphism, bool firstOnly) const;

        /**
         * Internal to findIsomorphisms().
         *
         * Examines properties of the given tetrahedra to find any
         * immediate evidence that \a src may not map to \a dest in a
         * boundary complete isomorphism (in which the vertices of \a src
         * are mapped to the vertices of \a dest according to the
         * permutation \a p).
         *
         * In particular, properties such as edge degrees and vertex links
         * are examined.
         *
         * @param src the first of the two tetrahedra to examine.
         * @param dest the second of the two tetrahedra to examine.
         * @param p the permutation under which the vertices of \a src
         * must map to the vertices of \a dest.
         * @return \c true if no immediate incompatibilities between the
         * tetrahedra were found, or \c false if properties of the
         * tetrahedra were found that differ between \a src and \a dest.
         */
        static bool compatibleTets(NTetrahedron* src, NTetrahedron* dest,
                NPerm4 p);

        /**
         * Internal to isoSig().
         *
         * Constructs a candidate isomorphism signature for a single
         * component of this triangulation.  This candidate signature
         * assumes that the given tetrahedron with the given labelling
         * of its vertices becomes tetrahedron zero with vertices 0,1,2,3
         * under the "canonical isomorphism".
         *
         * @param tet the index of some tetrahedron in this triangulation.
         * @param vertices some ordering of the four vertices of the
         * given tetrahedron.
         * @return the candidate isomorphism signature.
         */
        std::string isoSig(unsigned tet, const NPerm4& vertices) const;

        /**
         * Calculates all properties that can be deduced from an
         * examination of normal surfaces in standard tri-quad coordinates.
         */
        void calculateStandardSurfaceProperties();
        /**
         * Calculates all properties that can be deduced from an
         * examination of normal surfaces in quadrilateral-only coordinates.
         */
        void calculateQuadSurfaceProperties();

        void stretchBoundaryForestFromVertex(NVertex*,
                stdhash::hash_set<NEdge*, HashPointer>&,
                stdhash::hash_set<NVertex*, HashPointer>&) const;
            /**< Internal to maximalForestInBoundary(). */
        bool stretchForestFromVertex(NVertex*,
                stdhash::hash_set<NEdge*, HashPointer>&,
                stdhash::hash_set<NVertex*, HashPointer>&,
                stdhash::hash_set<NVertex*, HashPointer>&) const;
            /**< Internal to maximalForestInSkeleton(). */
        void stretchDualForestFromTet(NTetrahedron*,
                stdhash::hash_set<NFace*, HashPointer>&,
                stdhash::hash_set<NTetrahedron*, HashPointer>&) const;
            /**< Internal to maximalForestInDualSkeleton(). */

    friend class regina::NTetrahedron;
    friend class regina::NXMLTriangulationReader;
};

/*@}*/

} // namespace regina
// Some more headers that are required for inline functions:
#include "triangulation/ntetrahedron.h"
#include "triangulation/nface.h"
#include "triangulation/nedge.h"
#include "triangulation/nvertex.h"
#include "triangulation/ncomponent.h"
#include "triangulation/nboundarycomponent.h"
namespace regina {

// Inline functions for NTriangulation

inline NTriangulation::NTriangulation() : calculatedSkeleton(false) {
}

inline NTriangulation::NTriangulation(const NTriangulation& cloneMe) :
        NPacket(), NFilePropertyReader(), calculatedSkeleton(false) {
    cloneFrom(cloneMe);
}

inline NTriangulation::~NTriangulation() {
    clearAllProperties();
    deleteTetrahedra();
}

inline NPacket* NTriangulation::internalClonePacket(NPacket*) const {
    return new NTriangulation(*this);
}

inline bool NTriangulation::dependsOnParent() const {
    return false;
}

inline unsigned long NTriangulation::getNumberOfTetrahedra() const {
    return tetrahedra.size();
}

inline unsigned long NTriangulation::getNumberOfSimplices() const {
    return tetrahedra.size();
}

inline NTetrahedron* NTriangulation::getTetrahedron(unsigned long index) {
    return tetrahedra[index];
}

inline NTetrahedron* NTriangulation::getSimplex(unsigned long index) {
    return tetrahedra[index];
}

inline const NTetrahedron* NTriangulation::getTetrahedron(unsigned long index)
        const {
    return tetrahedra[index];
}

inline const NTetrahedron* NTriangulation::getSimplex(unsigned long index)
        const {
    return tetrahedra[index];
}

inline long NTriangulation::tetrahedronIndex(const NTetrahedron* tet) const {
    return tet->markedIndex();
}

inline long NTriangulation::simplexIndex(const NTetrahedron* tet) const {
    return tet->markedIndex();
}

inline NTetrahedron* NTriangulation::newTetrahedron() {
    ChangeEventSpan span(this);

    NTetrahedron* tet = new NTetrahedron();
    tet->tri = this;
    tetrahedra.push_back(tet);
    clearAllProperties();

    return tet;
}

inline NTetrahedron* NTriangulation::newSimplex() {
    return newTetrahedron();
}

inline NTetrahedron* NTriangulation::newTetrahedron(const std::string& desc) {
    ChangeEventSpan span(this);

    NTetrahedron* tet = new NTetrahedron(desc);
    tet->tri = this;
    tetrahedra.push_back(tet);
    clearAllProperties();

    return tet;
}

inline NTetrahedron* NTriangulation::newSimplex(const std::string& desc) {
    return newTetrahedron(desc);
}

inline long NTriangulation::getTetrahedronIndex(const NTetrahedron* tet) const {
    TetrahedronIterator pos =
        std::find(tetrahedra.begin(), tetrahedra.end(), tet);
    return (pos == tetrahedra.end() ? -1 : pos - tetrahedra.begin());
}

inline void NTriangulation::removeTetrahedronAt(unsigned long index) {
    ChangeEventSpan span(this);

    NTetrahedron* ans = tetrahedra[index];
    ans->isolate();
    tetrahedra.erase(tetrahedra.begin() + index);
    delete ans;

    clearAllProperties();
}

inline void NTriangulation::removeSimplexAt(unsigned long index) {
    removeTetrahedronAt(index);
}

inline void NTriangulation::removeTetrahedron(NTetrahedron* tet) {
    ChangeEventSpan span(this);

    tet->isolate();
    tetrahedra.erase(tetrahedra.begin() + tetrahedronIndex(tet));
    delete tet;

    clearAllProperties();
}

inline void NTriangulation::removeSimplex(NTetrahedron* tet) {
    removeTetrahedron(tet);
}

inline void NTriangulation::removeAllTetrahedra() {
    ChangeEventSpan span(this);
    deleteTetrahedra();
    clearAllProperties();
}

inline void NTriangulation::removeAllSimplices() {
    removeAllTetrahedra();
}

inline void NTriangulation::gluingsHaveChanged() {
}

inline unsigned long NTriangulation::getNumberOfBoundaryComponents() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return boundaryComponents.size();
}

inline unsigned long NTriangulation::getNumberOfComponents() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return components.size();
}

inline unsigned long NTriangulation::getNumberOfVertices() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return vertices.size();
}

inline unsigned long NTriangulation::getNumberOfEdges() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return edges.size();
}

inline unsigned long NTriangulation::getNumberOfFaces() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return faces.size();
}

inline long NTriangulation::getEulerCharTri() const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    // Cast away the unsignedness of std::vector::size().
    return static_cast<long>(vertices.size())
        - static_cast<long>(edges.size())
        + static_cast<long>(faces.size())
        - static_cast<long>(tetrahedra.size());
}

inline long NTriangulation::getEulerCharacteristic() const {
    return getEulerCharTri();
}

inline const std::vector<NTetrahedron*>& NTriangulation::getTetrahedra() const {
    return (const std::vector<NTetrahedron*>&)(tetrahedra);
}

inline const std::vector<NTetrahedron*>& NTriangulation::getSimplices() const {
    return (const std::vector<NTetrahedron*>&)(tetrahedra);
}

inline const std::vector<NBoundaryComponent*>&
        NTriangulation::getBoundaryComponents() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (const std::vector<NBoundaryComponent*>&)(boundaryComponents);
}

inline const std::vector<NComponent*>& NTriangulation::getComponents() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (const std::vector<NComponent*>&)(components);
}

inline const std::vector<NVertex*>& NTriangulation::getVertices() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (const std::vector<NVertex*>&)(vertices);
}

inline const std::vector<NEdge*>& NTriangulation::getEdges()
        const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (const std::vector<NEdge*>&)(edges);
}

inline const std::vector<NFace*>& NTriangulation::getFaces()
        const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (const std::vector<NFace*>&)(faces);
}

inline NComponent* NTriangulation::getComponent(unsigned long index) const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return components[index];
}

inline NBoundaryComponent* NTriangulation::getBoundaryComponent(
        unsigned long index) const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return boundaryComponents[index];
}

inline NVertex* NTriangulation::getVertex(unsigned long index) const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return vertices[index];
}

inline NEdge* NTriangulation::getEdge(unsigned long index) const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return edges[index];
}

inline NFace* NTriangulation::getFace(unsigned long index) const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return faces[index];
}

inline long NTriangulation::componentIndex(const NComponent* component) const {
    return component->markedIndex();
}

inline long NTriangulation::getComponentIndex(
        const NComponent* component) const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    ComponentIterator pos = std::find(components.begin(), components.end(),
        component);
    return (pos == components.end() ? -1 : pos - components.begin());
}

inline long NTriangulation::boundaryComponentIndex(
        const NBoundaryComponent* boundaryComponent) const {
    return boundaryComponent->markedIndex();
}

inline long NTriangulation::getBoundaryComponentIndex(
        const NBoundaryComponent* boundaryComponent) const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    BoundaryComponentIterator pos = std::find(boundaryComponents.begin(),
        boundaryComponents.end(), boundaryComponent);
    return (pos == boundaryComponents.end() ? -1 :
        pos - boundaryComponents.begin());
}

inline long NTriangulation::vertexIndex(const NVertex* vertex) const {
    return vertex->markedIndex();
}

inline long NTriangulation::getVertexIndex(const NVertex* vertex) const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    VertexIterator pos = std::find(vertices.begin(), vertices.end(), vertex);
    return (pos == vertices.end() ? -1 : pos - vertices.begin());
}

inline long NTriangulation::edgeIndex(const NEdge* edge) const {
    return edge->markedIndex();
}

inline long NTriangulation::getEdgeIndex(const NEdge* edge) const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    EdgeIterator pos = std::find(edges.begin(), edges.end(), edge);
    return (pos == edges.end() ? -1 : pos - edges.begin());
}

inline long NTriangulation::faceIndex(const NFace* face) const {
    return face->markedIndex();
}

inline long NTriangulation::getFaceIndex(const NFace* face) const {
    if (! calculatedSkeleton)
        calculateSkeleton();

    FaceIterator pos = std::find(faces.begin(), faces.end(), face);
    return (pos == faces.end() ? -1 : pos - faces.begin());
}

inline bool NTriangulation::hasTwoSphereBoundaryComponents() const {
    if (! twoSphereBoundaryComponents.known())
        calculateBoundaryProperties();
    return twoSphereBoundaryComponents.value();
}

inline bool NTriangulation::hasNegativeIdealBoundaryComponents() const {
    if (! negativeIdealBoundaryComponents.known())
        calculateBoundaryProperties();
    return negativeIdealBoundaryComponents.value();
}

inline bool NTriangulation::isValid() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return valid;
}

inline bool NTriangulation::isIdeal() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return ideal;
}

inline bool NTriangulation::isStandard() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return standard;
}

inline bool NTriangulation::hasBoundaryFaces() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (faces.size() > 2 * tetrahedra.size());
}

inline bool NTriangulation::isClosed() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return boundaryComponents.empty();
}

inline bool NTriangulation::isOrientable() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return orientable;
}

inline bool NTriangulation::isConnected() const {
    if (! calculatedSkeleton)
        calculateSkeleton();
    return (components.size() <= 1);
}

inline void NTriangulation::simplifiedFundamentalGroup(
        NGroupPresentation* newGroup) {
    fundamentalGroup = newGroup;
}

inline bool NTriangulation::knowsZeroEfficient() const {
    return zeroEfficient.known();
}

inline bool NTriangulation::knowsSplittingSurface() const {
    return splittingSurface.known();
}

inline unsigned long NTriangulation::getHomologyH2Z2() const {
    return getHomologyH1Rel().getRank() + getHomologyH1Rel().getTorsionRank(2);
}

inline const NTriangulation::TuraevViroSet&
        NTriangulation::allCalculatedTuraevViro() const {
    return turaevViroCache;
}

inline void NTriangulation::writeTextShort(std::ostream& out) const {
    out << "Triangulation with " << tetrahedra.size() << " tetrahedra.";
}

} // namespace regina

#endif

