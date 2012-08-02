
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
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

#include <boost/python.hpp>
#include "algebra/ngrouppresentation.h"
#include "triangulation/nisomorphism.h"
#include "triangulation/ntriangulation.h"

using namespace boost::python;
using regina::NTriangulation;

namespace {
    regina::NTetrahedron* (NTriangulation::*newTetrahedron_void)() =
        &NTriangulation::newTetrahedron;
    regina::NTetrahedron* (NTriangulation::*newTetrahedron_string)(
        const std::string&) = &NTriangulation::newTetrahedron;
    regina::NTetrahedron* (NTriangulation::*getTetrahedron_non_const)(
        unsigned long) = &NTriangulation::getTetrahedron;
    bool (NTriangulation::*twoZeroMove_vertex)(regina::NVertex*, bool, bool) =
        &NTriangulation::twoZeroMove;
    bool (NTriangulation::*twoZeroMove_edge)(regina::NEdge*, bool, bool) =
        &NTriangulation::twoZeroMove;
    std::string (NTriangulation::*isoSig_void)() const =
        &NTriangulation::isoSig;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_simplifyToLocalMinimum,
        NTriangulation::simplifyToLocalMinimum, 0, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_threeTwoMove,
        NTriangulation::threeTwoMove, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_twoThreeMove,
        NTriangulation::twoThreeMove, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_fourFourMove,
        NTriangulation::fourFourMove, 2, 4);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_twoZeroMove,
        NTriangulation::twoZeroMove, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_twoOneMove,
        NTriangulation::twoOneMove, 2, 4);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_openBook,
        NTriangulation::openBook, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_closeBook,
        NTriangulation::closeBook, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_shellBoundary,
        NTriangulation::shellBoundary, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_collapseEdge,
        NTriangulation::collapseEdge, 1, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_reorderTetrahedraBFS,
        NTriangulation::reorderTetrahedraBFS, 0, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_order,
        NTriangulation::order, 0, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_splitIntoComponents,
        NTriangulation::splitIntoComponents, 0, 2);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_connectedSumDecomposition,
        NTriangulation::connectedSumDecomposition, 0, 2);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_idealToFinite,
        NTriangulation::idealToFinite, 0, 1);

    void simplifiedFundamentalGroup_own(NTriangulation& tri,
            std::auto_ptr<regina::NGroupPresentation> group) {
        tri.simplifiedFundamentalGroup(group.release());
    }

    NTriangulation* enterTextTriangulation_stdio() {
        return NTriangulation::enterTextTriangulation(std::cin, std::cout);
    }

    void addTetrahedron_own(NTriangulation& tri,
            std::auto_ptr<regina::NTetrahedron> tet) {
        tri.addTetrahedron(tet.get());
        tet.release();
    }

    regina::NIsomorphism* isIsomorphicTo_ptr(NTriangulation& t,
            NTriangulation& s) {
        return t.isIsomorphicTo(s).release();
    }

    regina::NIsomorphism* isContainedIn_ptr(NTriangulation& t,
            NTriangulation& s) {
        return t.isContainedIn(s).release();
    }

    boost::python::list getTetrahedra_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::TetrahedronIterator it =
                t.getTetrahedra().begin(); it != t.getTetrahedra().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }

    boost::python::list getComponents_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::ComponentIterator it =
                t.getComponents().begin(); it != t.getComponents().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }

    boost::python::list getBoundaryComponents_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::BoundaryComponentIterator it =
                t.getBoundaryComponents().begin();
                it != t.getBoundaryComponents().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }

    boost::python::list getVertices_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::VertexIterator it =
                t.getVertices().begin(); it != t.getVertices().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }

    boost::python::list getEdges_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::EdgeIterator it =
                t.getEdges().begin(); it != t.getEdges().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }

    boost::python::list getFaces_list(NTriangulation& t) {
        boost::python::list ans;
        for (NTriangulation::FaceIterator it =
                t.getFaces().begin(); it != t.getFaces().end(); it++)
            ans.append(boost::python::ptr(*it));
        return ans;
    }
}

void addNTriangulation() {
    scope s = class_<NTriangulation, bases<regina::NPacket>,
            std::auto_ptr<NTriangulation>,
            boost::noncopyable>("NTriangulation")
        .def(init<const NTriangulation&>())
        .def("getNumberOfTetrahedra", &NTriangulation::getNumberOfTetrahedra)
        .def("getNumberOfSimplices", &NTriangulation::getNumberOfSimplices)
        .def("getTetrahedra", getTetrahedra_list)
        .def("getSimplices", getTetrahedra_list)
        .def("getTetrahedron", getTetrahedron_non_const,
            return_value_policy<reference_existing_object>())
        .def("getSimplex", getTetrahedron_non_const,
            return_value_policy<reference_existing_object>())
        .def("tetrahedronIndex", &NTriangulation::tetrahedronIndex)
        .def("simplexIndex", &NTriangulation::simplexIndex)
        .def("getTetrahedronIndex", &NTriangulation::getTetrahedronIndex)
        .def("newTetrahedron", newTetrahedron_void,
            return_value_policy<reference_existing_object>())
        .def("newSimplex", newTetrahedron_void,
            return_value_policy<reference_existing_object>())
        .def("newTetrahedron", newTetrahedron_string,
            return_value_policy<reference_existing_object>())
        .def("newSimplex", newTetrahedron_string,
            return_value_policy<reference_existing_object>())
        .def("addTetrahedron", addTetrahedron_own)
        .def("removeTetrahedron", &NTriangulation::removeTetrahedron)
        .def("removeSimplex", &NTriangulation::removeSimplex)
        .def("removeTetrahedronAt", &NTriangulation::removeTetrahedronAt)
        .def("removeSimplexAt", &NTriangulation::removeSimplexAt)
        .def("removeAllTetrahedra", &NTriangulation::removeAllTetrahedra)
        .def("removeAllSimplices", &NTriangulation::removeAllSimplices)
        .def("swapContents", &NTriangulation::swapContents)
        .def("moveContentsTo", &NTriangulation::moveContentsTo)
        .def("gluingsHaveChanged", &NTriangulation::gluingsHaveChanged)
        .def("getNumberOfComponents", &NTriangulation::getNumberOfComponents)
        .def("getNumberOfBoundaryComponents",
            &NTriangulation::getNumberOfBoundaryComponents)
        .def("getNumberOfVertices", &NTriangulation::getNumberOfVertices)
        .def("getNumberOfEdges", &NTriangulation::getNumberOfEdges)
        .def("getNumberOfFaces", &NTriangulation::getNumberOfFaces)
        .def("getComponents", getComponents_list)
        .def("getBoundaryComponents", getBoundaryComponents_list)
        .def("getVertices", getVertices_list)
        .def("getEdges", getEdges_list)
        .def("getFaces", getFaces_list)
        .def("getComponent", &NTriangulation::getComponent,
            return_value_policy<reference_existing_object>())
        .def("getBoundaryComponent", &NTriangulation::getBoundaryComponent,
            return_value_policy<reference_existing_object>())
        .def("getVertex", &NTriangulation::getVertex,
            return_value_policy<reference_existing_object>())
        .def("getEdge", &NTriangulation::getEdge,
            return_value_policy<reference_existing_object>())
        .def("getFace", &NTriangulation::getFace,
            return_value_policy<reference_existing_object>())
        .def("componentIndex", &NTriangulation::componentIndex)
        .def("getComponentIndex", &NTriangulation::getComponentIndex)
        .def("boundaryComponentIndex",
            &NTriangulation::boundaryComponentIndex)
        .def("getBoundaryComponentIndex",
            &NTriangulation::getBoundaryComponentIndex)
        .def("vertexIndex", &NTriangulation::vertexIndex)
        .def("getVertexIndex", &NTriangulation::getVertexIndex)
        .def("edgeIndex", &NTriangulation::edgeIndex)
        .def("getEdgeIndex", &NTriangulation::getEdgeIndex)
        .def("faceIndex", &NTriangulation::faceIndex)
        .def("getFaceIndex", &NTriangulation::getFaceIndex)
        .def("isIsomorphicTo", isIsomorphicTo_ptr,
            return_value_policy<manage_new_object>())
        .def("isContainedIn", isContainedIn_ptr,
            return_value_policy<manage_new_object>())
        .def("hasTwoSphereBoundaryComponents",
            &NTriangulation::hasTwoSphereBoundaryComponents)
        .def("hasNegativeIdealBoundaryComponents",
            &NTriangulation::hasNegativeIdealBoundaryComponents)
        .def("getEulerCharTri", &NTriangulation::getEulerCharTri)
        .def("getEulerCharManifold", &NTriangulation::getEulerCharManifold)
        .def("getEulerCharacteristic", &NTriangulation::getEulerCharacteristic)
        .def("isValid", &NTriangulation::isValid)
        .def("isIdeal", &NTriangulation::isIdeal)
        .def("isStandard", &NTriangulation::isStandard)
        .def("hasBoundaryFaces", &NTriangulation::hasBoundaryFaces)
        .def("isClosed", &NTriangulation::isClosed)
        .def("isOrientable", &NTriangulation::isOrientable)
        .def("isOriented", &NTriangulation::isOriented)
        .def("isOrdered", &NTriangulation::isOrdered)
        .def("isConnected", &NTriangulation::isConnected)
        .def("getFundamentalGroup", &NTriangulation::getFundamentalGroup,
            return_internal_reference<>())
        .def("simplifiedFundamentalGroup", simplifiedFundamentalGroup_own)
        .def("getHomologyH1", &NTriangulation::getHomologyH1,
            return_internal_reference<>())
        .def("getHomologyH1Rel", &NTriangulation::getHomologyH1Rel,
            return_internal_reference<>())
        .def("getHomologyH1Bdry", &NTriangulation::getHomologyH1Bdry,
            return_internal_reference<>())
        .def("getHomologyH2", &NTriangulation::getHomologyH2,
            return_internal_reference<>())
        .def("getHomologyH2Z2", &NTriangulation::getHomologyH2Z2)
        .def("turaevViro", &NTriangulation::turaevViro)
        .def("isZeroEfficient", &NTriangulation::isZeroEfficient)
        .def("knowsZeroEfficient", &NTriangulation::knowsZeroEfficient)
        .def("hasSplittingSurface", &NTriangulation::hasSplittingSurface)
        .def("knowsSplittingSurface", &NTriangulation::knowsSplittingSurface)
        .def("intelligentSimplify", &NTriangulation::intelligentSimplify)
        .def("simplifyToLocalMinimum", &NTriangulation::simplifyToLocalMinimum,
            OL_simplifyToLocalMinimum())
        .def("threeTwoMove", &NTriangulation::threeTwoMove, OL_threeTwoMove())
        .def("twoThreeMove", &NTriangulation::twoThreeMove, OL_twoThreeMove())
        .def("fourFourMove", &NTriangulation::fourFourMove, OL_fourFourMove())
        .def("twoZeroMove", twoZeroMove_vertex, OL_twoZeroMove())
        .def("twoZeroMove", twoZeroMove_edge, OL_twoZeroMove())
        .def("twoOneMove", &NTriangulation::twoOneMove, OL_twoOneMove())
        .def("openBook", &NTriangulation::openBook, OL_openBook())
        .def("closeBook", &NTriangulation::closeBook, OL_closeBook())
        .def("shellBoundary", &NTriangulation::shellBoundary,
            OL_shellBoundary())
        .def("collapseEdge", &NTriangulation::collapseEdge, OL_collapseEdge())
        .def("reorderTetrahedraBFS", &NTriangulation::reorderTetrahedraBFS,
            OL_reorderTetrahedraBFS())
        .def("orient", &NTriangulation::orient)
        .def("order", &NTriangulation::order, OL_order(args("force_oriented")))
        .def("splitIntoComponents", &NTriangulation::splitIntoComponents,
            OL_splitIntoComponents())
        .def("connectedSumDecomposition",
            &NTriangulation::connectedSumDecomposition,
            OL_connectedSumDecomposition())
        .def("isThreeSphere", &NTriangulation::isThreeSphere)
        .def("knowsThreeSphere", &NTriangulation::knowsThreeSphere)
        .def("isBall", &NTriangulation::isBall)
        .def("knowsBall", &NTriangulation::knowsBall)
        .def("isSolidTorus", &NTriangulation::isSolidTorus)
        .def("knowsSolidTorus", &NTriangulation::knowsSolidTorus)
        .def("hasCompressingDisc", &NTriangulation::hasCompressingDisc)
        .def("hasSimpleCompressingDisc",
            &NTriangulation::hasSimpleCompressingDisc)
        .def("makeDoubleCover", &NTriangulation::makeDoubleCover)
        .def("idealToFinite", &NTriangulation::idealToFinite,
            OL_idealToFinite())
        .def("finiteToIdeal", &NTriangulation::finiteToIdeal)
        .def("barycentricSubdivision", &NTriangulation::barycentricSubdivision)
        .def("drillEdge", &NTriangulation::drillEdge)
        .def("layerOn", &NTriangulation::layerOn,
            return_value_policy<reference_existing_object>())
        .def("insertLayeredSolidTorus",
            &NTriangulation::insertLayeredSolidTorus,
            return_value_policy<reference_existing_object>())
        .def("insertLayeredLensSpace",
            &NTriangulation::insertLayeredLensSpace)
        .def("insertLayeredLoop", &NTriangulation::insertLayeredLoop)
        .def("insertAugTriSolidTorus", &NTriangulation::insertAugTriSolidTorus)
        .def("insertSFSOverSphere", &NTriangulation::insertSFSOverSphere)
        .def("insertTriangulation", &NTriangulation::insertTriangulation)
        .def("insertRehydration", &NTriangulation::insertRehydration)
        .def("dehydrate", &NTriangulation::dehydrate)
        .def("rehydrate", &NTriangulation::rehydrate,
            return_value_policy<manage_new_object>())
        .def("isoSig", isoSig_void)
        .def("fromIsoSig", &NTriangulation::fromIsoSig,
            return_value_policy<manage_new_object>())
        .def("dumpConstruction", &NTriangulation::dumpConstruction)
        .def("enterTextTriangulation", enterTextTriangulation_stdio,
            return_value_policy<manage_new_object>())
        .staticmethod("rehydrate")
        .staticmethod("fromIsoSig")
        .staticmethod("enterTextTriangulation")
    ;

    s.attr("packetType") = NTriangulation::packetType;

    implicitly_convertible<std::auto_ptr<NTriangulation>,
        std::auto_ptr<regina::NPacket> >();
}

