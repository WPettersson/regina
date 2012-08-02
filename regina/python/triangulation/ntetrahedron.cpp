
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
#include "triangulation/ncomponent.h"
#include "triangulation/nedge.h"
#include "triangulation/nface.h"
#include "triangulation/ntetrahedron.h"
#include "triangulation/ntriangulation.h"
#include "triangulation/nvertex.h"
#include "../globalarray.h"

using namespace boost::python;
using regina::NTetrahedron;

void addNTetrahedron() {
    class_<NTetrahedron, bases<regina::ShareableObject>,
            std::auto_ptr<NTetrahedron>, boost::noncopyable>(
            "NTetrahedron", init<>())
        .def(init<const std::string&>())
        .def("getDescription", &NTetrahedron::getDescription,
            return_value_policy<return_by_value>())
        .def("setDescription", &NTetrahedron::setDescription)
        .def("adjacentTetrahedron", &NTetrahedron::adjacentTetrahedron,
            return_value_policy<reference_existing_object>())
        .def("adjacentSimplex", &NTetrahedron::adjacentSimplex,
            return_value_policy<reference_existing_object>())
        .def("getAdjacentTetrahedron", &NTetrahedron::getAdjacentTetrahedron,
            return_value_policy<reference_existing_object>())
        .def("adjacentGluing", &NTetrahedron::adjacentGluing)
        .def("getAdjacentTetrahedronGluing",
            &NTetrahedron::getAdjacentTetrahedronGluing)
        .def("adjacentFace", &NTetrahedron::adjacentFace)
        .def("adjacentFacet", &NTetrahedron::adjacentFacet)
        .def("getAdjacentFace", &NTetrahedron::getAdjacentFace)
        .def("hasBoundary", &NTetrahedron::hasBoundary)
        .def("joinTo", &NTetrahedron::joinTo)
        .def("unjoin", &NTetrahedron::unjoin,
            return_value_policy<reference_existing_object>())
        .def("isolate", &NTetrahedron::isolate)
        .def("getTriangulation", &NTetrahedron::getTriangulation,
            return_value_policy<reference_existing_object>())
        .def("getComponent", &NTetrahedron::getComponent,
            return_value_policy<reference_existing_object>())
        .def("getVertex", &NTetrahedron::getVertex,
            return_value_policy<reference_existing_object>())
        .def("getEdge", &NTetrahedron::getEdge,
            return_value_policy<reference_existing_object>())
        .def("getFace", &NTetrahedron::getFace,
            return_value_policy<reference_existing_object>())
        .def("getVertexMapping", &NTetrahedron::getVertexMapping)
        .def("getEdgeMapping", &NTetrahedron::getEdgeMapping)
        .def("getFaceMapping", &NTetrahedron::getFaceMapping)
        .def("orientation", &NTetrahedron::orientation)
    ;
}

