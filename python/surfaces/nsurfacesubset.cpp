
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
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

#include "surfaces/nsurfacefilter.h"
#include "surfaces/nsurfacesubset.h"
#include "triangulation/ntriangulation.h"
#include "../helpers.h"
#include "../safeheldtype.h"
#include <boost/python.hpp>

using namespace boost::python;
using namespace regina::python;
using regina::NSurfaceSubset;

namespace {
    void writeAllSurfaces_stdio(const regina::NNormalSurfaceList& s) {
        s.writeAllSurfaces(std::cout);
    }
}

void addNSurfaceSubset() {
    scope s = class_<NSurfaceSubset, std::auto_ptr<NSurfaceSubset>,
            boost::noncopyable>("NSurfaceSubset",
            init<const regina::NNormalSurfaceList&,
                const regina::NSurfaceFilter&>())
        .def("getFlavour", &NSurfaceSubset::coords)
        .def("coords", &NSurfaceSubset::coords)
        .def("allowsAlmostNormal", &NSurfaceSubset::allowsAlmostNormal)
        .def("allowsSpun", &NSurfaceSubset::allowsSpun)
        .def("allowsOriented", &NSurfaceSubset::allowsOriented)
        .def("isEmbeddedOnly", &NSurfaceSubset::isEmbeddedOnly)
        .def("triangulation", &NSurfaceSubset::triangulation,
            return_value_policy<to_held_type<> >())
        .def("getTriangulation", &NSurfaceSubset::triangulation,
            return_value_policy<to_held_type<> >())
        .def("size", &NSurfaceSubset::size)
        .def("getNumberOfSurfaces", &NSurfaceSubset::size)
        .def("surface", &NSurfaceSubset::surface,
            return_internal_reference<>())
        .def("getSurface", &NSurfaceSubset::surface,
            return_internal_reference<>())
        .def("writeAllSurfaces", writeAllSurfaces_stdio)
        .def(regina::python::add_output())
        .def(regina::python::add_eq_operators())
    ;
}

