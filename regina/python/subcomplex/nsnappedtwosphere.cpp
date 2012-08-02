
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
#include "subcomplex/nsnappedball.h"
#include "subcomplex/nsnappedtwosphere.h"
#include "triangulation/ntetrahedron.h"
#include "triangulation/ntriangulation.h"

using namespace boost::python;
using regina::NSnappedTwoSphere;

namespace {
    NSnappedTwoSphere* (*formsStructure_tets)
        (regina::NTetrahedron*, regina::NTetrahedron*) =
        &NSnappedTwoSphere::formsSnappedTwoSphere;
    NSnappedTwoSphere* (*formsStructure_balls)
        (regina::NSnappedBall*, regina::NSnappedBall*) =
        &NSnappedTwoSphere::formsSnappedTwoSphere;
}

void addNSnappedTwoSphere() {
    class_<NSnappedTwoSphere, bases<regina::ShareableObject>,
            std::auto_ptr<NSnappedTwoSphere>, boost::noncopyable>
            ("NSnappedTwoSphere", no_init)
        .def("clone", &NSnappedTwoSphere::clone,
            return_value_policy<manage_new_object>())
        .def("getSnappedBall", &NSnappedTwoSphere::getSnappedBall,
            return_internal_reference<>())
        .def("formsSnappedTwoSphere", formsStructure_tets,
            return_value_policy<manage_new_object>())
        .def("formsSnappedTwoSphere", formsStructure_balls,
            return_value_policy<manage_new_object>())
        .staticmethod("formsSnappedTwoSphere")
    ;
}

