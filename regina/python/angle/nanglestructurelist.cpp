
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
#include "angle/nanglestructurelist.h"
#include "progress/nprogressmanager.h"
#include "triangulation/ntriangulation.h"

using namespace boost::python;
using regina::NAngleStructureList;

namespace {
    // Write manual overload wrappers since this is a static member function.
    NAngleStructureList* enumerate_1(regina::NTriangulation* owner) {
        return NAngleStructureList::enumerate(owner);
    }
    NAngleStructureList* enumerate_2(regina::NTriangulation* owner,
            bool tautOnly) {
        return NAngleStructureList::enumerate(owner, tautOnly);
    }
    NAngleStructureList* enumerate_3(regina::NTriangulation* owner,
            bool tautOnly, regina::NProgressManager* manager) {
        return NAngleStructureList::enumerate(owner, tautOnly, manager);
    }
}

void addNAngleStructureList() {
    scope s = class_<NAngleStructureList, bases<regina::NPacket>,
            std::auto_ptr<NAngleStructureList>, boost::noncopyable>
            ("NAngleStructureList", no_init)
        .def("getTriangulation", &NAngleStructureList::getTriangulation,
            return_value_policy<reference_existing_object>())
        .def("isTautOnly", &NAngleStructureList::isTautOnly)
        .def("getNumberOfStructures",
            &NAngleStructureList::getNumberOfStructures)
        .def("getStructure", &NAngleStructureList::getStructure,
            return_internal_reference<>())
        .def("spansStrict", &NAngleStructureList::spansStrict)
        .def("spansTaut", &NAngleStructureList::spansTaut)
        .def("allowsStrict", &NAngleStructureList::allowsStrict)
        .def("allowsTaut", &NAngleStructureList::allowsTaut)
        .def("enumerate", enumerate_1,
            return_value_policy<reference_existing_object>())
        .def("enumerate", enumerate_2,
            return_value_policy<reference_existing_object>())
        .def("enumerate", enumerate_3,
            return_value_policy<reference_existing_object>())
        .staticmethod("enumerate")
    ;

    s.attr("packetType") = NAngleStructureList::packetType;

    implicitly_convertible<std::auto_ptr<NAngleStructureList>,
        std::auto_ptr<regina::NPacket> >();
}

