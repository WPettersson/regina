
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
#include "surfaces/nprism.h"
#include "surfaces/nnormalsurface.h"

using namespace boost::python;
using regina::NPrismSpec;
using regina::NPrismSetSurface;

void addNPrism() {
    class_<NPrismSpec>("NPrismSpec")
        .def(init<unsigned long, int>())
        .def(init<const NPrismSpec&>())
        .def_readwrite("tetIndex", &NPrismSpec::tetIndex)
        .def_readwrite("edge", &NPrismSpec::edge)
        .def(self == self)
        .def(self_ns::str(self))
    ;

    class_<NPrismSetSurface, std::auto_ptr<NPrismSetSurface>,
            boost::noncopyable>("NPrismSetSurface",
            init<const regina::NNormalSurface&>())
        .def("getQuadType", &NPrismSetSurface::getQuadType)
    ;
}

