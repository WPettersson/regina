
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
#include "packet/nscript.h"

using namespace boost::python;
using regina::NScript;

namespace {
    const std::string& (NScript::*getVariableValue_long)(unsigned long)
        const = &NScript::getVariableValue;
    const std::string& (NScript::*getVariableValue_string)(const std::string&)
        const = &NScript::getVariableValue;
}

void addNScript() {
    scope s = class_<NScript, bases<regina::NPacket>,
            std::auto_ptr<NScript>, boost::noncopyable>("NScript", init<>())
        .def("getNumberOfLines", &NScript::getNumberOfLines)
        .def("getLine", &NScript::getLine,
            return_value_policy<return_by_value>())
        .def("addFirst", &NScript::addFirst)
        .def("addLast", &NScript::addLast)
        .def("insertAtPosition", &NScript::insertAtPosition)
        .def("replaceAtPosition", &NScript::replaceAtPosition)
        .def("removeLineAt", &NScript::removeLineAt)
        .def("removeAllLines", &NScript::removeAllLines)
        .def("getNumberOfVariables", &NScript::getNumberOfVariables)
        .def("getVariableName", &NScript::getVariableName,
            return_value_policy<return_by_value>())
        .def("getVariableValue", getVariableValue_long,
            return_value_policy<return_by_value>())
        .def("getVariableValue", getVariableValue_string,
            return_value_policy<return_by_value>())
        .def("addVariable", &NScript::addVariable)
        .def("removeVariable", &NScript::removeVariable)
        .def("removeAllVariables", &NScript::removeAllVariables)
    ;

    s.attr("packetType") = NScript::packetType;

    implicitly_convertible<std::auto_ptr<NScript>,
        std::auto_ptr<regina::NPacket> >();
}

