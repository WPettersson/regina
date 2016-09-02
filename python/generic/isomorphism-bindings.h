
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

#include <boost/python.hpp>
#include "generic/isomorphism.h"
#include "generic/triangulation.h"
#include "maths/nperm5.h" // Specialisation needed for 4-D case.
#include "../helpers.h"
#include "../safeheldtype.h"

using namespace boost::python;
using namespace regina::python;
using regina::Isomorphism;

namespace {
    template <int dim>
    struct PyIsoHelper {
        typedef int (Isomorphism<dim>::*simpImage_const_type)(unsigned) const;

        typedef regina::NPerm<dim+1> (Isomorphism<dim>::*facetPerm_const_type)(
            unsigned) const;
    };
}

template <int dim>
void addIsomorphism(const char* name) {
    class_<Isomorphism<dim>, std::auto_ptr<Isomorphism<dim>>,
            boost::noncopyable>(name, init<const Isomorphism<dim>&>())
        .def("size", &Isomorphism<dim>::size)
        .def("getSourceSimplices", &Isomorphism<dim>::size)
        .def("simpImage", typename PyIsoHelper<dim>::simpImage_const_type(
            &Isomorphism<dim>::simpImage))
        .def("facetPerm", typename PyIsoHelper<dim>::facetPerm_const_type(
            &Isomorphism<dim>::facetPerm))
        .def("__getitem__", &Isomorphism<dim>::operator[])
        .def("isIdentity", &Isomorphism<dim>::isIdentity)
        .def("apply", &Isomorphism<dim>::apply,
            return_value_policy<to_held_type<>>())
        .def("applyInPlace", &Isomorphism<dim>::applyInPlace)
        .def("random", &Isomorphism<dim>::random,
            return_value_policy<manage_new_object>())
        .def("identity", &Isomorphism<dim>::identity,
            return_value_policy<manage_new_object>())
        .def(regina::python::add_output())
        .def(regina::python::add_eq_operators())
        .staticmethod("random")
        .staticmethod("identity")
    ;
}

