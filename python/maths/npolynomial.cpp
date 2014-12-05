
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2014, Ben Burton                                   *
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

/* end stub */

#include <boost/python.hpp>
#include "maths/npolynomial.h"
#include "maths/nrational.h"

using namespace boost::python;
using regina::NPolynomial;
using regina::NRational;

namespace {
    const regina::NRational& getItem(const NPolynomial<NRational>& p,
            size_t exp) {
        return p[exp];
    }
    void setItem(NPolynomial<NRational>& p, size_t exp,
            const regina::NRational& value) {
        p.set(exp, value);
    }

    void (NPolynomial<NRational>::*init_void)() =
        &NPolynomial<NRational>::init;
    void (NPolynomial<NRational>::*init_degree)(size_t) =
        &NPolynomial<NRational>::init;
}

void addNPolynomial() {
    scope s = class_<NPolynomial<NRational> >("NPolynomial")
        .def(init<size_t>())
        .def(init<const NPolynomial<NRational>&>())
        .def("init", init_void)
        .def("init", init_degree)
        .def("degree", &NPolynomial<NRational>::degree)
        .def("isZero", &NPolynomial<NRational>::isZero)
        .def("__getitem__", getItem, return_internal_reference<>())
        .def("__setitem__", setItem)
        .def("set", &NPolynomial<NRational>::set)
        .def(self == self)
        .def(self != self)
        .def(self *= NRational())
        .def(self /= NRational())
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(self /= self)
        .def("divisionAlg", &NPolynomial<NRational>::divisionAlg)
        .def("gcdWithCoeffs", &NPolynomial<NRational>::gcdWithCoeffs<NRational>)
        .def(self_ns::str(self))
        .def(self_ns::repr(self))
    ;
}

