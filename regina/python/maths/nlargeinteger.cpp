
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
#include "maths/nlargeinteger.h"

using namespace boost::python;
using regina::NLargeInteger;

namespace {
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_stringValue,
        NLargeInteger::stringValue, 0, 1);

    boost::python::tuple divisionAlg(const NLargeInteger& n,
            const NLargeInteger& divisor) {
        NLargeInteger remainder;
        NLargeInteger quotient = n.divisionAlg(divisor, remainder);
        return make_tuple(quotient, remainder);
    }
}

void addNLargeInteger() {
    scope s = class_<NLargeInteger>("NLargeInteger")
        .def(init<long>())
        .def(init<const NLargeInteger&>())
        .def(init<const regina::NLazyInteger&>())
        .def(init<const char*, optional<int> >())
        .def("isInfinite", &NLargeInteger::isInfinite)
        .def("longValue", &NLargeInteger::longValue)
        .def("stringValue", &NLargeInteger::stringValue, OL_stringValue())
        .def("swap", &NLargeInteger::swap)
        .def(self == self)
        .def(self == long())
        .def(self != self)
        .def(self != long())
        .def(self < self)
        .def(self < long())
        .def(self > self)
        .def(self > long())
        .def(self <= self)
        .def(self <= long())
        .def(self >= self)
        .def(self >= long())
        .def(self + self)
        .def(self + long())
        .def(self - self)
        .def(self - long())
        .def(self * self)
        .def(self * long())
        .def(self / self)
        .def(self / long())
        .def("divExact", &NLargeInteger::divExact)
        .def(self % self)
        .def(self % long())
        .def("divisionAlg", divisionAlg)
        .def(- self)
        .def(self += self)
        .def(self += long())
        .def(self -= self)
        .def(self -= long())
        .def(self *= self)
        .def(self *= long())
        .def(self /= self)
        .def(self /= long())
        .def("divByExact", &NLargeInteger::divByExact,
            return_internal_reference<>())
        .def(self %= self)
        .def(self %= long())
        .def("negate", &NLargeInteger::negate)
        .def("raiseToPower", &NLargeInteger::raiseToPower)
        .def("abs", &NLargeInteger::abs)
        .def("gcd", &NLargeInteger::gcd)
        .def("lcm", &NLargeInteger::lcm)
        .def("gcdWithCoeffs", &NLargeInteger::gcdWithCoeffs)
        .def("legendre", &NLargeInteger::legendre)
        .def("randomBoundedByThis", &NLargeInteger::randomBoundedByThis)
        .def("randomBinary", &NLargeInteger::randomBinary)
        .def("randomCornerBinary", &NLargeInteger::randomCornerBinary)
        .def(self_ns::str(self))
        .staticmethod("randomBinary")
        .staticmethod("randomCornerBinary")
    ;

    // Apparently there is no way in python to make a module attribute
    // read-only.
    s.attr("zero") = NLargeInteger::zero;
    s.attr("one") = NLargeInteger::one;
    s.attr("infinity") = NLargeInteger::infinity;

    boost::python::implicitly_convertible<long, NLargeInteger>();
    boost::python::implicitly_convertible<std::string, NLargeInteger>();
}

