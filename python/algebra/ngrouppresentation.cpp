
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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
#include "algebra/nabeliangroup.h"
#include "algebra/ngrouppresentation.h"
#include "algebra/nhomgrouppresentation.h"
#include "algebra/nmarkedabeliangroup.h"

using namespace boost::python;
using regina::NGroupExpressionTerm;
using regina::NGroupExpression;
using regina::NGroupPresentation;

namespace {
    void (NGroupExpression::*addTermFirst_term)(const NGroupExpressionTerm&) =
        &NGroupExpression::addTermFirst;
    void (NGroupExpression::*addTermFirst_pair)(unsigned long, long) =
        &NGroupExpression::addTermFirst;
    void (NGroupExpression::*addTermLast_term)(const NGroupExpressionTerm&) =
        &NGroupExpression::addTermLast;
    void (NGroupExpression::*addTermLast_pair)(unsigned long, long) =
        &NGroupExpression::addTermLast;
    NGroupExpressionTerm& (NGroupExpression::*getTerm_non_const)(
        unsigned long) = &NGroupExpression::getTerm;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_simplify,
        NGroupExpression::simplify, 0, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_substitute,
        NGroupExpression::substitute, 2, 3);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_addGenerator,
        NGroupPresentation::addGenerator, 0, 1);

    object getTerms_list(const NGroupExpression& e) {
        boost::python::list ans;
        for (std::list<NGroupExpressionTerm>::const_iterator it =
                e.getTerms().begin(); it != e.getTerms().end(); it++)
            ans.append(*it);
        return ans;
    }

    void expressionWriteText(const NGroupExpression& e, bool sw = false) {
        e.writeText(std::cout, sw);
    }

    void expressionWriteTeX(const NGroupExpression& e) {
        e.writeTeX(std::cout);
    }

    BOOST_PYTHON_FUNCTION_OVERLOADS(OL_expressionWriteText,
        expressionWriteText, 1, 2);

    void addRelation_own(NGroupPresentation& p,
            std::auto_ptr<NGroupExpression> e) {
        p.addRelation(e.get());
        e.release();
    }

    void presentationWriteTeX(const NGroupPresentation& p) {
        p.writeTeX(std::cout);
    }

    void presentationWriteTextCompact(const NGroupPresentation& p) {
        p.writeTextCompact(std::cout);
    }

    regina::NHomGroupPresentation* intelligentSimplifyDetail_ptr(
            NGroupPresentation& p) {
        return p.intelligentSimplifyDetail().release();
    }

    regina::NAbelianGroup* abelianisation_ptr(const NGroupPresentation& p) {
        return p.abelianisation().release();
    }

    regina::NMarkedAbelianGroup* markedAbelianisation_ptr(
            const NGroupPresentation& p) {
        return p.markedAbelianisation().release();
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(OL_proliferateRelators,
        NGroupPresentation::proliferateRelators, 0, 1);
}

void addNGroupPresentation() {
    class_<NGroupExpressionTerm>("NGroupExpressionTerm")
        .def_readwrite("generator", &NGroupExpressionTerm::generator)
        .def_readwrite("exponent", &NGroupExpressionTerm::exponent)
        .def(init<unsigned long, long>())
        .def(init<const NGroupExpressionTerm&>())
        .def(self == self)
        .def("inverse", &NGroupExpressionTerm::inverse)
        .def(self += self)
        .def(self_ns::str(self))
    ;

    class_<NGroupExpression, bases<regina::ShareableObject>,
            std::auto_ptr<NGroupExpression>, boost::noncopyable>
            ("NGroupExpression")
        .def(init<const NGroupExpression&>())
        .def("getTerms", getTerms_list)
        .def("getNumberOfTerms", &NGroupExpression::getNumberOfTerms)
        .def("wordLength", &NGroupExpression::wordLength)
        .def("erase", &NGroupExpression::erase)
        .def("getTerm", getTerm_non_const, return_internal_reference<>())
        .def("getGenerator", &NGroupExpression::getGenerator)
        .def("getExponent", &NGroupExpression::getExponent)
        .def("addTermFirst", addTermFirst_term)
        .def("addTermFirst", addTermFirst_pair)
        .def("addTermLast", addTermLast_term)
        .def("addTermLast", addTermLast_pair)
        .def("addTermsFirst", &NGroupExpression::addTermsFirst)
        .def("addTermsLast", &NGroupExpression::addTermsLast)
        .def("cycleLeft", &NGroupExpression::cycleLeft)
        .def("cycleRight", &NGroupExpression::cycleRight)
        .def("inverse", &NGroupExpression::inverse,
            return_value_policy<manage_new_object>())
        .def("invert", &NGroupExpression::invert)
        .def("power", &NGroupExpression::power,
            return_value_policy<manage_new_object>())
        .def("simplify", &NGroupExpression::simplify, OL_simplify())
        .def("substitute", &NGroupExpression::substitute, OL_substitute())
        .def("toTeX", &NGroupExpression::toTeX)
        .def("writeText", expressionWriteText, OL_expressionWriteText())
        .def("writeTeX", expressionWriteTeX)
    ;

    class_<NGroupPresentation, bases<regina::ShareableObject>,
            std::auto_ptr<NGroupPresentation>, boost::noncopyable>
            ("NGroupPresentation")
        .def(init<const NGroupPresentation&>())
        .def("addGenerator", &NGroupPresentation::addGenerator,
            OL_addGenerator())
        .def("addRelation", addRelation_own)
        .def("getNumberOfGenerators",
            &NGroupPresentation::getNumberOfGenerators)
        .def("getNumberOfRelations", &NGroupPresentation::getNumberOfRelations)
        .def("getRelation", &NGroupPresentation::getRelation,
            return_internal_reference<>())
        .def("intelligentSimplify", &NGroupPresentation::intelligentSimplify)
        .def("intelligentSimplifyDetail",
            intelligentSimplifyDetail_ptr,
            return_value_policy<manage_new_object>())
        .def("proliferateRelators", &NGroupPresentation::proliferateRelators,
            OL_proliferateRelators())
        .def("recogniseGroup", &NGroupPresentation::recogniseGroup)
        .def("relatorLength", &NGroupPresentation::relatorLength)
        .def("abelianisation", abelianisation_ptr,
            return_value_policy<manage_new_object>())
        .def("markedAbelianisation", markedAbelianisation_ptr,
            return_value_policy<manage_new_object>())
        .def("toTeX", &NGroupPresentation::toTeX)
        .def("toStringCompact", &NGroupPresentation::toStringCompact)
        .def("compact", &NGroupPresentation::compact)
        .def("writeTeX", presentationWriteTeX)
        .def("writeTextCompact", presentationWriteTextCompact)
    ;
}

