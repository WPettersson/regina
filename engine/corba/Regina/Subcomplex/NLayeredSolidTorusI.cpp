
/**************************************************************************
 *                                                                        *
 *  Regina - A normal surface theory calculator                           *
 *  Computational engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2001, Ben Burton                                   *
 *  For further details contact Ben Burton (benb@acm.org).                *
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
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

#include "NLayeredSolidTorusI.h"
#include "NTetrahedronI.h"

CORBA::Long NLayeredSolidTorus_i::getNumberOfTetrahedra() {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getNumberOfTetrahedra();
}
Regina::Triangulation::NTetrahedron_ptr NLayeredSolidTorus_i::getBase() {
	return NTetrahedron_i::newWrapper(
		GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->getBase());
}
CORBA::Long NLayeredSolidTorus_i::getBaseEdge(CORBA::Long group,
		CORBA::Long index) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getBaseEdge(group, index);
}
CORBA::Long NLayeredSolidTorus_i::getBaseEdgeGroup(CORBA::Long edge) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getBaseEdgeGroup(edge);
}
CORBA::Long NLayeredSolidTorus_i::getBaseFace(CORBA::Long index) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->getBaseFace(index);
}
Regina::Triangulation::NTetrahedron_ptr NLayeredSolidTorus_i::getTopLevel() {
	return NTetrahedron_i::newWrapper(
		GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->getTopLevel());
}
CORBA::Long NLayeredSolidTorus_i::getMeridinalCuts(CORBA::Long group) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getMeridinalCuts(group);
}
CORBA::Long NLayeredSolidTorus_i::getTopEdge(CORBA::Long group,
		CORBA::Long index) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getTopEdge(group, index);
}
CORBA::Long NLayeredSolidTorus_i::getTopEdgeGroup(CORBA::Long edge) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->
		getTopEdgeGroup(edge);
}
CORBA::Long NLayeredSolidTorus_i::getTopFace(CORBA::Long index) {
	return GET_ENGINE_OBJECT(NLayeredSolidTorus, this)->getTopFace(index);
}

