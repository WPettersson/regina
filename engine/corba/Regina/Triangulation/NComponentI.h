
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

#ifndef __NCOMPONENTI_H
#define __NCOMPONENTI_H

#include "config.h"

#ifdef __NO_INCLUDE_PATHS
    #include "ncomponent.h"
#else
    #include "engine/triangulation/ncomponent.h"
#endif

#include "NTetrahedronIDL.h"
#include "ShareableObjectI.h"

class NComponent_i : public virtual POA_Regina::Triangulation::NComponent,
        public ShareableObject_i {
	STANDARD_ENGINE_TYPEDEFS(NComponent_i, NComponent,
			Regina::Triangulation::NComponent)

    protected:
        NComponent_i(::NComponent* newCppPtr) : ShareableObject_i(newCppPtr) {
        }
    public:
        STANDARD_NEW_WRAPPER

        virtual CORBA::Boolean isIdeal();
        virtual CORBA::Boolean isOrientable();
        virtual CORBA::Boolean isClosed();
        virtual CORBA::Long getNumberOfTetrahedra();
        virtual CORBA::Long getNumberOfFaces();
        virtual CORBA::Long getNumberOfEdges();
        virtual CORBA::Long getNumberOfVertices();
        virtual CORBA::Long getNumberOfBoundaryComponents();
        virtual Regina::Triangulation::NTetrahedron_ptr getTetrahedron(
            CORBA::Long index);
        virtual Regina::Triangulation::NFace_ptr getFace(CORBA::Long index);
        virtual Regina::Triangulation::NEdge_ptr getEdge(CORBA::Long index);
        virtual Regina::Triangulation::NVertex_ptr getVertex(CORBA::Long index);
        virtual Regina::Triangulation::NBoundaryComponent_ptr
            getBoundaryComponent(CORBA::Long index);
};

#endif

