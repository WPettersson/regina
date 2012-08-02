
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include "triangulation/ntetrahedron.h"
#include "subcomplex/nsnappedtwosphere.h"
#include "subcomplex/nsnappedball.h"

namespace regina {

NSnappedTwoSphere* NSnappedTwoSphere::clone() const {
    NSnappedTwoSphere* ans = new NSnappedTwoSphere();
    ans->ball[0] = ball[0]->clone();
    ans->ball[1] = ball[1]->clone();
    return ans;
}

NSnappedTwoSphere* NSnappedTwoSphere::formsSnappedTwoSphere(
        NTetrahedron* tet1, NTetrahedron* tet2) {
    NSnappedBall* ball[2];
    if (! (ball[0] = NSnappedBall::formsSnappedBall(tet1)))
        return 0;
    if (! (ball[1] = NSnappedBall::formsSnappedBall(tet2))) {
        delete ball[0];
        return 0;
    }
    if (tet1->getEdge(ball[0]->getEquatorEdge()) !=
            tet2->getEdge(ball[1]->getEquatorEdge())) {
        delete ball[0];
        delete ball[1];
        return 0;
    }

    // This is it.
    NSnappedTwoSphere* ans = new NSnappedTwoSphere();
    ans->ball[0] = ball[0];
    ans->ball[1] = ball[1];
    return ans;
}

NSnappedTwoSphere* NSnappedTwoSphere::formsSnappedTwoSphere(
        NSnappedBall* ball1, NSnappedBall* ball2) {
    if (ball1->getTetrahedron()->getEdge(ball1->getEquatorEdge()) !=
            ball2->getTetrahedron()->getEdge(ball2->getEquatorEdge()))
        return 0;

    // This is it.
    NSnappedTwoSphere* ans = new NSnappedTwoSphere();
    ans->ball[0] = ball1->clone();
    ans->ball[1] = ball2->clone();
    return ans;
}

} // namespace regina

