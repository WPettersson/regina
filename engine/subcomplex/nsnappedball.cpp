
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include "algebra/nabeliangroup.h"
#include "manifold/nhandlebody.h"
#include "triangulation/ntetrahedron.h"
#include "subcomplex/nsnappedball.h"

namespace regina {

NSnappedBall* NSnappedBall::clone() const {
    NSnappedBall* ans = new NSnappedBall();
    ans->tet = tet;
    ans->equator = equator;
    return ans;
}

NSnappedBall* NSnappedBall::formsSnappedBall(NTetrahedron* tet) {
    int inFace1, inFace2;
    NPerm4 perm;
    for (inFace1 = 0; inFace1 < 3; inFace1++)
        if (tet->adjacentTetrahedron(inFace1) == tet) {
            perm = tet->adjacentGluing(inFace1);
            inFace2 = perm[inFace1];
            if (perm == NPerm4(inFace1, inFace2)) {
                // This is it!
                NSnappedBall* ans = new NSnappedBall();
                ans->tet = tet;
                ans->equator = NEdge::edgeNumber[inFace1][inFace2];
                return ans;
            }
        }

    return 0;
}

NManifold* NSnappedBall::manifold() const {
    return new NHandlebody(0, true);
}

NAbelianGroup* NSnappedBall::homology() const {
    return new NAbelianGroup();
}

} // namespace regina

