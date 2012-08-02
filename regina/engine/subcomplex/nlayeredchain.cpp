
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

#include "algebra/nabeliangroup.h"
#include "manifold/nhandlebody.h"
#include "triangulation/nedge.h"
#include "triangulation/ntetrahedron.h"
#include "subcomplex/nlayeredchain.h"

namespace regina {

bool NLayeredChain::extendAbove() {
    NTetrahedron* adj = top->adjacentTetrahedron(topVertexRoles[0]);
    if (adj == bottom || adj == top || adj == 0)
        return false;
    if (adj != top->adjacentTetrahedron(topVertexRoles[3]))
        return false;

    // Check the gluings.
    NPerm4 adjRoles = top->adjacentGluing(topVertexRoles[0]) *
        topVertexRoles * NPerm4(0, 1);
    if (adjRoles != top->adjacentGluing(topVertexRoles[3]) *
            topVertexRoles * NPerm4(2, 3))
        return false;

    // We can extend the layered chain.
    top = adj;
    topVertexRoles = adjRoles;
    index++;
    return true;
}

bool NLayeredChain::extendBelow() {
    NTetrahedron* adj = bottom->adjacentTetrahedron(bottomVertexRoles[1]);
    if (adj == bottom || adj == top || adj == 0)
        return false;
    if (adj != bottom->adjacentTetrahedron(bottomVertexRoles[2]))
        return false;

    // Check the gluings.
    NPerm4 adjRoles = bottom->adjacentGluing(bottomVertexRoles[1])
        * bottomVertexRoles * NPerm4(0, 1);
    if (adjRoles != bottom->adjacentGluing(bottomVertexRoles[2])
            * bottomVertexRoles * NPerm4(2, 3))
        return false;

    // We can extend the layered chain.
    bottom = adj;
    bottomVertexRoles = adjRoles;
    index++;
    return true;
}

bool NLayeredChain::extendMaximal() {
    bool changed = false;
    while (extendAbove())
        changed = true;
    while (extendBelow())
        changed = true;
    return changed;
}

void NLayeredChain::reverse() {
    NTetrahedron* tmp = top;
    top = bottom;
    bottom = tmp;

    NPerm4 pTmp = topVertexRoles * NPerm4(1, 0, 3, 2);
    topVertexRoles = bottomVertexRoles * NPerm4(1, 0, 3, 2);
    bottomVertexRoles = pTmp;
}

void NLayeredChain::invert() {
    topVertexRoles = topVertexRoles * NPerm4(3, 2, 1, 0);
    bottomVertexRoles = bottomVertexRoles * NPerm4(3, 2, 1, 0);
}

NManifold* NLayeredChain::getManifold() const {
    return new NHandlebody(index <= 1 ? 0 : 1, true);
}

NAbelianGroup* NLayeredChain::getHomologyH1() const {
    NAbelianGroup* ans = new NAbelianGroup();
    if (index > 1)
        ans->addRank();
    return ans;
}

} // namespace regina

