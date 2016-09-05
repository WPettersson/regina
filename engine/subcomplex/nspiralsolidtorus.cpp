
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

#include <set>
#include <vector>
#include "algebra/nabeliangroup.h"
#include "manifold/nhandlebody.h"
#include "subcomplex/nspiralsolidtorus.h"
#include "triangulation/ntriangulation.h"

namespace regina {

NSpiralSolidTorus* NSpiralSolidTorus::clone() const {
    NSpiralSolidTorus* ans = new NSpiralSolidTorus(nTet);
    for (size_t i = 0; i < nTet; i++) {
        ans->tet[i] = tet[i];
        ans->vertexRoles_[i] = vertexRoles_[i];
    }
    return ans;
}

void NSpiralSolidTorus::reverse() {
    NTetrahedron** newTet = new NTetrahedron*[nTet];
    NPerm4* newRoles = new NPerm4[nTet];

    NPerm4 switchPerm(3, 2, 1, 0);
    for (size_t i = 0; i < nTet; i++) {
        newTet[i] = tet[nTet - 1 - i];
        newRoles[i] = vertexRoles_[nTet - 1 - i] * switchPerm;
    }

    delete[] tet;
    delete[] vertexRoles_;
    tet = newTet;
    vertexRoles_ = newRoles;
}

void NSpiralSolidTorus::cycle(size_t k) {
    NTetrahedron** newTet = new NTetrahedron*[nTet];
    NPerm4* newRoles = new NPerm4[nTet];

    for (size_t i = 0; i < nTet; i++) {
        newTet[i] = tet[(i + k) % nTet];
        newRoles[i] = vertexRoles_[(i + k) % nTet];
    }

    delete[] tet;
    delete[] vertexRoles_;
    tet = newTet;
    vertexRoles_ = newRoles;
}

bool NSpiralSolidTorus::makeCanonical(const NTriangulation* tri) {
    size_t i, index;

    size_t baseTet = 0;
    size_t baseIndex = tet[0]->index();
    for (i = 1; i < nTet; i++) {
        index = tet[i]->index();
        if (index < baseIndex) {
            baseIndex = index;
            baseTet = i;
        }
    }

    bool reverseAlso = (vertexRoles_[baseTet][0] > vertexRoles_[baseTet][3]);

    if (baseTet == 0 && (! reverseAlso))
        return false;

    NTetrahedron** newTet = new NTetrahedron*[nTet];
    NPerm4* newRoles = new NPerm4[nTet];

    if (reverseAlso) {
        // Make baseTet into tetrahedron 0 and reverse.
        NPerm4 switchPerm(3, 2, 1, 0);
        for (size_t i = 0; i < nTet; i++) {
            newTet[i] = tet[(baseTet + nTet - i) % nTet];
            newRoles[i] = vertexRoles_[(baseTet + nTet - i) % nTet] *
                switchPerm;
        }
    } else {
        // Make baseTet into tetrahedron 0 but don't reverse.
        for (size_t i = 0; i < nTet; i++) {
            newTet[i] = tet[(i + baseTet) % nTet];
            newRoles[i] = vertexRoles_[(i + baseTet) % nTet];
        }
    }

    delete[] tet;
    delete[] vertexRoles_;
    tet = newTet;
    vertexRoles_ = newRoles;

    return true;
}

bool NSpiralSolidTorus::isCanonical(const NTriangulation* tri) const {
    if (vertexRoles_[0][0] > vertexRoles_[0][3])
        return false;

    long baseIndex = tet[0]->index();
    for (size_t i = 1; i < nTet; i++)
        if (static_cast<int>(tet[i]->index()) < baseIndex)
            return false;

    return true;
}

NSpiralSolidTorus* NSpiralSolidTorus::formsSpiralSolidTorus(NTetrahedron* tet,
        NPerm4 useVertexRoles) {
    NPerm4 invRoleMap(1, 2, 3, 0);  // Maps upper roles to lower roles.

    NTetrahedron* base = tet;
    NPerm4 baseRoles(useVertexRoles);

    std::vector<NTetrahedron*> tets;
    std::vector<NPerm4> roles;
    std::set<NTetrahedron*> usedTets;

    tets.push_back(tet);
    roles.push_back(useVertexRoles);
    usedTets.insert(tet);

    NTetrahedron* adjTet;
    NPerm4 adjRoles;

    while (1) {
        // Examine the tetrahedron beyond tet.
        adjTet = tet->adjacentTetrahedron(useVertexRoles[0]);
        adjRoles = tet->adjacentGluing(useVertexRoles[0]) *
            useVertexRoles * invRoleMap;

        // Check that we haven't hit the boundary.
        if (adjTet == 0)
            return 0;

        if (adjTet == base) {
            // We're back at the beginning of the loop.
            // Check that everything is glued up correctly.
            if (adjRoles != baseRoles)
                return 0;

            // Success!
            break;
        }

        if (usedTets.count(adjTet))
            return 0;

        // Move on to the next tetrahedron.
        tet = adjTet;
        useVertexRoles = adjRoles;

        tets.push_back(tet);
        roles.push_back(useVertexRoles);
        usedTets.insert(tet);
    }

    // We've found a spiralled solid torus.
    NSpiralSolidTorus* ans = new NSpiralSolidTorus(tets.size());
    copy(tets.begin(), tets.end(), ans->tet);
    copy(roles.begin(), roles.end(), ans->vertexRoles_);
    return ans;
}

NManifold* NSpiralSolidTorus::manifold() const {
    return new NHandlebody(1, true);
}

NAbelianGroup* NSpiralSolidTorus::homology() const {
    NAbelianGroup* ans = new NAbelianGroup();
    ans->addRank();
    return ans;
}

} // namespace regina

