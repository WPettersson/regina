
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

#include <algorithm>
#include "manifold/nsfs.h"
#include "triangulation/ntriangulation.h"
#include "subcomplex/nplugtrisolidtorus.h"

namespace regina {

const int NPlugTriSolidTorus::CHAIN_NONE = 0;
const int NPlugTriSolidTorus::CHAIN_MAJOR = 1;
const int NPlugTriSolidTorus::CHAIN_MINOR = 3;
const int NPlugTriSolidTorus::EQUATOR_MAJOR = 1;
const int NPlugTriSolidTorus::EQUATOR_MINOR = 3;

NPlugTriSolidTorus::~NPlugTriSolidTorus() {
    if (core)
        delete core;
    for (int i = 0; i < 3; i++)
        if (chain[i])
            delete chain[i];
}

NPlugTriSolidTorus* NPlugTriSolidTorus::clone() const {
    NPlugTriSolidTorus* ans = new NPlugTriSolidTorus();
    ans->core = core->clone();
    for (int i = 0; i < 3; i++) {
        if (chain[i])
            ans->chain[i] = new NLayeredChain(*chain[i]);
        ans->chainType[i] = chainType[i];
    }
    ans->equatorType = equatorType;
    return ans;
}

std::ostream& NPlugTriSolidTorus::writeName(std::ostream& out) const {
    long params[3];
    int nParams = 0;

    int i;
    for (i = 0; i < 3; i++)
        if (chainType[i] != CHAIN_NONE) {
            if (chainType[i] == CHAIN_MAJOR)
                params[nParams++] = chain[i]->getIndex();
            else
                params[nParams++] = -chain[i]->getIndex();
        }
    std::sort(params, params + nParams);

    out << (equatorType == EQUATOR_MAJOR ? "P(" : "P'(");
    if (nParams == 0)
        return out << "0)";
    for (i = 0; i < nParams; i++) {
        if (i > 0)
            out << ',';
        out << params[i];
    }
    return out << ')';
}

std::ostream& NPlugTriSolidTorus::writeTeXName(std::ostream& out) const {
    long params[3];
    int nParams = 0;

    int i;
    for (i = 0; i < 3; i++)
        if (chainType[i] != CHAIN_NONE) {
            if (chainType[i] == CHAIN_MAJOR)
                params[nParams++] = chain[i]->getIndex();
            else
                params[nParams++] = -chain[i]->getIndex();
        }
    std::sort(params, params + nParams);

    out << (equatorType == EQUATOR_MAJOR ? "P_{" : "P'_{");
    if (nParams == 0)
        return out << "0}";
    for (i = 0; i < nParams; i++) {
        if (i > 0)
            out << ',';
        out << params[i];
    }
    return out << '}';
}

void NPlugTriSolidTorus::writeTextLong(std::ostream& out) const {
    out << "Plugged triangular solid torus: ";
    writeName(out);
}

NManifold* NPlugTriSolidTorus::getManifold() const {
    NSFSpace* ans = new NSFSpace();
    ans->insertFibre(2, -1);
    ans->insertFibre(3, 1);

    long rot = (equatorType == EQUATOR_MAJOR ? 5 : 4);
    for (int i = 0; i < 3; i++)
        if (chainType[i] != CHAIN_NONE) {
            if (chainType[i] == equatorType)
                rot += chain[i]->getIndex();
            else
                rot -= chain[i]->getIndex();
        }
    if (rot != 0)
        ans->insertFibre(rot, 1);
    else {
        delete ans;
        return 0;
    }

    ans->reduce();
    return ans;
}

NPlugTriSolidTorus* NPlugTriSolidTorus::isPlugTriSolidTorus(
        NComponent* comp) {
    // TODO: Each triangular solid torus is tested three times since we
    // can't call getTetrahedronIndex() from within a component only.

    // Basic property checks.
    if ((! comp->isClosed()) || (! comp->isOrientable()))
        return 0;

    if (comp->getNumberOfVertices() > 1)
        return 0;

    unsigned long nTet = comp->getNumberOfTetrahedra();
    if (nTet < 5)
        return 0;

    // We have a 1-vertex closed orientable component with at least
    // 5 tetrahedra.

    // Hunt for a core.  Make sure we find each triangular solid torus
    // just once.
    unsigned long tetIndex;
    int coreIndex;
    NTriSolidTorus* core;
    NTetrahedron* coreTet[3];
    NEdge* axis[3];
    NPerm4 coreRoles[3];
    NTetrahedron* base[2];
    NPerm4 baseRoles[2];
    int i, j;
    bool error;

    NTetrahedron* plugTet[3][2];
    NPerm4 plugRoles[3][2];
    NPerm4 realPlugRoles[2];

    NLayeredChain* chain[3];
    int chainType[3];
    int equatorType = 0;

    chain[0] = chain[1] = chain[2] = 0;

    for (tetIndex = 0; tetIndex < nTet - 2; tetIndex++)
        for (coreIndex = 0; coreIndex < 24; coreIndex++) {
            coreRoles[0] = NPerm4::S4[coreIndex];
            if (coreRoles[0][0] > coreRoles[0][3])
                continue;

            core = NTriSolidTorus::formsTriSolidTorus(
                comp->getTetrahedron(tetIndex), coreRoles[0]);
            if (! core)
                continue;

            for (i = 0; i < 3; i++) {
                coreTet[i] = core->getTetrahedron(i);
                coreRoles[i] = core->getVertexRoles(i);
                axis[i] = coreTet[i]->getEdge(
                    NEdge::edgeNumber[coreRoles[i][0]][coreRoles[i][3]]);
            }

            if (axis[0] == axis[1] || axis[1] == axis[2] ||
                    axis[2] == axis[0]) {
                delete core;
                continue;
            }

            // We have the triangular solid torus and we know the three
            // axis edges are distinct.

            // Hunt for chains.
            for (i = 0; i < 3; i++) {
                base[0] = coreTet[(i + 1) % 3]->adjacentTetrahedron(
                    coreRoles[(i + 1) % 3][2]);
                base[1] = coreTet[(i + 2) % 3]->adjacentTetrahedron(
                    coreRoles[(i + 2) % 3][1]);
                if (base[0] != base[1]) {
                    // No chain.
                    chainType[i] = CHAIN_NONE;
                    continue;
                }

                // Have we layered over the major axis?
                baseRoles[0] = coreTet[(i + 1) % 3]->
                    adjacentGluing(coreRoles[(i + 1) % 3][2]) *
                    coreRoles[(i + 1) % 3] * NPerm4(0, 3, 2, 1);
                baseRoles[1] = coreTet[(i + 2) % 3]->
                    adjacentGluing(coreRoles[(i + 2) % 3][1]) *
                    coreRoles[(i + 2) % 3] * NPerm4(2, 1, 0, 3);
                if (baseRoles[0] == baseRoles[1]) {
                    chainType[i] = CHAIN_MAJOR;
                    chain[i] = new NLayeredChain(base[0], baseRoles[0]);
                    while (chain[i]->extendAbove())
                        ;
                    continue;
                }

                // Have we layered over the minor axis?
                baseRoles[0] = coreTet[(i + 1) % 3]->
                    adjacentGluing(coreRoles[(i + 1) % 3][2]) *
                    coreRoles[(i + 1) % 3] * NPerm4(3, 0, 2, 1);
                baseRoles[1] = coreTet[(i + 2) % 3]->
                    adjacentGluing(coreRoles[(i + 2) % 3][1]) *
                    coreRoles[(i + 2) % 3] * NPerm4(2, 1, 3, 0);
                if (baseRoles[0] == baseRoles[1]) {
                    chainType[i] = CHAIN_MINOR;
                    chain[i] = new NLayeredChain(base[0], baseRoles[0]);
                    while (chain[i]->extendAbove())
                        ;
                    continue;
                }

                // It's not a chain but it can't be a plug either.
                // We'll notice the error because i will be less than 3.
                break;
            }

            // Check whether we broke out of the previous loop with an error.
            // Check also whether one of the chains is another in
            // reverse, and that we've found the correct number of
            // tetrahedra in total.
            error = false;
            if (i < 3)
                error = true;
            else if (chain[0] && chain[1] &&
                    chain[0]->getBottom() == chain[1]->getTop())
                error = true;
            else if (chain[1] && chain[2] &&
                    chain[1]->getBottom() == chain[2]->getTop())
                error = true;
            else if (chain[2] && chain[0] &&
                    chain[2]->getBottom() == chain[0]->getTop())
                error = true;
            else if ((chain[0] ? chain[0]->getIndex() : 0) +
                    (chain[1] ? chain[1]->getIndex() : 0) +
                    (chain[2] ? chain[2]->getIndex() : 0) +
                    5 != nTet)
                error = true;

            if (error) {
                for (j = 0; j < 3; j++)
                    if (chain[j]) {
                        delete chain[j];
                        chain[j] = 0;
                    }
                delete core;
                continue;
            }

            // Still hanging in.
            // We know there's only 2 tetrahedra left.
            // Now we need to check the plug.
            error = false;

            for (i = 0; i < 3; i++) {
                if (chain[i]) {
                    plugTet[i][0] = chain[i]->getTop()->adjacentTetrahedron(
                        chain[i]->getTopVertexRoles()[3]);
                    plugTet[i][1] = chain[i]->getTop()->adjacentTetrahedron(
                        chain[i]->getTopVertexRoles()[0]);
                    plugRoles[i][0] = chain[i]->getTop()->
                        adjacentGluing(chain[i]->
                        getTopVertexRoles()[3]) *
                        chain[i]->getTopVertexRoles() *
                        (chainType[i] == CHAIN_MAJOR ? NPerm4(0, 1, 2, 3) :
                        NPerm4(1, 0, 2, 3));
                    plugRoles[i][1] = chain[i]->getTop()->
                        adjacentGluing(chain[i]->
                        getTopVertexRoles()[0]) *
                        chain[i]->getTopVertexRoles() *
                        (chainType[i] == CHAIN_MAJOR ? NPerm4(2, 3, 1, 0) :
                        NPerm4(3, 2, 1, 0));
                } else {
                    plugTet[i][0] = coreTet[(i + 1) % 3]->
                        adjacentTetrahedron(coreRoles[(i + 1) % 3][2]);
                    plugTet[i][1] = coreTet[(i + 2) % 3]->
                        adjacentTetrahedron(coreRoles[(i + 2) % 3][1]);
                    plugRoles[i][0] = coreTet[(i + 1) % 3]->
                        adjacentGluing(coreRoles[(i + 1) % 3][2])
                        * coreRoles[(i + 1) % 3] * NPerm4(0, 3, 1, 2);
                    plugRoles[i][1] = coreTet[(i + 2) % 3]->
                        adjacentGluing(coreRoles[(i + 2) % 3][1])
                        * coreRoles[(i + 2) % 3] * NPerm4(0, 3, 2, 1);
                }
            }

            // Make sure we meet precisely two tetrahedra, three times
            // each.  Note that this implies that the plug tetrahedra are
            // in fact thus far unseen.
            for (i = 0; i < 2; i++)
                if (plugTet[0][i] != plugTet[1][i] ||
                        plugTet[1][i] != plugTet[2][i]) {
                    error = true;
                    break;
                }

            // Make sure also that the gluing permutations for the plug
            // are correct.
            if (! error) {
                if (plugRoles[0][0][0] == plugRoles[1][0][0] &&
                        plugRoles[1][0][0] == plugRoles[2][0][0]) {
                    // Type EQUATOR_MINOR.
                    realPlugRoles[0] = plugRoles[0][0] * NPerm4(3, 2, 1, 0);
                    realPlugRoles[1] = plugRoles[0][1] * NPerm4(3, 0, 2, 1);

                    if (realPlugRoles[0] != plugRoles[1][0] *
                            NPerm4(1, 3, 2, 0))
                        error = true;
                    else if (realPlugRoles[0] != plugRoles[2][0] *
                            NPerm4(2, 1, 3, 0))
                        error = true;
                    else if (realPlugRoles[1] != plugRoles[1][1] *
                            NPerm4(2, 3, 0, 1))
                        error = true;
                    else if (realPlugRoles[1] != plugRoles[2][1] *
                            NPerm4(0, 2, 3, 1))
                        error = true;
                    else
                        equatorType = EQUATOR_MINOR;
                } else if (plugRoles[0][0][1] == plugRoles[1][0][1] &&
                        plugRoles[1][0][1] == plugRoles[2][0][1]) {
                    // Type EQUATOR_MAJOR.
                    realPlugRoles[0] = plugRoles[0][0] * NPerm4(3, 2, 0, 1);
                    realPlugRoles[1] = plugRoles[0][1] * NPerm4(3, 1, 2, 0);

                    if (realPlugRoles[0] != plugRoles[1][0] *
                            NPerm4(0, 3, 2, 1))
                        error = true;
                    else if (realPlugRoles[0] != plugRoles[2][0] *
                            NPerm4(2, 0, 3, 1))
                        error = true;
                    else if (realPlugRoles[1] != plugRoles[1][1] *
                            NPerm4(2, 3, 1, 0))
                        error = true;
                    else if (realPlugRoles[1] != plugRoles[2][1] *
                            NPerm4(1, 2, 3, 0))
                        error = true;
                    else
                        equatorType = EQUATOR_MAJOR;
                } else
                    error = true;
            }

            // Finally check the internal face of the plug.
            if (! error) {
                if (plugTet[0][0]->adjacentTetrahedron(realPlugRoles[0][3])
                        != plugTet[0][1])
                    error = true;
                else if (plugTet[0][0]->adjacentGluing(
                        realPlugRoles[0][3]) * realPlugRoles[0] !=
                        realPlugRoles[1])
                    error = true;
            }

            if (error) {
                for (j = 0; j < 3; j++)
                    if (chain[j]) {
                        delete chain[j];
                        chain[j] = 0;
                    }
                delete core;
                continue;
            }

            // Success!
            NPlugTriSolidTorus* plug = new NPlugTriSolidTorus();
            plug->core = core;
            for (i = 0; i < 3; i++) {
                plug->chain[i] = chain[i];
                plug->chainType[i] = chainType[i];
            }
            plug->equatorType = equatorType;
            return plug;
        }

    // Nothing was found.
    return 0;
}

} // namespace regina
