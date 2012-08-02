
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

#include "split/nsignature.h"
#include "triangulation/nexampletriangulation.h"
#include "triangulation/ntriangulation.h"

namespace {
    static const int poincareAdj[5][4] = {
        { 1, 2, 3, 4},
        { 0, 2, 3, 4},
        { 0, 1, 3, 4},
        { 0, 1, 2, 4},
        { 0, 1, 2, 3}
    };

    static const int poincareGluings[5][4][4] = {
        { { 0, 2, 3, 1 }, { 3, 0, 1, 2 }, { 2, 3, 0, 1 }, { 3, 1, 2, 0 } },
        { { 0, 3, 1, 2 }, { 2, 1, 3, 0 }, { 3, 0, 1, 2 }, { 2, 3, 0, 1 } },
        { { 1, 2, 3, 0 }, { 3, 1, 0, 2 }, { 1, 3, 2, 0 }, { 3, 0, 1, 2 } },
        { { 2, 3, 0, 1 }, { 1, 2, 3, 0 }, { 3, 0, 2, 1 }, { 1, 2, 0, 3 } },
        { { 3, 1, 2, 0 }, { 2, 3, 0, 1 }, { 1, 2, 3, 0 }, { 2, 0, 1, 3 } }
    };

    static const int weeksAdj[9][4] = {
        { 0, 0, 1, 2},
        { 0, 3, 4, 5},
        { 0, 3, 4, 6},
        { 1, 2, 5, 7},
        { 1, 2, 7, 8},
        { 1, 3, 6, 8},
        { 2, 5, 7, 8},
        { 3, 4, 6, 8},
        { 4, 5, 6, 7}
    };

    static const int weeksGluings[9][4][4] = {
        { { 1, 2, 3, 0 }, { 3, 0, 1, 2 }, { 3, 2, 0, 1 }, { 2, 3, 1, 0 } },
        { { 2, 3, 1, 0 }, { 1, 0, 2, 3 }, { 1, 3, 0, 2 }, { 2, 3, 1, 0 } },
        { { 3, 2, 0, 1 }, { 0, 1, 3, 2 }, { 0, 2, 1, 3 }, { 1, 3, 2, 0 } },
        { { 1, 0, 2, 3 }, { 0, 1, 3, 2 }, { 2, 3, 1, 0 }, { 3, 2, 1, 0 } },
        { { 2, 0, 3, 1 }, { 0, 2, 1, 3 }, { 0, 3, 1, 2 }, { 2, 3, 1, 0 } },
        { { 3, 2, 0, 1 }, { 3, 2, 0, 1 }, { 0, 3, 1, 2 }, { 3, 2, 0, 1 } },
        { { 3, 0, 2, 1 }, { 0, 2, 3, 1 }, { 3, 1, 2, 0 }, { 1, 0, 3, 2 } },
        { { 3, 2, 1, 0 }, { 0, 2, 3, 1 }, { 3, 1, 2, 0 }, { 1, 2, 0, 3 } },
        { { 3, 2, 0, 1 }, { 2, 3, 1, 0 }, { 1, 0, 3, 2 }, { 2, 0, 1, 3 } }
    };

    static const int closedOrHypAdj[9][4] = {
        { 6, 8, 2, 8 },
        { 6, 8, 3, 7 },
        { 7, 0, 3, 4 },
        { 1, 5, 5, 2 },
        { 2, 6, 5, 7 },
        { 3, 8, 3, 4 },
        { 0, 4, 7, 1 },
        { 1, 4, 2, 6 },
        { 1, 0, 5, 0 }
    };

    static const int closedOrHypGluings[9][4][4] = {
        { { 0, 1, 3, 2 }, { 3, 1, 2, 0 }, { 0, 2, 1, 3 }, { 0, 2, 1, 3 } },
        { { 3, 1, 2, 0 }, { 1, 0, 2, 3 }, { 3, 2, 0, 1 }, { 2, 3, 1, 0 } },
        { { 2, 0, 3, 1 }, { 0, 2, 1, 3 }, { 0, 1, 3, 2 }, { 3, 1, 2, 0 } },
        { { 2, 3, 1, 0 }, { 3, 2, 0, 1 }, { 2, 1, 0, 3 }, { 0, 1, 3, 2 } },
        { { 3, 1, 2, 0 }, { 0, 1, 3, 2 }, { 0, 1, 3, 2 }, { 3, 2, 0, 1 } },
        { { 2, 1, 0, 3 }, { 0, 2, 1, 3 }, { 2, 3, 1, 0 }, { 0, 1, 3, 2 } },
        { { 0, 1, 3, 2 }, { 0, 1, 3, 2 }, { 0, 1, 3, 2 }, { 3, 1, 2, 0 } },
        { { 3, 2, 0, 1 }, { 2, 3, 1, 0 }, { 1, 3, 0, 2 }, { 0, 1, 3, 2 } },
        { { 1, 0, 2, 3 }, { 3, 1, 2, 0 }, { 0, 2, 1, 3 }, { 0, 2, 1, 3 } }
    };

    static const int closedNorHypAdj[11][4] = {
        { 8, 2, 8, 2 },
        { 5, 3, 2, 9 },
        { 1, 4, 0, 0 },
        { 6, 1, 4, 6 },
        { 10, 2, 10, 3 },
        { 7, 7, 6, 1 },
        { 8, 3, 3, 5 },
        { 5, 9, 8, 5 },
        { 0, 0, 6, 7 },
        { 10, 10, 1, 7 },
        { 9, 4, 4, 9 }
    };

    static const int closedNorHypGluings[11][4][4] = {
        { { 1, 3, 2, 0 }, { 0, 3, 2, 1 }, { 2, 1, 0, 3 }, { 3, 1, 0, 2 } },
        { { 3, 0, 1, 2 }, { 3, 1, 0, 2 }, { 2, 1, 0, 3 }, { 1, 0, 3, 2 } },
        { { 2, 1, 0, 3 }, { 3, 1, 2, 0 }, { 2, 1, 3, 0 }, { 0, 3, 2, 1 } },
        { { 2, 1, 3, 0 }, { 2, 1, 3, 0 }, { 2, 0, 3, 1 }, { 0, 3, 2, 1 } },
        { { 2, 1, 0, 3 }, { 3, 1, 2, 0 }, { 3, 2, 1, 0 }, { 1, 3, 0, 2 } },
        { { 3, 1, 2, 0 }, { 1, 0, 3, 2 }, { 0, 1, 3, 2 }, { 1, 2, 3, 0 } },
        { { 2, 1, 0, 3 }, { 0, 3, 2, 1 }, { 3, 1, 0, 2 }, { 0, 1, 3, 2 } },
        { { 1, 0, 3, 2 }, { 0, 3, 2, 1 }, { 0, 1, 3, 2 }, { 3, 1, 2, 0 } },
        { { 2, 1, 0, 3 }, { 3, 0, 2, 1 }, { 2, 1, 0, 3 }, { 0, 1, 3, 2 } },
        { { 3, 1, 2, 0 }, { 2, 0, 1, 3 }, { 1, 0, 3, 2 }, { 0, 3, 2, 1 } },
        { { 1, 2, 0, 3 }, { 3, 2, 1, 0 }, { 2, 1, 0, 3 }, { 3, 1, 2, 0 } }
    };

    static const int whiteheadAdj[4][4] = {
        { 3, 2, 1, 3},
        { 3, 2, 2, 0},
        { 1, 3, 0, 1},
        { 2, 0, 0, 1}
    };

    static const int whiteheadGluings[4][4][4] = {
        { { 2, 3, 1, 0 }, { 3, 2, 0, 1 }, { 0, 1, 3, 2 }, { 3, 2, 0, 1 } },
        { { 3, 2, 0, 1 }, { 2, 3, 1, 0 }, { 3, 2, 0, 1 }, { 0, 1, 3, 2 } },
        { { 2, 3, 1, 0 }, { 1, 0, 2, 3 }, { 2, 3, 1, 0 }, { 3, 2, 0, 1 } },
        { { 1, 0, 2, 3 }, { 2, 3, 1, 0 }, { 3, 2, 0, 1 }, { 2, 3, 1, 0 } }
    };
}

namespace regina {

NTriangulation* NExampleTriangulation::threeSphere() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("3-sphere");

    ans->insertLayeredLensSpace(1, 0);

    return ans;
}

NTriangulation* NExampleTriangulation::bingsHouse() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Bing's house with two rooms");

    NTetrahedron* r = ans->newTetrahedron();
    NTetrahedron* s = ans->newTetrahedron();
    r->joinTo(0, r, NPerm4(0, 1));
    s->joinTo(0, s, NPerm4(0, 1));
    r->joinTo(3, s, NPerm4(3, 1, 0, 2));
    s->joinTo(3, r, NPerm4(3, 1, 0, 2));

    return ans;
}

NTriangulation* NExampleTriangulation::s2xs1() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("S2 x S1");

    ans->insertLayeredLensSpace(0, 1);

    return ans;
}

NTriangulation* NExampleTriangulation::rp2xs1() {
    // Section 3.5.1 of Benjamin Burton's PhD thesis describes how to
    // construct RP^2 x S^1 by identifying the boundary faces of a
    // solid Klein bottle.
    NTriangulation* ans = solidKleinBottle();
    ans->setPacketLabel("RP2 x S1");

    NTetrahedron* r = ans->getTetrahedron(0);
    NTetrahedron* t = ans->getTetrahedron(2);
    r->joinTo(1, t, NPerm4(2, 3, 0, 1));
    r->joinTo(3, t, NPerm4(2, 3, 0, 1));

    return ans;
}

NTriangulation* NExampleTriangulation::rp3rp3() {
    // This can be generated as the enclosing triangulation of a splitting
    // surface, as described in chapter 4 of Benjamin Burton's PhD thesis.
    NSignature* sig = NSignature::parse("aabccd.b.d");
    NTriangulation* ans = sig->triangulate();
    ans->setPacketLabel("RP3 # RP3");
    delete sig;

    return ans;
}

NTriangulation* NExampleTriangulation::lens8_3() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("L(8,3)");

    ans->insertLayeredLensSpace(8, 3);

    return ans;
}

NTriangulation* NExampleTriangulation::poincareHomologySphere() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Poincare homology sphere");

    ans->insertConstruction(5, poincareAdj, poincareGluings);

    return ans;
}

NTriangulation* NExampleTriangulation::weeks() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Weeks manifold");

    ans->insertConstruction(9, weeksAdj, weeksGluings);

    return ans;
}

NTriangulation* NExampleTriangulation::weberSeifert() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Weber-Seifert dodecahedral space");

    // Bah.  Dehydration strings are somewhat impenetrable,
    // but the alternative is 23 lines of hard-coded tetrahedron gluings.
    //
    // This triangulation was constructed by building a 60-tetrahedron
    // dodecahedron and identifying opposite faces with a 3/10 twist,
    // and then simplifying down to one vertex and 23 tetrahedra.
    ans->insertRehydration(
        "xppphocgaeaaahimmnkontspmuuqrsvuwtvwwxwjjsvvcxxjjqattdwworrko");

    return ans;
}

NTriangulation* NExampleTriangulation::seifertWeber() {
    // Kept for backward compatibility.
    // Use the old name in the packet label.
    NTriangulation* ans = weberSeifert();
    ans->setPacketLabel("Seifert-Weber dodecahedral space");
    return ans;
}

NTriangulation* NExampleTriangulation::smallClosedOrblHyperbolic() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Closed orientable hyperbolic 3-manifold");

    ans->insertConstruction(9, closedOrHypAdj, closedOrHypGluings);

    return ans;
}

NTriangulation* NExampleTriangulation::smallClosedNonOrblHyperbolic() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Closed non-orientable hyperbolic 3-manifold");

    ans->insertConstruction(11, closedNorHypAdj, closedNorHypGluings);

    return ans;
}

NTriangulation* NExampleTriangulation::lst3_4_7() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Layered solid torus");

    ans->insertLayeredSolidTorus(3, 4);

    return ans;
}

NTriangulation* NExampleTriangulation::solidKleinBottle() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Solid Klein bottle");

    // A three-tetrahedron solid Klein bottle is described in section
    // 3.5.1 of Benjamin Burton's PhD thesis.
    NTetrahedron* r = ans->newTetrahedron();
    NTetrahedron* s = ans->newTetrahedron();
    NTetrahedron* t = ans->newTetrahedron();
    s->joinTo(0, r, NPerm4(0, 1, 2, 3));
    s->joinTo(3, r, NPerm4(3, 0, 1, 2));
    s->joinTo(1, t, NPerm4(3, 0, 1, 2));
    s->joinTo(2, t, NPerm4(0, 1, 2, 3));

    return ans;
}

NTriangulation* NExampleTriangulation::figureEightKnotComplement() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Figure eight knot complement");

    // The two-tetrahedron figure eight knot complement is described at
    // the beginning of chapter 8 of Richard Rannard's PhD thesis.
    NTetrahedron* r = ans->newTetrahedron();
    NTetrahedron* s = ans->newTetrahedron();
    r->joinTo(0, s, NPerm4(1, 3, 0, 2));
    r->joinTo(1, s, NPerm4(2, 0, 3, 1));
    r->joinTo(2, s, NPerm4(0, 3, 2, 1));
    r->joinTo(3, s, NPerm4(2, 1, 0, 3));

    return ans;
}

NTriangulation* NExampleTriangulation::whiteheadLinkComplement() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Whitehead link complement");

    ans->insertConstruction(4, whiteheadAdj, whiteheadGluings);

    return ans;
}

NTriangulation* NExampleTriangulation::gieseking() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Gieseking manifold");

    NTetrahedron* r = ans->newTetrahedron();
    r->joinTo(0, r, NPerm4(1, 2, 0, 3));
    r->joinTo(2, r, NPerm4(0, 2, 3, 1));

    return ans;
}

NTriangulation* NExampleTriangulation::cuspedGenusTwoTorus() {
    NTriangulation* ans = new NTriangulation();
    ans->setPacketLabel("Cusped genus two solid torus");

    // We create this by first constructing an ordinary solid genus two
    // torus and then converting the real boundary to an ideal vertex.
    NTetrahedron* r = ans->newTetrahedron();
    NTetrahedron* s = ans->newTetrahedron();
    NTetrahedron* t = ans->newTetrahedron();
    NTetrahedron* u = ans->newTetrahedron();
    r->joinTo(0, s, NPerm4());
    r->joinTo(1, t, NPerm4(1, 2, 3, 0));
    r->joinTo(2, u, NPerm4(1, 0, 3, 2));
    s->joinTo(3, t, NPerm4());
    t->joinTo(1, u, NPerm4());
    ans->finiteToIdeal();

    return ans;
}

} // namespace regina

