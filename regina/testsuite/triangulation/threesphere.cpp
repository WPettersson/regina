
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Test Suite                                                            *
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

#include <cmath>
#include <memory>
#include <sstream>
#include <cppunit/extensions/HelperMacros.h>
#include "census/ncensus.h"
#include "packet/ncontainer.h"
#include "split/nsignature.h"
#include "triangulation/nexampletriangulation.h"
#include "triangulation/ntriangulation.h"
#include "testsuite/triangulation/testtriangulation.h"
#include "utilities/nbooleans.h"

using regina::NBoolSet;
using regina::NCensus;
using regina::NContainer;
using regina::NExampleTriangulation;
using regina::NPerm4;
using regina::NSignature;
using regina::NTetrahedron;
using regina::NTriangulation;

class ThreeSphereTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(ThreeSphereTest);

    CPPUNIT_TEST(threeSphereRecognition);
    CPPUNIT_TEST(threeSphereCensus);
    CPPUNIT_TEST(threeBallRecognition);
    CPPUNIT_TEST(solidTorusRecognition);

    CPPUNIT_TEST_SUITE_END();

    public:
        NTriangulation* generateFromSig(const std::string& sigStr) {
            NSignature* sig = NSignature::parse(sigStr);
            if (sig == 0)
                return 0;

            NTriangulation* triNew = sig->triangulate();
            delete sig;
            return triNew;
        }

        void setUp() {
        }

        void tearDown() {
        }

        NTriangulation* verifyThreeSphere(NTriangulation* tri,
                const std::string& triName) {
            if (! tri->isThreeSphere()) {
                CPPUNIT_FAIL(("The 3-sphere " + triName +
                    " is not recognised as such.").c_str());
            }

            // Try again with a barycentric subdivision.
            NTriangulation big(*tri);
            big.barycentricSubdivision();
            if (! big.isThreeSphere()) {
                CPPUNIT_FAIL(("The barycentric subdivision of the 3-sphere "
                    + triName + " is not recognised as such.").c_str());
            }

            return tri;
        }

        NTriangulation* verifyNotThreeSphere(NTriangulation* tri,
                const std::string& triName) {
            if (tri->isThreeSphere()) {
                CPPUNIT_FAIL(("The non-3-sphere " + triName +
                    " is recognised as a 3-sphere.").c_str());
            }

            // Try again with a barycentric subdivision.
            NTriangulation big(*tri);
            big.barycentricSubdivision();
            if (big.isThreeSphere()) {
                CPPUNIT_FAIL(("The barycentric subdivision of the non-3-sphere "
                    + triName + " is recognised as a 3-sphere.").c_str());
            }

            return tri;
        }

        void verifySigThreeSphere(const std::string& sigStr) {
            delete verifyThreeSphere(generateFromSig(sigStr), sigStr);
        }

        void verifySigNotThreeSphere(const std::string& sigStr) {
            delete verifyNotThreeSphere(generateFromSig(sigStr), sigStr);
        }

        void threeSphereRecognition() {
            NTriangulation* tri;

            // 3-spheres obtained from splitting surface signatures:
            verifySigThreeSphere("(a)(a)");
            verifySigThreeSphere("(ab)(a)(b)");
            verifySigThreeSphere("(ab)(a)(b)");
            verifySigThreeSphere("(abc)(a)(b)(c)");
            verifySigThreeSphere("(ab)(ac)(b)(c)");
            verifySigThreeSphere("(abcd)(a)(b)(c)(d)");
            verifySigThreeSphere("(abc)(abd)(c)(d)");
            verifySigThreeSphere("(abc)(acd)(b)(d)");
            verifySigThreeSphere("(abc)(ad)(b)(c)(d)");
            verifySigThreeSphere("(ab)(ac)(bd)(cd)");
            verifySigThreeSphere("(ab)(ac)(bd)(c)(d)");
            verifySigThreeSphere("(abcd)(aefg)(b)(c)(d)(e)(f)(g)");

            // 3-spheres obtained as Lens spaces:
            tri = new NTriangulation();
            tri->insertLayeredLensSpace(1,0);
            delete verifyThreeSphere(tri, "L(1,0)");

            // Non-3-spheres obtained from splitting surface signatures:
            verifySigNotThreeSphere("(aab)(b)");
            verifySigNotThreeSphere("(ab)(ab)");
            verifySigNotThreeSphere("(aabcb)(c)");
            verifySigNotThreeSphere("(aabc)(b)(c)");
            verifySigNotThreeSphere("(aabcdcb)(d)");
            verifySigNotThreeSphere("(aabcdb)(c)(d)");
            verifySigNotThreeSphere("(aabcd)(b)(c)(d)");
            verifySigNotThreeSphere("(aabc)(bd)(c)(d)");
            verifySigNotThreeSphere("(abac)(bd)(cd)");
            verifySigNotThreeSphere("(abac)(bd)(c)(d)");
            verifySigNotThreeSphere("(abcd)(ac)(bd)");
            verifySigNotThreeSphere("(aab)(bc)(cd)(d)");
            verifySigNotThreeSphere("(abc)(ab)(cd)(d)");
            verifySigNotThreeSphere("(abc)(ad)(bd)(c)");
            verifySigNotThreeSphere("(abac)(b)(c)");
            verifySigNotThreeSphere("(abacdc)(b)(d)");
            verifySigNotThreeSphere("(abcabd)(c)(d)");
            verifySigNotThreeSphere("(abacd)(b)(c)(d)");
            verifySigNotThreeSphere("(aabc)(bd)(cd)");
            verifySigNotThreeSphere("(abacde)(cf)(fg)(b)(d)(e)(g)");
            verifySigNotThreeSphere("(abc)(abc)");
            verifySigNotThreeSphere("(abc)(acb)");

            verifySigNotThreeSphere("(aa)");
            verifySigNotThreeSphere("(aabb)");
            verifySigNotThreeSphere("(abab)");
            verifySigNotThreeSphere("(aabccb)");
            verifySigNotThreeSphere("(abcabc)");
            verifySigNotThreeSphere("(aab)(bcc)");
            verifySigNotThreeSphere("(aab)(bc)(c)");
            verifySigNotThreeSphere("(ab)(ac)(bc)");
            verifySigNotThreeSphere("(aabcddcb)");
            verifySigNotThreeSphere("(abcdabcd)");
            verifySigNotThreeSphere("(aabbc)(cd)(d)");
            verifySigNotThreeSphere("(aabcb)(cdd)");
            verifySigNotThreeSphere("(aabcb)(cd)(d)");
            verifySigNotThreeSphere("(ababc)(cd)(d)");
            verifySigNotThreeSphere("(abac)(bdcd)");
            verifySigNotThreeSphere("(abac)(bcd)(d)");
            verifySigNotThreeSphere("(abac)(bdd)(c)");
            verifySigNotThreeSphere("(abcd)(abcd)");
            verifySigNotThreeSphere("(abcd)(adcb)");
            verifySigNotThreeSphere("(aab)(bcd)(c)(d)");
            verifySigNotThreeSphere("(abc)(abd)(cd)");
            verifySigNotThreeSphere("(abc)(acd)(bd)");
            verifySigNotThreeSphere("(abcdefgh)(abcdefgh)");

            verifySigNotThreeSphere("(aabccd)(b)(d)");;
            verifySigNotThreeSphere("(abacbd)(cd)");;
            verifySigNotThreeSphere("(aabcdecb)(dfeg)(fg)");
            verifySigNotThreeSphere("(aabbc)(cdef)(egg)(d)(f)");
            verifySigNotThreeSphere("(ababc)(cdef)(egg)(d)(f)");
            verifySigNotThreeSphere("(abcd)(acbe)(dfeg)(f)(g)");
            verifySigNotThreeSphere("(abcd)(adce)(befg)(f)(g)");
            verifySigNotThreeSphere("(abcde)(aff)(cgg)(b)(d)(e)");
            verifySigNotThreeSphere("(abcde)(adf)(bfg)(ce)(g)");
            verifySigNotThreeSphere("(ababc)(cde)(dfg)(fg)(e)");
            verifySigNotThreeSphere("(abcdef)(acegg)(bdf)");
            verifySigNotThreeSphere("(abacde)(dffgg)(b)(c)(e)");
            verifySigNotThreeSphere("(aabccdeffeg)(b)(d)(g)");
            verifySigNotThreeSphere("(aabcde)(cfg)(dgf)(b)(e)");

            // And of course the Poincare homology sphere(S3/P120).
            // We'll build this a few different ways.

            // First, one out of the can:
            tri = NExampleTriangulation::poincareHomologySphere();
            delete verifyNotThreeSphere(tri,
                "Poincare homology sphere (example)");

            // Poincare homology sphere as a plugged triangular solid torus:
            tri = new NTriangulation;
            NTetrahedron* tet[5];
            int i;
            for (i = 0; i < 5; i++)
                tet[i] = tri->newTetrahedron();
            tet[0]->joinTo(0, tet[4], NPerm4(1,0,2,3));
            tet[0]->joinTo(1, tet[3], NPerm4(0,2,3,1));
            tet[0]->joinTo(2, tet[1], NPerm4(0,1,3,2));
            tet[0]->joinTo(3, tet[2], NPerm4(2,1,3,0));
            tet[1]->joinTo(0, tet[3], NPerm4(1,3,2,0));
            tet[1]->joinTo(1, tet[2], NPerm4(0,2,3,1));
            tet[1]->joinTo(2, tet[4], NPerm4(2,1,0,3));
            tet[2]->joinTo(1, tet[4], NPerm4(0,2,3,1));
            tet[2]->joinTo(3, tet[3], NPerm4(3,1,2,0));
            tet[3]->joinTo(3, tet[4], NPerm4(0,1,2,3));
            delete verifyNotThreeSphere(tri,
                "Poincare homology sphere (plugged)");

            // Poincare homology sphere as an augmented triangular solid
            // torus:
            tri = new NTriangulation();
            tri->insertAugTriSolidTorus(2, -1, 3, 1, 5, -4);
            delete verifyNotThreeSphere(tri,
                "Poincare homology sphere (aug I)");

            // Poincare homology sphere as another augmented triangular solid
            // torus:
            tri = new NTriangulation();
            tri->insertAugTriSolidTorus(2, -1, 3, -2, 5, 1);
            delete verifyNotThreeSphere(tri,
                "Poincare homology sphere (aug II)");

            // Let's make sure silly things like balls aren't picked up.
            tri = new NTriangulation();
            tri->newTetrahedron();
            delete verifyNotThreeSphere(tri, "Standalone tetrahedron");

            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[0]->joinTo(0, tet[0], NPerm4(3, 1, 2, 0));
            delete verifyNotThreeSphere(tri, "Snapped tetrahedron");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(1, 2);
            delete verifyNotThreeSphere(tri, "LST(1,2,3)");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(3, 4);
            delete verifyNotThreeSphere(tri, "LST(3,4,7)");

            tri = new NTriangulation();
            delete verifyNotThreeSphere(tri, "Empty triangulation");
        }

        static bool verifyThreeSphereCensus(NTriangulation* tri, void*) {
            // PRECONDITION: tri has <= 6 tetrahedra.
            //
            // We know from the closed orientable census data that the
            // only homology 3-sphere with <= 6 tetrahedra is the
            // Poincare homology 3-sphere.  This can be distinguished
            // from the real 3-sphere using the (5,1) Turaev-Viro
            // invariant -- for S^3 the invariant is 0.138197, and for
            // the Poincare homology sphere it is 0.947214.

            // What did 3-sphere recognition say?
            bool found = tri->isThreeSphere();

            bool real = (tri->isValid() &&
                         tri->isClosed() &&
                         tri->isOrientable() &&
                         tri->isConnected() &&
                         tri->getHomologyH1().isTrivial() &&
                         tri->turaevViro(5, 1) < 0.5);

            if (real && ! found)
                CPPUNIT_FAIL("A census 3-sphere was not recognised "
                    "as a 3-sphere.");
            else if (found && ! real)
                CPPUNIT_FAIL("A census non-3-sphere was recognised "
                    "as a 3-sphere.");

            return false;
        }

        void threeSphereCensus() {
            NContainer parent;

            // Enumerate orientable triangulations.
            NCensus::formCensus(&parent, 4 /* tetrahedra */,
                NBoolSet::sTrue /* finite */,
                NBoolSet::sTrue /* orientable */,
                NBoolSet::sFalse /* bounded */,
                -1, 0, &verifyThreeSphereCensus, 0);

            // Enumerate non-orientable triangulations also, for good measure.
            NCensus::formCensus(&parent, 3 /* tetrahedra */,
                NBoolSet::sTrue /* finite */,
                NBoolSet::sTrue /* orientable */,
                NBoolSet::sFalse /* bounded */,
                -1, 0, &verifyThreeSphereCensus, 0);
        }

        NTriangulation* verifyThreeBall(NTriangulation* tri,
                const std::string& triName) {
            if (! tri->isBall()) {
                CPPUNIT_FAIL(("The 3-ball " + triName +
                    " is not recognised as such.").c_str());
            }

            // Try again with a barycentric subdivision.
            NTriangulation big(*tri);
            big.barycentricSubdivision();
            if (! big.isBall()) {
                CPPUNIT_FAIL(("The barycentric subdivision of the 3-ball "
                    + triName + " is not recognised as such.").c_str());
            }

            return tri;
        }

        NTriangulation* verifyNotThreeBall(NTriangulation* tri,
                const std::string& triName) {
            if (tri->isBall()) {
                CPPUNIT_FAIL(("The non-3-ball" + triName +
                    " is recognised as a 3-ball.").c_str());
            }

            // Try again with a barycentric subdivision.
            NTriangulation big(*tri);
            big.barycentricSubdivision();
            if (big.isBall()) {
                CPPUNIT_FAIL(("The barycentric subdivision of the non-3-ball "
                    + triName + " is recognised as a 3-ball.").c_str());
            }

            return tri;
        }

        void threeBallRecognition() {
            NTriangulation* tri;
            NTetrahedron* tet[4];

            // Balls:
            tri = new NTriangulation();
            tri->newTetrahedron();
            delete verifyThreeBall(tri, "Standalone tetrahedron");

            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[0]->joinTo(0, tet[0], NPerm4(3, 1, 2, 0));
            delete verifyThreeBall(tri, "Snapped tetrahedron");

            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[1] = tri->newTetrahedron();
            tet[0]->joinTo(0, tet[1], NPerm4());
            tet[0]->joinTo(1, tet[1], NPerm4());
            tet[0]->joinTo(2, tet[1], NPerm4());
            delete verifyThreeBall(tri, "Triangular pillow");

            // This ball used to crash the simplification routines once
            // upon a time.  Throw it into the test suite for good measure.
            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[1] = tri->newTetrahedron();
            tet[2] = tri->newTetrahedron();
            tet[3] = tri->newTetrahedron();
            tet[0]->joinTo(2, tet[0], NPerm4(0,2));
            tet[0]->joinTo(1, tet[1], NPerm4(2,0,1,3));
            tet[1]->joinTo(2, tet[2], NPerm4());
            tet[1]->joinTo(1, tet[2], NPerm4(2,0,1,3));
            tet[2]->joinTo(1, tet[3], NPerm4(2,0,1,3));
            tet[3]->joinTo(2, tet[3], NPerm4(1,2));
            delete verifyThreeBall(tri, "4-tetrahedron ball");

            // Non-balls:
            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(1, 2);
            delete verifyNotThreeBall(tri, "LST(1,2,3)");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(3, 4);
            delete verifyNotThreeBall(tri, "LST(3,4,7)");

            tri = new NTriangulation();
            delete verifyNotThreeBall(tri, "Empty triangulation");

            // Make a punctured Poincare homology sphere.
            tri = NExampleTriangulation::poincareHomologySphere();
            tri->barycentricSubdivision();
            tri->removeTetrahedronAt(0);
            tri->intelligentSimplify();
            delete verifyNotThreeBall(tri,
                "Punctured Poincare homology sphere");

            // Throw in a couple of closed manifolds for good measure.
            tri = new NTriangulation();
            tri->insertLayeredLensSpace(1,0);
            delete verifyNotThreeBall(tri, "L(1,0)");

            tri = new NTriangulation();
            tri->insertLayeredLensSpace(2,1);
            delete verifyNotThreeBall(tri, "L(2,1)");

            tri = NExampleTriangulation::poincareHomologySphere();
            delete verifyNotThreeBall(tri, "Poincare homology sphere");
        }

        NTriangulation* verifySolidTorus(NTriangulation* tri,
                const std::string& triName) {
            NTriangulation bounded(*tri);
            if (bounded.isIdeal())
                bounded.idealToFinite();

            NTriangulation ideal(*tri);
            if (ideal.hasBoundaryFaces())
                ideal.finiteToIdeal();

            NTriangulation boundedBig(bounded);
            boundedBig.barycentricSubdivision();

            NTriangulation idealBig(ideal);
            idealBig.barycentricSubdivision();

            if (! bounded.isSolidTorus()) {
                CPPUNIT_FAIL(("The real solid torus " + triName +
                    " is not recognised as such.").c_str());
            }
            if (! ideal.isSolidTorus()) {
                CPPUNIT_FAIL(("The ideal solid torus " + triName +
                    " is not recognised as such.").c_str());
            }
            if (! boundedBig.isSolidTorus()) {
                CPPUNIT_FAIL(("The subdivided real solid torus " + triName +
                    " is not recognised as such.").c_str());
            }
            if (! idealBig.isSolidTorus()) {
                CPPUNIT_FAIL(("The subdivided ideal solid torus " + triName +
                    " is not recognised as such.").c_str());
            }

            return tri;
        }

        NTriangulation* verifyNotSolidTorus(NTriangulation* tri,
                const std::string& triName) {
            NTriangulation bounded(*tri);
            if (bounded.isIdeal())
                bounded.idealToFinite();

            NTriangulation ideal(*tri);
            if (ideal.hasBoundaryFaces())
                ideal.finiteToIdeal();

            NTriangulation boundedBig(bounded);
            boundedBig.barycentricSubdivision();

            NTriangulation idealBig(ideal);
            idealBig.barycentricSubdivision();

            if (bounded.isSolidTorus()) {
                CPPUNIT_FAIL(("The real non-solid-torus " + triName +
                    " is recognised as a solid torus.").c_str());
            }
            if (ideal.isSolidTorus()) {
                CPPUNIT_FAIL(("The ideal non-solid-torus " + triName +
                    " is recognised as a solid torus.").c_str());
            }
            if (boundedBig.isSolidTorus()) {
                CPPUNIT_FAIL(("The subdivided real non-solid-torus " +
                    triName + " is recognised as a solid torus.").c_str());
            }
            if (idealBig.isSolidTorus()) {
                CPPUNIT_FAIL(("The subdivided ideal non-solid-torus " +
                    triName + " is recognised as a solid torus.").c_str());
            }

            return tri;
        }

        void solidTorusRecognition() {
            NTriangulation* tri;

            tri = new NTriangulation();
            delete verifyNotSolidTorus(tri, "Empty triangulation");

            tri = new NTriangulation();
            tri->newTetrahedron();
            delete verifyNotSolidTorus(tri, "Single tetrahedron");

            tri = new NTriangulation();
            NTetrahedron* tet = tri->newTetrahedron();
            tet->joinTo(0, tet, NPerm4(3, 1, 2, 0));
            delete verifyNotSolidTorus(tri, "Snapped tetrahedron");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(1, 2);
            delete verifySolidTorus(tri, "LST(1,2,3)");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(1, 20);
            delete verifySolidTorus(tri, "LST(1,20,21)");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(1, 1);
            delete verifySolidTorus(tri, "LST(1,1,2)");

            tri = new NTriangulation();
            tri->insertLayeredSolidTorus(0, 1);
            delete verifySolidTorus(tri, "LST(0,1,1)");

            tri = NTriangulation::fromIsoSig("cMcabbgds");
            delete verifySolidTorus(tri, "Ideal solid torus");

            tri = NExampleTriangulation::figureEightKnotComplement();
            delete verifyNotSolidTorus(tri, "Figure 8 Knot Complement");

            // Throw in a couple of closed manifolds for good measure.
            tri = new NTriangulation();
            tri->insertLayeredLensSpace(1,0);
            delete verifyNotSolidTorus(tri, "L(1,0)");

            tri = new NTriangulation();
            tri->insertLayeredLensSpace(2,1);
            delete verifyNotSolidTorus(tri, "L(2,1)");

            tri = NExampleTriangulation::poincareHomologySphere();
            delete verifyNotSolidTorus(tri, "Poincare homology sphere");
        }
};

void addThreeSphere(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(ThreeSphereTest::suite());
}

