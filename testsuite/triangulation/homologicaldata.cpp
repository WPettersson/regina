
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Test Suite                                                            *
 *                                                                        *
 *  Copyright (c) 1999-2006, Ben Burton                                   *
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

#include <sstream>
#include <cppunit/extensions/HelperMacros.h>
#include "triangulation/nhomologicaldata.h"
#include "triangulation/ntriangulation.h"
#include "testsuite/triangulation/testtriangulation.h"

using regina::NHomologicalData;
using regina::NTriangulation;

/**
 * Thanks to Ryan Budney for supplying the questions and answers for
 * these tests!
 */
class HomologicalDataTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(HomologicalDataTest);

    CPPUNIT_TEST(bdryManifoldMapH1);
    CPPUNIT_TEST(standardCells);
    CPPUNIT_TEST(dualCells);
    CPPUNIT_TEST(torsionRankVector);
    CPPUNIT_TEST(torsionSigmaVector);
    CPPUNIT_TEST(torsionLegendreSymbolVector);

    CPPUNIT_TEST_SUITE_END();

    private:
        // Closed orientable:
        NTriangulation lens7_1;
            /**< The layered lens space L(7,1). */
        NTriangulation d88xz15;
            /**< The orbit manifold S^3 / D_88 x Z_15, also known
                 as SFS [S2: (2,1) (2,1) (11,19)]. */
        NTriangulation closedHypA;
            /**< The manifold with volume 1.01494161 from the
                 Hodgson-Weeks closed orientable census. */
        NTriangulation closedHypB;
            /**< The manifold with volume 2.45402944 from the
                 Hodgson-Weeks closed orientable census. */

        // Closed non-orientable:
        NTriangulation norA;
            /**< The non-orientable manifold SFS [M_/n2: (2,1)]. */

        // Ideal orientable:
        NTriangulation s028;
            /**< The manifold s028 from the cusped hyperbolic census of
                 Callahan, Hildebrand and Weeks, as shipped with SnapPea. */
        NTriangulation s955;
            /**< The manifold s955 from the cusped hyperbolic census of
                 Callahan, Hildebrand and Weeks, as shipped with SnapPea. */

        // Bounded:
        NTriangulation lst3_4_7;
            /**< The layered solid torus LST(3,4,7). */

    public:
        void setUp() {
            // Deal with the easy ones first.
            lens7_1.insertLayeredLensSpace(7, 1);
            d88xz15.insertAugTriSolidTorus(2, -1, 2, -1, 11, -30);
            lst3_4_7.insertLayeredSolidTorus(3, 4);

            // The others we recreate using dehydration strings.
            // The dehydration strings were obtained by running
            // NTriangulation::dehydrate() over triangulations from
            // the relevant census data files.
            closedHypA.insertRehydration("jgpadaaeffghfiihirmxitdagbj");
            closedHypB.insertRehydration("mnnjjcabehfgjijkijllaaqabhoehrtab");
            norA.insertRehydration("jofbdaabccfhgihiiffhofoxohx");
            s028.insertRehydration("gkfacaccdeffffohhhf");
            s955.insertRehydration("gbpaabcfdffefohfxhf");
        }

        void tearDown() {
        }

        void verifyBdryManifoldMapH1(NTriangulation& tri, const char* name,
                const char* ans) {
            NHomologicalData dat(tri);
            std::string val = dat.getBMmapH(1).toString();
            if (val != ans) {
                std::ostringstream msg;
                msg << name << ": Map from H1(bdry) to H1(mfd) is "
                    << val << ", not " << ans << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void bdryManifoldMapH1() {
            verifyBdryManifoldMapH1(lens7_1, "L(7,1)", "zero map");
            verifyBdryManifoldMapH1(d88xz15, "S^3 / D_88 x Z_15", "zero map");
            verifyBdryManifoldMapH1(closedHypA, "Closed Hyp (vol=1.01494161)",
                "zero map");
            verifyBdryManifoldMapH1(closedHypB, "Closed Hyp (vol=2.45402944)",
                "zero map");
            verifyBdryManifoldMapH1(norA, "SFS [M_/n2: (2,1)]", "zero map");
            verifyBdryManifoldMapH1(s028, "SnapPea s028",
                "kernel Z | cokernel Z_16 | image Z");
            verifyBdryManifoldMapH1(s955, "SnapPea s955",
                "kernel Z | cokernel Z_20 | image Z + Z_2");
            verifyBdryManifoldMapH1(lst3_4_7, "LST(3,4,7)",
                "epic, with kernel Z");
        }

        void verifyStandardCells(NTriangulation& tri, const char* name,
                unsigned long c0, unsigned long c1, unsigned long c2,
                unsigned long c3) {
            NHomologicalData dat(tri);
            std::vector<unsigned long> vals = dat.getNumStandardCells();

            if (vals.size() != 4) {
                std::ostringstream msg;
                msg << name << ": Standard cell count contains "
                    << vals.size() << " element(s), not four.";
                CPPUNIT_FAIL(msg.str());
            }

            unsigned long ans[4] = { c0, c1, c2, c3 };
            for (int i = 0; i < 4; i++)
                if (vals[i] != ans[i]) {
                    std::ostringstream msg;
                    msg << name << ": Number of standard cells of dimension "
                        << i << " is " << vals[i] << ", not " << ans[i] << ".";
                    CPPUNIT_FAIL(msg.str());
                }
        }

        void standardCells() {
            verifyStandardCells(lens7_1, "L(7,1)", 1, 5, 8, 4);
            verifyStandardCells(d88xz15, "S^3 / D_88 x Z_15", 1, 10, 18, 9);
            verifyStandardCells(closedHypA, "Closed Hyp (vol=1.01494161)",
                1, 10, 18, 9);
            verifyStandardCells(closedHypB, "Closed Hyp (vol=2.45402944)",
                1, 13, 24, 12);
            verifyStandardCells(norA, "SFS [M_/n2: (2,1)]", 1, 10, 18, 9);
            verifyStandardCells(s028, "SnapPea s028", 12, 42, 36, 6);
            verifyStandardCells(s955, "SnapPea s955", 12, 42, 36, 6);
            verifyStandardCells(lst3_4_7, "LST(3,4,7)", 1, 5, 7, 3);
        }

        void verifyDualCells(NTriangulation& tri, const char* name,
                unsigned long c0, unsigned long c1, unsigned long c2,
                unsigned long c3) {
            NHomologicalData dat(tri);
            std::vector<unsigned long> vals = dat.getNumDualCells();

            if (vals.size() != 4) {
                std::ostringstream msg;
                msg << name << ": Dual cell count contains "
                    << vals.size() << " element(s), not four.";
                CPPUNIT_FAIL(msg.str());
            }

            unsigned long ans[4] = { c0, c1, c2, c3 };
            for (int i = 0; i < 4; i++)
                if (vals[i] != ans[i]) {
                    std::ostringstream msg;
                    msg << name << ": Number of dual cells of dimension "
                        << i << " is " << vals[i] << ", not " << ans[i] << ".";
                    CPPUNIT_FAIL(msg.str());
                }
        }

        void dualCells() {
            verifyDualCells(lens7_1, "L(7,1)", 4, 8, 5, 1);
            verifyDualCells(d88xz15, "S^3 / D_88 x Z_15", 9, 18, 10, 1);
            verifyDualCells(closedHypA, "Closed Hyp (vol=1.01494161)",
                9, 18, 10, 1);
            verifyDualCells(closedHypB, "Closed Hyp (vol=2.45402944)",
                12, 24, 13, 1);
            verifyDualCells(norA, "SFS [M_/n2: (2,1)]", 9, 18, 10, 1);
            verifyDualCells(s028, "SnapPea s028", 6, 12, 6, 0);
            verifyDualCells(s955, "SnapPea s955", 6, 12, 6, 0);
            verifyDualCells(lst3_4_7, "LST(3,4,7)", 3, 5, 2, 0);
        }

        void verifyTorsionRankVector(NTriangulation& tri,
                const char* name, const char* ans) {
            NHomologicalData dat(tri);
            std::string val = dat.getTorsionRankVectorString();
            if (val != ans) {
                std::ostringstream msg;
                msg << name << ": Torsion subgroup rank vector is "
                    << val << ", not " << ans << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void torsionRankVector() {
            verifyTorsionRankVector(lens7_1, "L(7,1)", "7(1)");
            verifyTorsionRankVector(d88xz15, "S^3 / D_88 x Z_15",
                "2(0 0 1) 3(1) 5(1)");
            verifyTorsionRankVector(closedHypA, "Closed Hyp (vol=1.01494161)",
                "2(1) 3(2)");
            verifyTorsionRankVector(closedHypB, "Closed Hyp (vol=2.45402944)",
                "7(2)");
            verifyTorsionRankVector(norA, "SFS [M_/n2: (2,1)]", "2(1 1)");
            verifyTorsionRankVector(s028, "SnapPea s028", "2(0 0 0 1)");
            verifyTorsionRankVector(s955, "SnapPea s955", "2(0 1) 5(1)");
            verifyTorsionRankVector(lst3_4_7, "LST(3,4,7)", "no torsion");
        }

        void verifyTorsionSigmaVector(NTriangulation& tri,
                const char* name, const char* ans) {
            NHomologicalData dat(tri);
            std::string val = dat.getTorsionSigmaVectorString();
            if (val != ans) {
                std::ostringstream msg;
                msg << name << ": 2-torsion sigma vector is "
                    << val << ", not " << ans << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void torsionSigmaVector() {
            verifyTorsionSigmaVector(lens7_1, "L(7,1)", "no 2-torsion");
            verifyTorsionSigmaVector(d88xz15, "S^3 / D_88 x Z_15", "7 7 inf");
            verifyTorsionSigmaVector(closedHypA, "Closed Hyp (vol=1.01494161)",
                "inf");
            verifyTorsionSigmaVector(closedHypB, "Closed Hyp (vol=2.45402944)",
                "no 2-torsion");
            verifyTorsionSigmaVector(norA, "SFS [M_/n2: (2,1)]",
                "manifold is non-orientable");
            verifyTorsionSigmaVector(s028, "SnapPea s028", "7 7 7 inf");
            verifyTorsionSigmaVector(s955, "SnapPea s955", "inf 0");
            verifyTorsionSigmaVector(lst3_4_7, "LST(3,4,7)", "no 2-torsion");
        }

        void verifyTorsionLegendreSymbolVector(NTriangulation& tri,
                const char* name, const char* ans) {
            NHomologicalData dat(tri);
            std::string val = dat.getTorsionLegendreSymbolVectorString();
            if (val != ans) {
                std::ostringstream msg;
                msg << name << ": Odd p-torsion Legendre symbol vector is "
                    << val << ", not " << ans << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void torsionLegendreSymbolVector() {
            verifyTorsionLegendreSymbolVector(lens7_1, "L(7,1)", "7(1)");
            verifyTorsionLegendreSymbolVector(d88xz15, "S^3 / D_88 x Z_15",
                "3(1) 5(1)");
            verifyTorsionLegendreSymbolVector(closedHypA,
                "Closed Hyp (vol=1.01494161)", "3(1)");
            verifyTorsionLegendreSymbolVector(closedHypB,
                "Closed Hyp (vol=2.45402944)", "7(-1)");
            verifyTorsionLegendreSymbolVector(norA, "SFS [M_/n2: (2,1)]",
                "manifold is non-orientable");
            verifyTorsionLegendreSymbolVector(s028, "SnapPea s028",
                "no odd p-torsion");
            verifyTorsionLegendreSymbolVector(s955, "SnapPea s955", "5(-1)");
            verifyTorsionLegendreSymbolVector(lst3_4_7, "LST(3,4,7)",
                "no odd p-torsion");
        }
};

void addHomologicalData(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(HomologicalDataTest::suite());
}
