
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

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>
#include <cppunit/extensions/HelperMacros.h>
#include "algebra/nabeliangroup.h"
#include "algebra/ngrouppresentation.h"
#include "maths/approx.h"
#include "maths/numbertheory.h"
#include "packet/ncontainer.h"
#include "split/nsignature.h"
#include "subcomplex/nstandardtri.h"
#include "triangulation/nexampletriangulation.h"
#include "triangulation/nisomorphism.h"
#include "triangulation/ntriangulation.h"
#include "triangulation/nvertex.h"
#include "testsuite/triangulation/testtriangulation.h"

using regina::NAbelianGroup;
using regina::NExampleTriangulation;
using regina::NGroupPresentation;
using regina::NIsomorphism;
using regina::NPerm4;
using regina::NSignature;
using regina::NStandardTriangulation;
using regina::NTetrahedron;
using regina::NTriangulation;
using regina::NVertex;

class NTriangulationTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(NTriangulationTest);

    CPPUNIT_TEST(validity);
    CPPUNIT_TEST(standardness);
    CPPUNIT_TEST(orientability);
    CPPUNIT_TEST(boundaryComponents);
    CPPUNIT_TEST(vertexLinks);
    CPPUNIT_TEST(eulerCharacteristic);
    CPPUNIT_TEST(homologyH1);
    CPPUNIT_TEST(homologyH1Bdry);
    CPPUNIT_TEST(fundGroup);
    CPPUNIT_TEST(zeroEfficiency);
    CPPUNIT_TEST(turaevViro);
    CPPUNIT_TEST(doubleCover);
    CPPUNIT_TEST(barycentricSubdivision);
    CPPUNIT_TEST(idealToFinite);
    CPPUNIT_TEST(finiteToIdeal);
    CPPUNIT_TEST(drillEdge);
    CPPUNIT_TEST(dehydration);
    CPPUNIT_TEST(isomorphismSignature);
    CPPUNIT_TEST(simplification);
    CPPUNIT_TEST(reordering);
    CPPUNIT_TEST(propertyUpdates);

    CPPUNIT_TEST_SUITE_END();

    private:
        // Trivial:
        NTriangulation empty;
            /**< An empty triangulation. */
        NTriangulation singleTet;
            /**< A single tetrahedron with no face gluings. */

        // Closed orientable:
        NTriangulation s3;
            /**< A one-vertex 3-sphere. */
        NTriangulation s2xs1;
            /**< The product space S^2 x S^1. */
        NTriangulation rp3_1;
            /**< A one-vertex triangulation of RP^3. */
        NTriangulation rp3_2;
            /**< A two-vertex triangulation of RP^3. */
        NTriangulation lens3_1;
            /**< A two-vertex lens space L(3,1). */
        NTriangulation lens8_3;
            /**< The layered lens space L(8,3). */
        NTriangulation lens7_1_loop;
            /**< An untwisted layered loop representing L(7,1). */
        NTriangulation rp3rp3;
            /**< The connected sum RP^3 # RP^3. */
        NTriangulation q32xz3;
            /**< The orbit manifold S^3 / Q_32 x Z_3. */
        NTriangulation q28;
            /**< A twisted layered loop representing S^3 / Q_28. */
        NTriangulation weberSeifert;
            /**< The Weber-Seifert dodecahedral space.  With 23 tetrahedra,
                 this is too large for running normal surface algorithms.*/

        // Closed orientable, very large:
        NTriangulation lens100_1;
            /**< The lens space L(100,1). */

        // Bounded orientable:
        NTriangulation lst3_4_7;
            /**< The layered solid torus LST(3,4,7). */

        // Ideal orientable:
        NTriangulation figure8;
            /**< The figure eight knot complement. */

        // Closed non-orientable:
        NTriangulation rp2xs1;
            /**< The product RP^2 x S^1. */

        // Bounded non-orientable:
        NTriangulation solidKB;
            /**< A solid Klein bottle. */

        // Ideal non-orientable:
        NTriangulation gieseking;
            /**< The Gieseking manifold. */

        // Invalid edges:
        NTriangulation invalidEdges;
            /**< Contains two invalid edges, but all vertices have
                 2-sphere links. */

        // Non-standard vertex links:
        NTriangulation twoProjPlaneCusps;
            /**< A subdivision of invalidEdges, resulting in all edges
                 valid but two projective plane cusps.  Note that this
                 triangulation has a 3-sphere orientable double cover. */
        NTriangulation cuspedGenusTwoTorus;
            /**< A solid genus two torus with a cusped boundary. */
        NTriangulation pinchedSolidTorus;
            /**< A solid torus with one longitude pinched to a point. */
        NTriangulation pinchedSolidKB;
            /**< A solid Klein bottle with one longitude pinched to a point. */

        // Non-minimal:
        NTriangulation s3_large;
            /**< A non-minimal 3-sphere triangulation. */
        NTriangulation rp3_large;
            /**< The manifold RP^3. */
        NTriangulation lens8_3_large;
            /**< The lens space L(8,3). */
        NTriangulation q20_large;
            /**< The manifold S^3 / Q_20. */
        NTriangulation ball_large;
            /**< A four tetrahedron ball that used to cause crashes in
                 the simplification routines once upon a time. */
        NTriangulation ball_large_pillows;
            /**< A four tetrahedron ball formed by joining together two
                 triangular pillows. */
        NTriangulation ball_large_snapped;
            /**< A ball formed by joining together three snapped balls. */
        NTriangulation singleTet_bary;
            /**< The barycentric subdivision of a single tetrahedron. */
        NTriangulation fig8_bary;
            /**< The barycentric subdivision of the figure eight knot
                 complement. */

    public:
        void copyAndDelete(NTriangulation& dest, NTriangulation* source) {
            dest.insertTriangulation(*source);
            delete source;
        }

        void generateFromSig(NTriangulation& tri, const std::string& sigStr) {
            NSignature* sig = NSignature::parse(sigStr);
            if (sig == 0)
                return;

            NTriangulation* triNew = sig->triangulate();
            delete sig;
            if (triNew == 0)
                return;

            tri.insertTriangulation(*triNew);
            delete triNew;
        }

        void setUp() {
            // Begin with trivial cases.
            // The empty triangulation needs no initialisation whatsoever.
            singleTet.newTetrahedron();

            // Some of our triangulations can be constructed automatically.
            s3.insertLayeredLensSpace(1, 0);
            s2xs1.insertLayeredLensSpace(0, 1);
            rp3_1.insertLayeredLensSpace(2, 1);
            rp3_2.insertLayeredLoop(2, false);
            lens8_3.insertLayeredLensSpace(8, 3);
            lens100_1.insertLayeredLensSpace(100, 1);
            lst3_4_7.insertLayeredSolidTorus(3, 4);
            q28.insertLayeredLoop(7, true);
            lens7_1_loop.insertLayeredLoop(7, false);

            // Some of our triangulations can be generated from
            // splitting surfaces.
            generateFromSig(rp3rp3, "aabccd.b.d");
            generateFromSig(q32xz3, "aabcdb.cedfef");
            generateFromSig(s3_large, "abc.abd.cef.de.fg.g");
            generateFromSig(lens8_3_large, "aabcb.cd.d");
            generateFromSig(rp3_large, "aabcdedcfb.fg.e.g");
            generateFromSig(q20_large, "abcdeabcdef.fg.g");

            // Some are hard-coded in the calculation engine as sample
            // triangulations.
            copyAndDelete(weberSeifert, NExampleTriangulation::weberSeifert());
            copyAndDelete(figure8,
                NExampleTriangulation::figureEightKnotComplement());
            copyAndDelete(solidKB, NExampleTriangulation::solidKleinBottle());
            copyAndDelete(rp2xs1, NExampleTriangulation::rp2xs1());
            copyAndDelete(gieseking, NExampleTriangulation::gieseking());
            copyAndDelete(cuspedGenusTwoTorus,
                NExampleTriangulation::cuspedGenusTwoTorus());

            singleTet_bary.newTetrahedron();
            singleTet_bary.barycentricSubdivision();

            copyAndDelete(fig8_bary,
                NExampleTriangulation::figureEightKnotComplement());
            fig8_bary.barycentricSubdivision();

            // The rest alas must be done manually.
            NTetrahedron* r;
            NTetrahedron* s;
            NTetrahedron* t;
            NTetrahedron* u;

            // A two-tetrahedron two-vertex L(3,1) is straightforward to
            // construct using a vertex of degree two.
            r = lens3_1.newTetrahedron();
            s = lens3_1.newTetrahedron();
            r->joinTo(0, s, NPerm4(0, 2, 3, 1));
            r->joinTo(1, s, NPerm4());
            r->joinTo(2, s, NPerm4());
            r->joinTo(3, s, NPerm4());

            // For a triangulation with invalid edges, we simply fold
            // the faces of a tetrahedron together in pairs (as in a
            // 3-sphere triangulation) but apply a reflection to each fold.
            r = invalidEdges.newTetrahedron();
            r->joinTo(0, r, NPerm4(1, 0, 3, 2));
            r->joinTo(2, r, NPerm4(1, 0, 3, 2));

            twoProjPlaneCusps.insertTriangulation(invalidEdges);
            twoProjPlaneCusps.barycentricSubdivision();

            // To construct a solid torus with a pinched longitude, we
            // identify two opposite faces of a square pyramid.
            r = pinchedSolidTorus.newTetrahedron();
            s = pinchedSolidTorus.newTetrahedron();
            r->joinTo(3, s, NPerm4(0, 1, 2, 3));
            r->joinTo(2, s, NPerm4(0, 3, 1, 2));

            // The pinched solid Klein bottle is much the same, except
            // for a twist before the opposite faces are identified.
            r = pinchedSolidKB.newTetrahedron();
            s = pinchedSolidKB.newTetrahedron();
            r->joinTo(3, s, NPerm4(0, 1, 2, 3));
            r->joinTo(2, s, NPerm4(0, 2, 1, 3));

            // This ball used to cause a crash once upon a time.
            // Throw it into the test suite for good measure.
            r = ball_large.newTetrahedron();
            s = ball_large.newTetrahedron();
            t = ball_large.newTetrahedron();
            u = ball_large.newTetrahedron();
            r->joinTo(2, r, NPerm4(0,2));
            r->joinTo(1, s, NPerm4(2,0,1,3));
            s->joinTo(2, t, NPerm4());
            s->joinTo(1, t, NPerm4(2,0,1,3));
            t->joinTo(1, u, NPerm4(2,0,1,3));
            u->joinTo(2, u, NPerm4(1,2));

            // Make two triangular pillows, then join them together.
            // This crashed with 2-0 vertex moves once upon a time.
            r = ball_large_pillows.newTetrahedron();
            s = ball_large_pillows.newTetrahedron();
            t = ball_large_pillows.newTetrahedron();
            u = ball_large_pillows.newTetrahedron();
            r->joinTo(0, s, NPerm4());
            r->joinTo(1, s, NPerm4());
            r->joinTo(2, s, NPerm4());
            t->joinTo(0, u, NPerm4());
            t->joinTo(1, u, NPerm4());
            t->joinTo(2, u, NPerm4());
            r->joinTo(3, t, NPerm4());

            // Make three snapped balls and join them together.
            r = ball_large_snapped.newTetrahedron();
            s = ball_large_snapped.newTetrahedron();
            t = ball_large_snapped.newTetrahedron();
            r->joinTo(2, r, NPerm4(2, 3));
            s->joinTo(2, s, NPerm4(2, 3));
            t->joinTo(2, t, NPerm4(2, 1));
            r->joinTo(1, s, NPerm4());
            s->joinTo(0, t, NPerm4());
        }

        void tearDown() {
        }

        void validity() {
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not valid.",
                empty.isValid());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron is not valid.",
                singleTet.isValid());
            CPPUNIT_ASSERT_MESSAGE("S^3 is not valid.",
                s3.isValid());
            CPPUNIT_ASSERT_MESSAGE("S^2 x S^1 is not valid.",
                s2xs1.isValid());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (1 vtx) is not valid.",
                rp3_1.isValid());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (2 vtx) is not valid.",
                rp3_2.isValid());
            CPPUNIT_ASSERT_MESSAGE("L(3,1) is not valid.",
                lens3_1.isValid());
            CPPUNIT_ASSERT_MESSAGE("Layered loop L(7,1) is not valid.",
                lens7_1_loop.isValid());
            CPPUNIT_ASSERT_MESSAGE("L(8,3) is not valid.",
                lens8_3.isValid());
            CPPUNIT_ASSERT_MESSAGE("Large L(8,3) is not valid.",
                lens8_3_large.isValid());
            CPPUNIT_ASSERT_MESSAGE("RP^3 # RP^3 is not valid.",
                rp3rp3.isValid());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_28 is not valid.",
                q28.isValid());
            CPPUNIT_ASSERT_MESSAGE(
                "The Weber-Seifert dodecahedral space is not valid.",
                weberSeifert.isValid());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_32 x Z_3 is not valid.",
                q32xz3.isValid());
            CPPUNIT_ASSERT_MESSAGE("L(100,1) is not valid.",
                lens100_1.isValid());
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron ball is not valid.",
                ball_large.isValid());
            CPPUNIT_ASSERT_MESSAGE(
                "The 4-tetrahedron pillow ball is not valid.",
                ball_large_pillows.isValid());
            CPPUNIT_ASSERT_MESSAGE(
                "The 3-tetrahedron snapped ball is not valid.",
                ball_large_snapped.isValid());
            CPPUNIT_ASSERT_MESSAGE("LST(3,4,7) is not valid.",
                lst3_4_7.isValid());
            CPPUNIT_ASSERT_MESSAGE("The figure eight knot complement "
                "is not valid.",
                figure8.isValid());
            CPPUNIT_ASSERT_MESSAGE("RP^2 x S^1 is not valid.",
                rp2xs1.isValid());
            CPPUNIT_ASSERT_MESSAGE("The solid Klein bottle is not valid.",
                solidKB.isValid());
            CPPUNIT_ASSERT_MESSAGE("The Gieseking manifold is not valid.",
                gieseking.isValid());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with invalid edges "
                "is reported as valid.",
                ! invalidEdges.isValid());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with projective plane "
                "cusps is not valid.",
                twoProjPlaneCusps.isValid());
            CPPUNIT_ASSERT_MESSAGE("The cusped solid genus two torus "
                "is not valid.",
                cuspedGenusTwoTorus.isValid());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid torus "
                "is reported as valid.",
                ! pinchedSolidTorus.isValid());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid Klein bottle "
                "is reported as valid.",
                ! pinchedSolidKB.isValid());
        }

        void standardness() {
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not standard.",
                empty.isStandard());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron is not standard.",
                singleTet.isStandard());
            CPPUNIT_ASSERT_MESSAGE("S^3 is not standard.",
                s3.isStandard());
            CPPUNIT_ASSERT_MESSAGE("S^2 x S^1 is not standard.",
                s2xs1.isStandard());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (1 vtx) is not standard.",
                rp3_1.isStandard());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (2 vtx) is not standard.",
                rp3_2.isStandard());
            CPPUNIT_ASSERT_MESSAGE("L(3,1) is not standard.",
                lens3_1.isStandard());
            CPPUNIT_ASSERT_MESSAGE("Layered loop L(7,1) is not standard.",
                lens7_1_loop.isStandard());
            CPPUNIT_ASSERT_MESSAGE("L(8,3) is not standard.",
                lens8_3.isStandard());
            CPPUNIT_ASSERT_MESSAGE("Large L(8,3) is not standard.",
                lens8_3_large.isStandard());
            CPPUNIT_ASSERT_MESSAGE("RP^3 # RP^3 is not standard.",
                rp3rp3.isStandard());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_28 is not standard.",
                q28.isStandard());
            CPPUNIT_ASSERT_MESSAGE(
                "The Weber-Seifert dodecahedral space is not standard.",
                weberSeifert.isStandard());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_32 x Z_3 is not standard.",
                q32xz3.isStandard());
            CPPUNIT_ASSERT_MESSAGE("L(100,1) is not standard.",
                lens100_1.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron ball is not standard.",
                ball_large.isStandard());
            CPPUNIT_ASSERT_MESSAGE(
                "The 4-tetrahedron pillow ball is not standard.",
                ball_large_pillows.isStandard());
            CPPUNIT_ASSERT_MESSAGE(
                "The 3-tetrahedron snapped ball is not standard.",
                ball_large_snapped.isStandard());
            CPPUNIT_ASSERT_MESSAGE("LST(3,4,7) is not standard.",
                lst3_4_7.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The figure eight knot complement "
                "is not standard.",
                figure8.isStandard());
            CPPUNIT_ASSERT_MESSAGE("RP^2 x S^1 is not standard.",
                rp2xs1.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The solid Klein bottle is not standard.",
                solidKB.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The Gieseking manifold is not standard.",
                gieseking.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with invalid edges "
                "is not standard.",
                invalidEdges.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with projective plane "
                "cusps is standard.",
                ! twoProjPlaneCusps.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The cusped solid genus two torus "
                "is standard.",
                ! cuspedGenusTwoTorus.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid torus is standard.",
                ! pinchedSolidTorus.isStandard());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid Klein bottle "
                "is standard.",
                ! pinchedSolidKB.isStandard());
        }

        void orientability() {
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not orientable.",
                empty.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron is not orientable.",
                singleTet.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("S^3 is not orientable.",
                s3.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("S^2 x S^1 is not orientable.",
                s2xs1.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (1 vtx) is not orientable.",
                rp3_1.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (2 vtx) is not orientable.",
                rp3_2.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("L(3,1) is not orientable.",
                lens3_1.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("Layered loop L(7,1) is not orientable.",
                lens7_1_loop.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("L(8,3) is not orientable.",
                lens8_3.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("Large L(8,3) is not orientable.",
                lens8_3_large.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("RP^3 # RP^3 is not orientable.",
                rp3rp3.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_28 is not orientable.",
                q28.isOrientable());
            CPPUNIT_ASSERT_MESSAGE(
                "The Weber-Seifert dodecahedral space is not orientable.",
                weberSeifert.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_32 x Z_3 is not orientable.",
                q32xz3.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("L(100,1) is not orientable.",
                lens100_1.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron ball is not orientable.",
                ball_large.isOrientable());
            CPPUNIT_ASSERT_MESSAGE(
                "The 4-tetrahedron pillow ball is not orientable.",
                ball_large_pillows.isOrientable());
            CPPUNIT_ASSERT_MESSAGE(
                "The 3-tetrahedron snapped ball is not orientable.",
                ball_large_snapped.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("LST(3,4,7) is not orientable.",
                lst3_4_7.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The figure eight knot complement "
                "is not orientable.",
                figure8.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("RP^2 x S^1 is orientable.",
                ! rp2xs1.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The solid Klein bottle is orientable.",
                ! solidKB.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The Gieseking manifold is orientable.",
                ! gieseking.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with invalid edges "
                "is orientable.",
                ! invalidEdges.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The triangulation with projective plane "
                "cusps is orientable.",
                ! twoProjPlaneCusps.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The cusped solid genus two torus "
                "is not orientable.",
                cuspedGenusTwoTorus.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid torus is not orientable.",
                pinchedSolidTorus.isOrientable());
            CPPUNIT_ASSERT_MESSAGE("The pinched solid Klein bottle "
                "is orientable.",
                ! pinchedSolidKB.isOrientable());
        }

        void boundaryComponents() {
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation has boundary "
                "components.", empty.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron has no boundary "
                "components.", singleTet.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("S^3 has boundary components.",
                s3.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("S^2 x S^1 has boundary components.",
                s2xs1.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("RP^3 (1 vtx) has boundary components.",
                rp3_1.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("RP^3 (2 vtx) has boundary components.",
                rp3_2.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("L(3,1) has boundary components.",
                lens3_1.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("Layered loop L(7,1) "
                "has boundary components.",
                lens7_1_loop.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("L(8,3) has boundary components.",
                lens8_3.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("Large L(8,3) has boundary components.",
                lens8_3_large.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("RP^3 # RP^3 has boundary components.",
                rp3rp3.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_28 has boundary components.",
                q28.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE(
                "The Weber-Seifert dodecahedral space has boundary components.",
                weberSeifert.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_32 x Z_3 has boundary components.",
                q32xz3.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("L(100,1) has boundary components.",
                lens100_1.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron ball has no boundary "
                "components.", ball_large.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron pillow ball has no "
                "boundary components.",
                ball_large_pillows.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The 3-tetrahedron snapped ball has no "
                "boundary components.",
                ball_large_snapped.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("LST(3,4,7) has no boundary components.",
                lst3_4_7.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The figure eight knot complement "
                "has no boundary components.",
                figure8.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("RP^2 x S^1 has boundary components.",
                rp2xs1.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("The solid Klein bottle "
                "has no boundary components.",
                solidKB.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The Gieseking manifold "
                "has no boundary components.",
                gieseking.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The triangulation with invalid edges "
                "has boundary components.",
                invalidEdges.getNumberOfBoundaryComponents() == 0);
            CPPUNIT_ASSERT_MESSAGE("The triangulation with projective plane "
                "cusps has no boundary components.",
                twoProjPlaneCusps.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The cusped solid genus two torus "
                "has no boundary components.",
                cuspedGenusTwoTorus.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The pinched solid torus "
                "has no boundary components.",
                pinchedSolidTorus.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The pinched solid torus "
                "has too many boundary components.",
                pinchedSolidTorus.getNumberOfBoundaryComponents() < 2);
            CPPUNIT_ASSERT_MESSAGE("The pinched solid Klein bottle "
                "has no boundary components.",
                pinchedSolidKB.getNumberOfBoundaryComponents() > 0);
            CPPUNIT_ASSERT_MESSAGE("The pinched solid Klein bottle "
                "has too many boundary components.",
                pinchedSolidKB.getNumberOfBoundaryComponents() < 2);

            // TODO: Test the individual boundary components.
            // TODO: Check that nobody has too many boundary components.
        }

        void verifyVertexCount(NTriangulation& tri, unsigned nVertices,
                const char* triName) {
            if (tri.getNumberOfVertices() != nVertices) {
                std::ostringstream msg;
                msg << triName << " has " << tri.getNumberOfVertices()
                    << " vertices, not " << nVertices << '.';
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexDisc(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::DISC) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as DISC.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 1) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 1.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexSphere(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::SPHERE) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as SPHERE.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 2) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 2.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexTorus(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::TORUS) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as TORUS.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 0) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 0.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexKB(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::KLEIN_BOTTLE) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as KLEIN_BOTTLE.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 0) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 0.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexTorusG2(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::NON_STANDARD_CUSP) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as NON_STANDARD_CUSP.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != -2) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not -2.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexProjPlane(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::NON_STANDARD_CUSP) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as NON_STANDARD_CUSP.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 1) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 1.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexAnnulus(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::NON_STANDARD_BDRY) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as NON_STANDARD_BDRY.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is non-orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 0) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 0.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyVertexMobius(NTriangulation& tri, unsigned vertex,
                const char* triName) {
            if (vertex >= tri.getNumberOfVertices()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " does not exist.  Only " << tri.getNumberOfVertices()
                    << " vertices are available.";
                CPPUNIT_FAIL(msg.str());
            }

            NVertex* v = tri.getVertex(vertex);

            if (v->getLink() != NVertex::NON_STANDARD_BDRY) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is not listed as NON_STANDARD_BDRY.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkClosed()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is closed.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isIdeal()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is ideal.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! v->isBoundary()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " is not a boundary vertex.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isStandard()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is standard.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->isLinkOrientable()) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link is orientable.";
                CPPUNIT_FAIL(msg.str());
            }
            if (v->getLinkEulerCharacteristic() != 0) {
                std::ostringstream msg;
                msg << triName << ", vertex " << vertex
                    << " link has Euler characteristic "
                    << v->getLinkEulerCharacteristic() << ", not 0.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void vertexLinks() {
            verifyVertexCount(singleTet, 4, "Single tetrahedron");
            verifyVertexDisc(singleTet, 0, "Single tetrahedron");
            verifyVertexDisc(singleTet, 1, "Single tetrahedron");
            verifyVertexDisc(singleTet, 2, "Single tetrahedron");
            verifyVertexDisc(singleTet, 3, "Single tetrahedron");

            verifyVertexCount(s3, 1, "S^3");
            verifyVertexSphere(s3, 0, "S^3");

            verifyVertexCount(s2xs1, 1, "S^2 x S^1");
            verifyVertexSphere(s2xs1, 0, "S^2 x S^1");

            verifyVertexCount(rp3_1, 1, "RP^3 (1 vtx)");
            verifyVertexSphere(rp3_1, 0, "RP^3 (1 vtx)");

            verifyVertexCount(rp3_2, 2, "RP^3 (2 vtx)");
            verifyVertexSphere(rp3_2, 0, "RP^3 (2 vtx)");
            verifyVertexSphere(rp3_2, 1, "RP^3 (2 vtx)");

            verifyVertexCount(lens3_1, 2, "L(3,1)");
            verifyVertexSphere(lens3_1, 0, "L(3,1)");
            verifyVertexSphere(lens3_1, 1, "L(3,1)");

            verifyVertexCount(lens8_3, 1, "L(8,3)");
            verifyVertexSphere(lens8_3, 0, "L(8,3)");

            verifyVertexCount(lens8_3_large, 1, "Large L(8,3)");
            verifyVertexSphere(lens8_3_large, 0, "Large L(8,3)");

            verifyVertexCount(lens7_1_loop, 2, "Layered loop L(7,1)");
            verifyVertexSphere(lens7_1_loop, 0, "Layered loop L(7,1)");
            verifyVertexSphere(lens7_1_loop, 1, "Layered loop L(7,1)");

            verifyVertexCount(rp3rp3, 1, "RP^3 # RP^3");
            verifyVertexSphere(rp3rp3, 0, "RP^3 # RP^3");

            verifyVertexCount(q32xz3, 1, "S^3 / Q_32 x Z_3");
            verifyVertexSphere(q32xz3, 0, "S^3 / Q_32 x Z_3");

            verifyVertexCount(q28, 1, "S^3 / Q_28");
            verifyVertexSphere(q28, 0, "S^3 / Q_28");

            verifyVertexCount(weberSeifert, 1,
                "Weber-Seifert dodecahedral space");
            verifyVertexSphere(weberSeifert, 0,
                "Weber-Seifert dodecahedral space");

            verifyVertexCount(lens100_1, 1, "L(100,1)");
            verifyVertexSphere(lens100_1, 0, "L(100,1)");

            verifyVertexCount(ball_large, 5, "4-tetrahedron ball");
            verifyVertexDisc(ball_large, 0, "4-tetrahedron ball");
            verifyVertexDisc(ball_large, 1, "4-tetrahedron ball");
            verifyVertexSphere(ball_large, 2, "4-tetrahedron ball");
            verifyVertexDisc(ball_large, 3, "4-tetrahedron ball");
            verifyVertexDisc(ball_large, 4, "4-tetrahedron ball");

            verifyVertexCount(ball_large_pillows, 5,
                "4-tetrahedron pillow ball");
            verifyVertexDisc(ball_large_pillows, 0,
                "4-tetrahedron pillow ball");
            verifyVertexDisc(ball_large_pillows, 1,
                "4-tetrahedron pillow ball");
            verifyVertexDisc(ball_large_pillows, 2,
                "4-tetrahedron pillow ball");
            verifyVertexSphere(ball_large_pillows, 3,
                "4-tetrahedron pillow ball");
            verifyVertexSphere(ball_large_pillows, 4,
                "4-tetrahedron pillow ball");

            verifyVertexCount(ball_large_snapped, 4,
                "3-tetrahedron snapped ball");
            verifyVertexSphere(ball_large_snapped, 0,
                "3-tetrahedron snapped ball");
            verifyVertexDisc(ball_large_snapped, 1,
                "3-tetrahedron snapped ball");
            verifyVertexDisc(ball_large_snapped, 2,
                "3-tetrahedron snapped ball");
            verifyVertexDisc(ball_large_snapped, 3,
                "3-tetrahedron snapped ball");

            verifyVertexCount(lst3_4_7, 1, "LST(3,4,7)");
            verifyVertexDisc(lst3_4_7, 0, "LST(3,4,7)");

            verifyVertexCount(figure8, 1, "Figure eight knot complement");
            verifyVertexTorus(figure8, 0, "Figure eight knot complement");

            verifyVertexCount(rp2xs1, 1, "RP^2 x S^1");
            verifyVertexSphere(rp2xs1, 0, "RP^2 x S^1");

            verifyVertexCount(solidKB, 2, "Solid Klein bottle");
            verifyVertexDisc(solidKB, 0, "Solid Klein bottle");
            verifyVertexDisc(solidKB, 1, "Solid Klein bottle");

            verifyVertexCount(gieseking, 1, "Gieseking manifold");
            verifyVertexKB(gieseking, 0, "Gieseking manifold");

            verifyVertexCount(invalidEdges, 2,
                "Triangulation with invalid edges");
            verifyVertexSphere(invalidEdges, 0,
                "Triangulation with invalid edges");
            verifyVertexSphere(invalidEdges, 1,
                "Triangulation with invalid edges");

            verifyVertexCount(twoProjPlaneCusps, 9,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 0,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 1,
                "Triangulation with RP^2 cusps");
            verifyVertexProjPlane(twoProjPlaneCusps, 2,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 3,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 4,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 5,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 6,
                "Triangulation with RP^2 cusps");
            verifyVertexSphere(twoProjPlaneCusps, 7,
                "Triangulation with RP^2 cusps");
            verifyVertexProjPlane(twoProjPlaneCusps, 8,
                "Triangulation with RP^2 cusps");

            verifyVertexCount(cuspedGenusTwoTorus, 2,
                "Cusped solid genus two torus");
            verifyVertexSphere(cuspedGenusTwoTorus, 0,
                "Cusped solid genus two torus");
            verifyVertexTorusG2(cuspedGenusTwoTorus, 1,
                "Cusped solid genus two torus");

            verifyVertexCount(pinchedSolidTorus, 3,
                "Pinched solid torus");
            verifyVertexAnnulus(pinchedSolidTorus, 0,
                "Pinched solid torus");
            verifyVertexDisc(pinchedSolidTorus, 1,
                "Pinched solid torus");
            verifyVertexDisc(pinchedSolidTorus, 2,
                "Pinched solid torus");

            verifyVertexCount(pinchedSolidKB, 3,
                "Pinched solid Klein bottle");
            verifyVertexMobius(pinchedSolidKB, 0,
                "Pinched solid Klein bottle");
            verifyVertexDisc(pinchedSolidKB, 1,
                "Pinched solid Klein bottle");
            verifyVertexDisc(pinchedSolidKB, 2,
                "Pinched solid Klein bottle");
        }

        void verifyEuler(const NTriangulation& tri,
                long expectedManifold, long expectedTri, const char* triName) {
            long eulerManifold = tri.getEulerCharManifold();
            long eulerTri = tri.getEulerCharTri();

            if (eulerManifold != expectedManifold) {
                std::ostringstream msg;
                msg << triName << " gives a manifold Euler characteristic of "
                    << eulerManifold << " instead of "
                    << expectedManifold << ".";
                CPPUNIT_FAIL(msg.str());
            }
            if (eulerTri != expectedTri) {
                std::ostringstream msg;
                msg << triName
                    << " gives a triangulation Euler characteristic of "
                    << eulerTri << " instead of " << expectedTri << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void eulerCharacteristic() {
            verifyEuler(empty, 0, 0, "Empty triangulation");
            verifyEuler(singleTet, 1, 1, "Single tetrahedron");
            verifyEuler(s3, 0, 0, "S^3");
            verifyEuler(s2xs1, 0, 0, "S^2 x S^1");
            verifyEuler(rp3_1, 0, 0, "RP^3 (1 vtx)");
            verifyEuler(rp3_2, 0, 0, "RP^3 (2 vtx)");
            verifyEuler(lens3_1, 0, 0, "L(3,1)");
            verifyEuler(lens8_3, 0, 0, "L(8,3)");
            verifyEuler(lens7_1_loop, 0, 0, "Layered loop L(7,1)");
            verifyEuler(rp3rp3, 0, 0, "RP^3 # RP^3");
            verifyEuler(q32xz3, 0, 0, "S^3 / Q_32 x Z_3");
            verifyEuler(q28, 0, 0, "S^3 / Q_28");
            verifyEuler(weberSeifert, 0, 0, "Weber-Seifert dodecahedral space");
            verifyEuler(lens100_1, 0, 0, "L(100,1)");
            verifyEuler(lst3_4_7, 0, 0, "LST(3,4,7)");
            verifyEuler(figure8, 0, 1, "Figure eight knot complement");
            verifyEuler(rp2xs1, 0, 0, "RP^2 x S^1");
            verifyEuler(solidKB, 0, 0, "Solid Klein bottle");
            verifyEuler(gieseking, 0, 1, "Gieseking manifold");

            verifyEuler(invalidEdges, 1, -1,
                "Triangulation with invalid edges");
            verifyEuler(twoProjPlaneCusps, 1, 1,
                "Triangulation with RP^2 cusps");
            verifyEuler(cuspedGenusTwoTorus, -1, 2,
                "Cusped solid genus two torus");
            verifyEuler(pinchedSolidTorus, 0, 1, "Pinched solid torus");
            verifyEuler(pinchedSolidKB, 0, 1, "Pinched solid Klein bottle");

            {
                NTriangulation t(twoProjPlaneCusps);
                t.idealToFinite();
                verifyEuler(t, 1, 1, "Triangulation with RP^2 boundaries");
            }
            {
                NTriangulation t(cuspedGenusTwoTorus);
                t.idealToFinite();
                verifyEuler(t, -1, -1, "Solid genus two torus");
            }
            {
                NTriangulation t(pinchedSolidTorus);
                t.idealToFinite();
                verifyEuler(t, 0, 0, "Unpinched solid torus");
            }
            {
                NTriangulation t(pinchedSolidKB);
                t.idealToFinite();
                verifyEuler(t, 0, 0, "Unpinched solid Klein bottle");
            }

            verifyEuler(s3_large, 0, 0, "Large S^3");
            verifyEuler(rp3_large, 0, 0, "Large RP^3");
            verifyEuler(lens8_3_large, 0, 0, "Large L(8,3)");
            verifyEuler(q20_large, 0, 0, "Large S^3 / Q_20");
            verifyEuler(ball_large, 1, 1, "4-tetrahedron ball");
            verifyEuler(ball_large_pillows, 1, 1, "4-tetrahedron pillow ball");
            verifyEuler(ball_large_snapped, 1, 1, "3-tetrahedron snapped ball");
            verifyEuler(fig8_bary, 0, 1, "Large figure eight knot complement");
        }

        void verifyGroup(const NAbelianGroup& g, const std::string& grpName,
                unsigned rank) {
            // Construct the error message.
            std::ostringstream msg;
            msg << grpName << " is " << g.toString() << ", not ";
            if (rank == 0)
                msg << "0";
            else if (rank == 1)
                msg << "Z";
            else
                msg << rank << " Z";
            msg << '.';

            // Check the group.
            if (g.getRank() != rank)
                CPPUNIT_FAIL(msg.str());
            else if (g.getNumberOfInvariantFactors() != 0)
                CPPUNIT_FAIL(msg.str());
        }

        void verifyGroup(const NAbelianGroup& g, const std::string& grpName,
                unsigned rank, unsigned long torsionDegree) {
            // Construct the error message.
            std::ostringstream msg;
            msg << grpName << " is " << g.toString() << ", not ";
            if (rank == 1)
                msg << "Z + ";
            else if (rank > 1)
                msg << rank << " Z + ";
            msg << "Z_" << torsionDegree << '.';

            // Check the group.
            if (g.getRank() != rank)
                CPPUNIT_FAIL(msg.str());
            else {
                // Ranks match.
                if (g.getNumberOfInvariantFactors() != 1 ||
                        g.getInvariantFactor(0) != torsionDegree)
                    CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyGroup(const NAbelianGroup& g, const std::string& grpName,
                unsigned rank, unsigned long torsionDegree1,
                unsigned long torsionDegree2) {
            // Construct the error message.
            std::ostringstream msg;
            msg << grpName << " is " << g.toString() << ", not ";
            if (rank == 1)
                msg << "Z + ";
            else if (rank > 1)
                msg << rank << " Z + ";
            msg << "Z_" << torsionDegree1 << " + Z_" << torsionDegree2 << '.';

            // Check the group.
            if (g.getRank() != rank)
                CPPUNIT_FAIL(msg.str());
            else {
                // Ranks match.
                if (g.getNumberOfInvariantFactors() != 2 ||
                        g.getInvariantFactor(0) != torsionDegree1 ||
                        g.getInvariantFactor(1) != torsionDegree2)
                    CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyGroup(const NAbelianGroup& g, const std::string& grpName,
                unsigned rank, unsigned long torsionDegree1,
                unsigned long torsionDegree2, unsigned long torsionDegree3) {
            // Construct the error message.
            std::ostringstream msg;
            msg << grpName << " is " << g.toString() << ", not ";
            if (rank == 1)
                msg << "Z + ";
            else if (rank > 1)
                msg << rank << " Z + ";
            msg << "Z_" << torsionDegree1 << " + Z_" << torsionDegree2
                << " + Z_" << torsionDegree3 << '.';

            // Check the group.
            if (g.getRank() != rank)
                CPPUNIT_FAIL(msg.str());
            else {
                // Ranks match.
                if (g.getNumberOfInvariantFactors() != 3 ||
                        g.getInvariantFactor(0) != torsionDegree1 ||
                        g.getInvariantFactor(1) != torsionDegree2 ||
                        g.getInvariantFactor(2) != torsionDegree3)
                    CPPUNIT_FAIL(msg.str());
            }
        }

        void homologyH1() {
            verifyGroup(empty.getHomologyH1(),
                "H1(empty triangulation)", 0);
            verifyGroup(singleTet.getHomologyH1(),
                "H1(single tetrahedron)", 0);
            verifyGroup(s3.getHomologyH1(),
                "H1(S^3)", 0);
            verifyGroup(s2xs1.getHomologyH1(),
                "H1(S^2 x S^1)", 1);
            verifyGroup(rp3_1.getHomologyH1(),
                "H1(RP^3, 1 vtx)", 0, 2);
            verifyGroup(rp3_2.getHomologyH1(),
                "H1(RP^3, 2 vtx)", 0, 2);
            verifyGroup(lens3_1.getHomologyH1(),
                "H1(L(3,1))", 0, 3);
            verifyGroup(lens7_1_loop.getHomologyH1(),
                "H1(Loop L(7,1))", 0, 7);
            verifyGroup(lens8_3.getHomologyH1(),
                "H1(L(8,3))", 0, 8);
            verifyGroup(lens8_3_large.getHomologyH1(),
                "H1(Large L(8,3))", 0, 8);
            verifyGroup(rp3rp3.getHomologyH1(),
                "H1(RP^3 # RP^3)", 0, 2, 2);
            verifyGroup(q28.getHomologyH1(),
                "H1(S^3 / Q_28)", 0, 4);
            verifyGroup(weberSeifert.getHomologyH1(),
                "H1(SeifertWeber)", 0, 5, 5, 5);
            verifyGroup(q32xz3.getHomologyH1(),
                "H1(S^3 / Q_32 x Z_3)", 0, 2, 6);
            verifyGroup(lens100_1.getHomologyH1(),
                "H1(L(100,1))", 0, 100);
            verifyGroup(ball_large.getHomologyH1(),
                "H1(4-tetrahedron ball)", 0);
            verifyGroup(ball_large_pillows.getHomologyH1(),
                "H1(4-tetrahedron pillow ball)", 0);
            verifyGroup(ball_large_snapped.getHomologyH1(),
                "H1(3-tetrahedron snapped ball)", 0);
            verifyGroup(lst3_4_7.getHomologyH1(),
                "H1(LST(3,4,7))", 1);
            verifyGroup(figure8.getHomologyH1(),
                "H1(figure eight knot complement)", 1);
            verifyGroup(rp2xs1.getHomologyH1(),
                "H1(RP^2 x S^1)", 1, 2);
            verifyGroup(solidKB.getHomologyH1(),
                "H1(solid Klein bottle)", 1);
            verifyGroup(gieseking.getHomologyH1(),
                "H1(Gieseking manifold)", 1);
            verifyGroup(invalidEdges.getHomologyH1(),
                "H1(tri with invalid edges)", 0);
            verifyGroup(twoProjPlaneCusps.getHomologyH1(),
                "H1(tri with projective plane cusps)", 0, 2);
            verifyGroup(cuspedGenusTwoTorus.getHomologyH1(),
                "H1(cusped solid genus two torus)", 2);
            verifyGroup(pinchedSolidTorus.getHomologyH1(),
                "H1(pinched solid torus)", 1);
            verifyGroup(pinchedSolidKB.getHomologyH1(),
                "H1(pinched solid Klein bottle)", 1);
        }

        void homologyH1Bdry() {
            verifyGroup(empty.getHomologyH1Bdry(),
                "Boundary H1(empty triangulation)", 0);
            verifyGroup(singleTet.getHomologyH1Bdry(),
                "Boundary H1(single tetrahedron)", 0);
            verifyGroup(s3.getHomologyH1Bdry(),
                "Boundary H1(S^3)", 0);
            verifyGroup(s2xs1.getHomologyH1Bdry(),
                "Boundary H1(S^2 x S^1)", 0);
            verifyGroup(rp3_1.getHomologyH1Bdry(),
                "Boundary H1(RP^3, 1 vtx)", 0);
            verifyGroup(rp3_2.getHomologyH1Bdry(),
                "Boundary H1(RP^3, 2 vtx)", 0);
            verifyGroup(lens3_1.getHomologyH1Bdry(),
                "Boundary H1(L(3,1))", 0);
            verifyGroup(lens7_1_loop.getHomologyH1Bdry(),
                "Boundary H1(Loop L(7,1))", 0);
            verifyGroup(lens8_3.getHomologyH1Bdry(),
                "Boundary H1(L(8,3))", 0);
            verifyGroup(lens8_3_large.getHomologyH1Bdry(),
                "Boundary H1(Large L(8,3))", 0);
            verifyGroup(rp3rp3.getHomologyH1Bdry(),
                "Boundary H1(RP^3 # RP^3)", 0);
            verifyGroup(q28.getHomologyH1Bdry(),
                "Boundary H1(S^3 / Q_28)", 0);
            verifyGroup(weberSeifert.getHomologyH1Bdry(),
                "Boundary H1(Weber-Seifert)", 0);
            verifyGroup(q32xz3.getHomologyH1Bdry(),
                "Boundary H1(S^3 / Q_32 x Z_3)", 0);
            verifyGroup(lens100_1.getHomologyH1Bdry(),
                "Boundary H1(L(100,1))", 0);
            verifyGroup(ball_large.getHomologyH1Bdry(),
                "Boundary H1(4-tetrahedron ball)", 0);
            verifyGroup(ball_large_pillows.getHomologyH1Bdry(),
                "Boundary H1(4-tetrahedron pillow ball)", 0);
            verifyGroup(ball_large_snapped.getHomologyH1Bdry(),
                "Boundary H1(3-tetrahedron snapped ball)", 0);
            verifyGroup(lst3_4_7.getHomologyH1Bdry(),
                "Boundary H1(LST(3,4,7))", 2);
            verifyGroup(figure8.getHomologyH1Bdry(),
                "Boundary H1(figure eight knot complement)", 2);
            verifyGroup(rp2xs1.getHomologyH1Bdry(),
                "Boundary H1(RP^2 x S^1)", 0);
            verifyGroup(solidKB.getHomologyH1Bdry(),
                "Boundary H1(solid Klein bottle)", 1, 2);
            verifyGroup(gieseking.getHomologyH1Bdry(),
                "Boundary H1(Gieseking manifold)", 1, 2);
            verifyGroup(twoProjPlaneCusps.getHomologyH1Bdry(),
                "Boundary H1(tri with projective plane cusps)", 0, 2, 2);
            verifyGroup(cuspedGenusTwoTorus.getHomologyH1Bdry(),
                "Boundary H1(cusped solid genus two torus)", 4);
        }

        void verifyFundGroup(const NGroupPresentation& g,
                const std::string& grpName, const std::string& expected) {
            std::string actual = g.recogniseGroup();
            if (actual.empty())
                actual = "unknown";

            // Construct the error message.
            std::ostringstream msg;
            msg << grpName << " is " << actual << ", not " << expected << '.';

            // Check the group.
            if (expected != actual)
                CPPUNIT_FAIL(msg.str());
        }

        void fundGroup() {
            verifyFundGroup(empty.getFundamentalGroup(),
                "Fund(empty triangulation)", "0");
            verifyFundGroup(singleTet.getFundamentalGroup(),
                "Fund(single tetrahedron)", "0");
            verifyFundGroup(s3.getFundamentalGroup(),
                "Fund(S^3)", "0");
            verifyFundGroup(s2xs1.getFundamentalGroup(),
                "Fund(S^2 x S^1)", "Z");
            verifyFundGroup(rp3_1.getFundamentalGroup(),
                "Fund(RP^3, 1 vtx)", "Z_2");
            verifyFundGroup(rp3_2.getFundamentalGroup(),
                "Fund(RP^3, 2 vtx)", "Z_2");
            verifyFundGroup(lens3_1.getFundamentalGroup(),
                "Fund(L(3,1))", "Z_3");
            verifyFundGroup(lens7_1_loop.getFundamentalGroup(),
                "Fund(Loop L(7,1))", "Z_7");
            verifyFundGroup(lens8_3.getFundamentalGroup(),
                "Fund(L(8,3))", "Z_8");
            verifyFundGroup(lens8_3_large.getFundamentalGroup(),
                "Fund(Large L(8,3))", "Z_8");
            //verifyFundGroup(rp3rp3.getFundamentalGroup(),
            //    "Fund(RP^3 # RP^3)", 0, 2, 2);
            //verifyFundGroup(q28.getFundamentalGroup(),
            //    "Fund(S^3 / Q_28)", 0, 4);
            //verifyGroup(weberSeifert.getHomologyH1(),
            //    "Fund(SeifertWeber)", 0, 5, 5, 5);
            //verifyFundGroup(q32xz3.getFundamentalGroup(),
            //    "Fund(S^3 / Q_32 x Z_3)", 0, 2, 6);
            verifyFundGroup(lens100_1.getFundamentalGroup(),
                "Fund(L(100,1))", "Z_100");
            verifyFundGroup(ball_large.getFundamentalGroup(),
                "Fund(4-tetrahedron ball)", "0");
            verifyFundGroup(ball_large_pillows.getFundamentalGroup(),
                "Fund(4-tetrahedron pillow ball)", "0");
            verifyFundGroup(ball_large_snapped.getFundamentalGroup(),
                "Fund(3-tetrahedron snapped ball)", "0");
            verifyFundGroup(lst3_4_7.getFundamentalGroup(),
                "Fund(LST(3,4,7))", "Z");
            //verifyFundGroup(figure8.getFundamentalGroup(),
            //    "Fund(figure eight knot complement)", 1);
            //verifyFundGroup(rp2xs1.getFundamentalGroup(),
            //    "Fund(RP^2 x S^1)", 1, 2);
            verifyFundGroup(solidKB.getFundamentalGroup(),
                "Fund(solid Klein bottle)", "Z");
            //verifyFundGroup(gieseking.getFundamentalGroup(),
            //    "Fund(Gieseking manifold)", 1);
            verifyFundGroup(invalidEdges.getFundamentalGroup(),
                "Fund(tri with invalid edges)", "0");
            verifyFundGroup(twoProjPlaneCusps.getFundamentalGroup(),
                "Fund(tri with projective plane cusps)", "Z_2");
            verifyFundGroup(cuspedGenusTwoTorus.getFundamentalGroup(),
                "Fund(cusped solid genus two torus)", "Free (2 generators)");
            verifyFundGroup(pinchedSolidTorus.getFundamentalGroup(),
                "Fund(pinched solid torus)", "Z");
            verifyFundGroup(pinchedSolidKB.getFundamentalGroup(),
                "Fund(pinched solid Klein bottle)", "Z");
        }

        void zeroEfficiency() {
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not "
                "0-efficient.", empty.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron is 0-efficient.",
                ! singleTet.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("S^3 is not 0-efficient.",
                s3.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("S^2 x S^1 is 0-efficient.",
                ! s2xs1.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (1 vtx) is 0-efficient.",
                ! rp3_1.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("RP^3 (2 vtx) is 0-efficient.",
                ! rp3_2.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("L(3,1) is 0-efficient.",
                ! lens3_1.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("Layered loop L(7,1) is 0-efficient.",
                ! lens7_1_loop.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("L(8,3) is not 0-efficient.",
                lens8_3.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("Large L(8,3) is 0-efficient.",
                ! lens8_3_large.isZeroEfficient());
                // Contains a non-vertex-linking 2-sphere.
            CPPUNIT_ASSERT_MESSAGE("RP^3 # RP^3 is 0-efficient.",
                ! rp3rp3.isZeroEfficient());
                // Contains a pair of one-sided projective planes.
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_28 is not 0-efficient.",
                q28.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("S^3 / Q_32 x Z_3 is not 0-efficient.",
                q32xz3.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("The 4-tetrahedron ball is 0-efficient.",
                ! ball_large.isZeroEfficient());
                // Contains a non-trivial disc.
            CPPUNIT_ASSERT_MESSAGE(
                "The 4-tetrahedron pillow ball is 0-efficient.",
                ! ball_large_pillows.isZeroEfficient());
                // Contains a non-trivial disc.
            CPPUNIT_ASSERT_MESSAGE(
                "The 3-tetrahedron snapped ball is 0-efficient.",
                ! ball_large_snapped.isZeroEfficient());
                // Contains a non-trivial disc.
            CPPUNIT_ASSERT_MESSAGE("LST(3,4,7) is 0-efficient.",
                ! lst3_4_7.isZeroEfficient());
                // Contains a non-trivial disc.
            CPPUNIT_ASSERT_MESSAGE("The figure eight knot complement "
                "is not 0-efficient.",
                figure8.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("RP^2 x S^1 is not 0-efficient.",
                rp2xs1.isZeroEfficient());
                // Contains a two-sided projective plane.
            CPPUNIT_ASSERT_MESSAGE("The solid Klein bottle is 0-efficient.",
                ! solidKB.isZeroEfficient());
                // Contains a non-trivial disc.
            CPPUNIT_ASSERT_MESSAGE("The Gieseking manifold is not 0-efficient.",
                gieseking.isZeroEfficient());
            // Don't run this one, it takes too long.
            // CPPUNIT_ASSERT_MESSAGE("The cusped solid genus two torus "
            //     "is 0-efficient.",
            //     ! cuspedGenusTwoTorus.isZeroEfficient());
        }

        void verifyTV3(NTriangulation& t, const std::string& triName) {
            // Verify the Turaev-Viro invariants for r=3.
            // The expected values are described in the paper of Turaev
            // and Viro.
            // For the time being we will only use epsilon == -1 since
            // the expected value is easier to calculate.
            for (unsigned q0 = 2; q0 <= 4; q0++) {
                if (q0 == 3)
                    continue;

                double tv = t.turaevViro(3, q0);
                double expectedTV = 0.5;
                for (unsigned long i = 0; i < t.getHomologyH2Z2(); i++)
                    expectedTV += expectedTV;

                if (regina::isNonZero(tv - expectedTV)) {
                    std::ostringstream msg;
                    msg << "Turaev-Viro(" << triName << ", r = 3, root = "
                        << q0 << ") is " << tv << ", not " << expectedTV
                        << ".";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void verifyTVS3(unsigned r) {
            // Verify all Turaev-Viro invariants for S^3 for the given
            // value of r.
            // The expected values are described in the paper of Turaev
            // and Viro.
            for (unsigned q0 = 1; q0 < 2 * r; q0++) {
                if (regina::gcd(q0, r) > 1)
                    continue;

                double tv = s3.turaevViro(r, q0);

                double expectedTV = 2 * sin(M_PI * q0 / static_cast<double>(r));
                expectedTV = (expectedTV * expectedTV) /
                    (2 * static_cast<double>(r));

                if (regina::isNonZero(tv - expectedTV)) {
                    std::ostringstream msg;
                    msg << "Turaev-Viro(S^3, r = " << r << ", root = " << q0
                        << ") is " << tv << ", not " << expectedTV << ".";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void verifyTVRP3(unsigned r) {
            // Verify all Turaev-Viro invariants for RP^3 for the given
            // value of r.
            // The expected values are described in the paper of Turaev
            // and Viro.
            for (unsigned q0 = 1; q0 < 2 * r; q0++) {
                if (regina::gcd(q0, r) > 1)
                    continue;

                double tv = rp3_2.turaevViro(r, q0);

                double expectedTV;
                if (q0 % 2 == r % 2)
                    expectedTV = 0.0;
                else {
                    double real = cos(M_PI * q0 / static_cast<double>(r)) -
                        1.0;
                    double imag = sin(M_PI * q0 / static_cast<double>(r));
                    expectedTV = (real * real) + (imag * imag);
                    expectedTV /= static_cast<double>(r);
                }

                if (regina::isNonZero(tv - expectedTV)) {
                    std::ostringstream msg;
                    msg << "Turaev-Viro(RP^3, r = " << r << ", root = " << q0
                        << ") is " << tv << ", not " << expectedTV << ".";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void verifyTVL31(unsigned r) {
            // Verify all Turaev-Viro invariants for L(3,1) for the given
            // value of r.
            // The expected values are described in the paper of Turaev
            // and Viro.
            for (unsigned q0 = 1; q0 < 2 * r; q0++) {
                if (regina::gcd(q0, r) > 1)
                    continue;

                double tv = lens3_1.turaevViro(r, q0);

                unsigned q0pow = ((r - 2) / 3) + 1;
                double expectedTV = 2 * sin(M_PI * q0 * q0pow /
                    static_cast<double>(r));
                expectedTV = (expectedTV * expectedTV) /
                    (2 * static_cast<double>(r));

                if (regina::isNonZero(tv - expectedTV)) {
                    std::ostringstream msg;
                    msg << "Turaev-Viro(L(3,1), r = " << r << ", root = "
                        << q0 << ") is " << tv << ", not " << expectedTV
                        << ".";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void verifyTVS2xS1(unsigned r) {
            // Verify all Turaev-Viro invariants for S^2 x S^1 for the given
            // value of r.
            // The expected values are described in the paper of Turaev
            // and Viro.
            for (unsigned q0 = 1; q0 < 2 * r; q0++) {
                if (regina::gcd(q0, r) > 1)
                    continue;

                double tv = s2xs1.turaevViro(r, q0);

                double expectedTV = 1.0;

                if (regina::isNonZero(tv - expectedTV)) {
                    std::ostringstream msg;
                    msg << "Turaev-Viro(S^2 x S^1, r = " << r << ", root = "
                        << q0 << ") is " << tv << ", not " << expectedTV
                        << ".";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void turaevViro() {
            verifyTV3(s3, "S^3");
            verifyTV3(s2xs1, "S^2 x S^1");
            verifyTV3(rp3_1, "RP^3 (1 vtx)");
            verifyTV3(rp3_2, "RP^3 (2 vtx)");
            verifyTV3(lens3_1, "L(3,1)");
            verifyTV3(lens7_1_loop, "Layered loop L(7,1)");
            verifyTV3(lens8_3, "L(8,3)");
            verifyTV3(lens8_3_large, "Large L(8,3)");
            verifyTV3(rp3rp3, "RP^3 # RP^3");
            verifyTV3(q28, "S^3 / Q_28");
            verifyTV3(weberSeifert, "Weber-Seifert");
            verifyTV3(q32xz3, "S^3 / Q_32 x Z_3");
            verifyTV3(rp2xs1, "RP^2 x S^1");

            verifyTVS3(4); verifyTVS3(5); verifyTVS3(6);
            verifyTVS3(7); verifyTVS3(8);

            verifyTVRP3(4); verifyTVRP3(5); verifyTVRP3(6);
            verifyTVRP3(7); verifyTVRP3(8);

            verifyTVL31(4); verifyTVL31(5); verifyTVL31(6);
            verifyTVL31(7); verifyTVL31(8);

            verifyTVS2xS1(4); verifyTVS2xS1(5); verifyTVS2xS1(6);
            verifyTVS2xS1(7); verifyTVS2xS1(8);
        }

        void verifyDoubleCover(const NTriangulation& tri, const char* triName) {
            // PRE: tri is either empty or connected.

            NTriangulation cover(tri);
            cover.makeDoubleCover();

            if (tri.getNumberOfTetrahedra() == 0) {
                if (cover.getNumberOfTetrahedra() != 0)
                    CPPUNIT_FAIL("Empty triangulation: "
                        "Double cover is non-empty.");
                return;
            }

            // We have a non-empty connected triangulation.
            if (tri.isOrientable()) {
                // We should simply come away with two identical copies
                // of tri.
                regina::NContainer parent;
                if (cover.splitIntoComponents(&parent) != 2) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely two components.";
                    CPPUNIT_FAIL(msg.str());
                }

                NTriangulation* child = static_cast<NTriangulation*>(
                    parent.getFirstTreeChild());
                while (child) {
                    if (! tri.isIsomorphicTo(*child).get()) {
                        std::ostringstream msg;
                        msg << triName << ": Orientable double cover "
                            "contains a component not isomorphic to the "
                            "original.";
                        CPPUNIT_FAIL(msg.str());
                    }

                    child = static_cast<NTriangulation*>(
                        child->getNextTreeSibling());
                }
            } else {
                // We should come away with a proper connected double cover.
                if (cover.getNumberOfComponents() != 1) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely one component.";
                    CPPUNIT_FAIL(msg.str());
                }

                if (! cover.isOrientable()) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover is not "
                        "orientable.";
                    CPPUNIT_FAIL(msg.str());
                }

                if (cover.getNumberOfTetrahedra() !=
                        2 * tri.getNumberOfTetrahedra()) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely twice as many tetrahedra.";
                    CPPUNIT_FAIL(msg.str());
                }

                if (cover.getNumberOfFaces() != 2 * tri.getNumberOfFaces()) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely twice as many faces.";
                    CPPUNIT_FAIL(msg.str());
                }

                if (tri.isValid() && cover.getNumberOfEdges() !=
                        2 * tri.getNumberOfEdges()) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely twice as many edges.";
                    CPPUNIT_FAIL(msg.str());
                }

                if (tri.isValid() && (! tri.isIdeal()) &&
                        cover.getNumberOfVertices() !=
                        2 * tri.getNumberOfVertices()) {
                    std::ostringstream msg;
                    msg << triName << ": Orientable double cover does not "
                        "contain precisely twice as many vertices.";
                    CPPUNIT_FAIL(msg.str());
                }

                // We expect the first homology group to be identical,
                // or to be missing a copy of Z_2.
                if (! (tri.getHomologyH1() == cover.getHomologyH1())) {
                    NAbelianGroup hCover(cover.getHomologyH1());
                    hCover.addTorsionElement(2);
                    if (! (tri.getHomologyH1() == hCover)) {
                        std::ostringstream msg;
                        msg << triName << ": Orientable double cover has H1 = "
                            << cover.getHomologyH1().toString()
                            << ", which does not match the original H1 = "
                            << tri.getHomologyH1().toString() << '.';
                        CPPUNIT_FAIL(msg.str());
                    }
                }
            }
        }

        void doubleCover() {
            verifyDoubleCover(empty, "Empty triangulation");
            verifyDoubleCover(singleTet, "Single tetrahedron");
            verifyDoubleCover(s3, "S^3");
            verifyDoubleCover(s2xs1, "S^2 x S^1");
            verifyDoubleCover(rp3_1, "RP^3 (1 vtx)");
            verifyDoubleCover(rp3_2, "RP^3 (2 vtx)");
            verifyDoubleCover(lens3_1, "L(3,1)");
            verifyDoubleCover(lens8_3, "L(8,3)");
            verifyDoubleCover(lens8_3_large, "Large L(8,3)");
            verifyDoubleCover(lens7_1_loop, "Layered loop L(7,1)");
            verifyDoubleCover(rp3rp3, "RP^3 # RP^3");
            verifyDoubleCover(q32xz3, "S^3 / Q_32 x Z_3");
            verifyDoubleCover(q28, "S^3 / Q_28");
            verifyDoubleCover(weberSeifert, "Weber-Seifert");
            verifyDoubleCover(lens100_1, "L(100,1)");
            verifyDoubleCover(ball_large, "4-tetrahedron ball");
            verifyDoubleCover(ball_large_pillows, "4-tetrahedron pillow ball");
            verifyDoubleCover(ball_large_snapped, "3-tetrahedron snapped ball");
            verifyDoubleCover(lst3_4_7, "LST(3,4,7)");
            verifyDoubleCover(figure8, "Figure eight knot complement");
            verifyDoubleCover(rp2xs1, "RP^2 x S^1");
            verifyDoubleCover(solidKB, "Solid Klein bottle");
            verifyDoubleCover(gieseking, "Gieseking manifold");
            verifyDoubleCover(invalidEdges, "Triangulation with invalid edges");
            verifyDoubleCover(twoProjPlaneCusps,
                "Triangulation with RP^2 cusps");
            verifyDoubleCover(cuspedGenusTwoTorus,
                "Cusped solid genus 2 torus");
            verifyDoubleCover(pinchedSolidTorus, "Pinched solid torus");
            verifyDoubleCover(pinchedSolidKB, "Pinched solid Klein bottle");
        }

        void verifyBary(const NTriangulation& tri, const char* name) {
            NTriangulation b(tri);
            b.barycentricSubdivision();

            // Note that subdivisions can turn invalid into valid, but
            // they can never turn valid into invalid.
            if (tri.isValid() && ! b.isValid()) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks validity.";
                CPPUNIT_FAIL(msg.str());
            }

            // Subdivisions can also turn invalid into ideal.
            // Only consider the valid -> valid case here.
            if (tri.isValid() && (tri.isIdeal() != b.isIdeal())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks idealness.";
                CPPUNIT_FAIL(msg.str());
            }

            if (tri.hasBoundaryFaces() != b.hasBoundaryFaces()) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks boundary faces.";
                CPPUNIT_FAIL(msg.str());
            }

            // As with ideal, consider valid inputs only.
            if (tri.isValid() && (tri.isClosed() != b.isClosed())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks closedness.";
                CPPUNIT_FAIL(msg.str());
            }

            if (tri.isOrientable() != b.isOrientable()) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks orientability.";
                CPPUNIT_FAIL(msg.str());
            }

            if (tri.isConnected() != b.isConnected()) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks connectedness.";
                CPPUNIT_FAIL(msg.str());
            }

            if (tri.getNumberOfComponents() != b.getNumberOfComponents()) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks connected components.";
                CPPUNIT_FAIL(msg.str());
            }

            // Invalid vertices can become new boundary components.
            // Only consider valid inputs here.
            if (tri.isValid() && (tri.getNumberOfBoundaryComponents() !=
                    b.getNumberOfBoundaryComponents())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks boundary components.";
                CPPUNIT_FAIL(msg.str());
            }

            // The same problem with invalid triangulations and boundary
            // components bites us with Euler characteristic also.
            if (tri.isValid() &&
                    (tri.getEulerCharTri() != b.getEulerCharTri())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks Euler char (tri).";
                CPPUNIT_FAIL(msg.str());
            }

            if (tri.isValid() &&
                    (tri.getEulerCharManifold() != b.getEulerCharManifold())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks Euler char (mfd).";
                CPPUNIT_FAIL(msg.str());
            }

            // Now run more expensive tests that will be better with
            // *small* triangulations.
            if (! tri.isValid())
                return;

            b.intelligentSimplify();

            if (! (tri.getHomologyH1() == b.getHomologyH1())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks H1.";
                CPPUNIT_FAIL(msg.str());
            }

            if (! (tri.getHomologyH2() == b.getHomologyH2())) {
                std::ostringstream msg;
                msg << name
                    << ": Barycentric subdivision breaks H2.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void barycentricSubdivision() {
            verifyBary(empty, "Empty triangulation");
            verifyBary(singleTet, "Single tetrahedron");
            verifyBary(s3, "S^3");
            verifyBary(s2xs1, "S^2 x S^1");
            verifyBary(rp3_1, "RP^3 (1 vtx)");
            verifyBary(rp3_2, "RP^3 (2 vtx)");
            verifyBary(lens3_1, "L(3,1)");
            verifyBary(lens8_3, "L(8,3)");
            // (too large) verifyBary(lens8_3_large, "Large L(8,3)");
            verifyBary(lens7_1_loop, "Layered loop L(7,1)");
            verifyBary(rp3rp3, "RP^3 # RP^3");
            verifyBary(q32xz3, "S^3 / Q_32 x Z_3");
            verifyBary(q28, "S^3 / Q_28");
            // (too large) verifyBary(weberSeifert, "Weber-Seifert");
            verifyBary(lens100_1, "L(100,1)");
            verifyBary(ball_large, "4-tetrahedron ball");
            verifyBary(ball_large_pillows, "4-tetrahedron pillow ball");
            verifyBary(ball_large_snapped, "3-tetrahedron snapped ball");
            verifyBary(lst3_4_7, "LST(3,4,7)");
            verifyBary(figure8, "Figure eight knot complement");
            verifyBary(rp2xs1, "RP^2 x S^1");
            verifyBary(solidKB, "Solid Klein bottle");
            verifyBary(gieseking, "Gieseking manifold");
            verifyBary(invalidEdges, "Triangulation with invalid edges");
            verifyBary(twoProjPlaneCusps,
                "Triangulation with RP^2 cusps");
            verifyBary(cuspedGenusTwoTorus,
                "Cusped solid genus 2 torus");
            verifyBary(pinchedSolidTorus, "Pinched solid torus");
            verifyBary(pinchedSolidKB, "Pinched solid Klein bottle");
        }

        void verifyIdealToFinite(const NTriangulation& tri,
                const char* triName) {
            NTriangulation finite(tri);
            finite.idealToFinite();

            // Are there any ideal vertices remaining?
            if (finite.isIdeal()) {
                std::ostringstream msg;
                msg << triName << ": idealToFinite() leaves ideal vertices.";
                CPPUNIT_FAIL(msg.str());
            }

            // Are there any invalid vertices remaining?
            for (NTriangulation::VertexIterator vit =
                    finite.getVertices().begin(); vit !=
                    finite.getVertices().end(); ++vit)
                if ((*vit)->isBoundary() && ! (*vit)->isStandard()) {
                    std::ostringstream msg;
                    msg << triName << ": idealToFinite() leaves "
                        "invalid vertices .";
                    CPPUNIT_FAIL(msg.str());
                }

            // Make sure the invalid edges are left alone.
            unsigned oldInvEdges = 0, newInvEdges = 0;
            NTriangulation::EdgeIterator eit;
            for (eit = tri.getEdges().begin();
                    eit != tri.getEdges().end(); ++eit)
                if (! (*eit)->isValid())
                    ++oldInvEdges;
            for (eit = finite.getEdges().begin();
                    eit != finite.getEdges().end(); ++eit)
                if (! (*eit)->isValid())
                    ++newInvEdges;
            if (oldInvEdges != newInvEdges) {
                std::ostringstream msg;
                msg << triName << ": idealToFinite() changes "
                    "invalid edges .";
                CPPUNIT_FAIL(msg.str());
            }

            // Make sure we don't change the number of boundary components.
            if (tri.getNumberOfBoundaryComponents() !=
                    finite.getNumberOfBoundaryComponents()) {
                std::ostringstream msg;
                msg << triName << ": idealToFinite() changes "
                    "the number of boundary components.";
                CPPUNIT_FAIL(msg.str());
            }

            // In the case of a valid triangulation, ensure that the
            // boundary components are topologically unchanged.
            if (tri.isValid()) {
                typedef std::pair<long, bool> BCSpec;
                NTriangulation::BoundaryComponentIterator bcit;

                std::vector<BCSpec> bcOld;
                for (bcit = tri.getBoundaryComponents().begin();
                        bcit != tri.getBoundaryComponents().end(); ++bcit)
                    bcOld.push_back(BCSpec((*bcit)->getEulerCharacteristic(),
                        (*bcit)->isOrientable()));
                std::sort(bcOld.begin(), bcOld.end());

                std::vector<BCSpec> bcNew;
                for (bcit = finite.getBoundaryComponents().begin();
                        bcit != finite.getBoundaryComponents().end(); ++bcit)
                    bcNew.push_back(BCSpec((*bcit)->getEulerCharacteristic(),
                        (*bcit)->isOrientable()));
                std::sort(bcNew.begin(), bcNew.end());

                if (bcOld != bcNew) {
                    std::ostringstream msg;
                    msg << triName << ": idealToFinite() changes "
                        "the topology of one or more boundary components.";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void idealToFinite() {
            verifyIdealToFinite(empty, "Empty triangulation");
            verifyIdealToFinite(singleTet, "Single tetrahedron");
            verifyIdealToFinite(singleTet_bary,
                "Single tetrahedron subdivided");
            verifyIdealToFinite(s3, "S^3");
            verifyIdealToFinite(s2xs1, "S^2 x S^1");
            verifyIdealToFinite(rp3_1, "RP^3 (1 vtx)");
            verifyIdealToFinite(rp3_2, "RP^3 (2 vtx)");
            verifyIdealToFinite(lens3_1, "L(3,1)");
            verifyIdealToFinite(lens8_3, "L(8,3)");
            verifyIdealToFinite(lens8_3_large, "Large L(8,3)");
            verifyIdealToFinite(lens7_1_loop, "Layered loop L(7,1)");
            verifyIdealToFinite(rp3rp3, "RP^3 # RP^3");
            verifyIdealToFinite(q32xz3, "S^3 / Q_32 x Z_3");
            verifyIdealToFinite(q28, "S^3 / Q_28");
            verifyIdealToFinite(weberSeifert, "Weber-Seifert");
            verifyIdealToFinite(lens100_1, "L(100,1)");
            verifyIdealToFinite(ball_large, "4-tetrahedron ball");
            verifyIdealToFinite(ball_large_pillows,
                "4-tetrahedron pillow ball");
            verifyIdealToFinite(ball_large_snapped,
                "3-tetrahedron snapped ball");
            verifyIdealToFinite(lst3_4_7, "LST(3,4,7)");
            verifyIdealToFinite(figure8, "Figure eight knot complement");
            verifyIdealToFinite(fig8_bary,
                "Figure eight knot complement subdivided");
            verifyIdealToFinite(rp2xs1, "RP^2 x S^1");
            verifyIdealToFinite(solidKB, "Solid Klein bottle");
            verifyIdealToFinite(gieseking, "Gieseking manifold");
            verifyIdealToFinite(invalidEdges,
                "Triangulation with invalid edges");
            verifyIdealToFinite(twoProjPlaneCusps,
                "Triangulation with RP^2 cusps");
            verifyIdealToFinite(cuspedGenusTwoTorus,
                "Cusped solid genus 2 torus");
            verifyIdealToFinite(pinchedSolidTorus, "Pinched solid torus");
            verifyIdealToFinite(pinchedSolidKB, "Pinched solid Klein bottle");
        }

        void verifyFiniteToIdeal(const NTriangulation& tri,
                const char* triName) {
            NTriangulation ideal(tri);
            ideal.finiteToIdeal();

            // Are there any boundary faces remaining?
            if (ideal.hasBoundaryFaces()) {
                std::ostringstream msg;
                msg << triName << ": finiteToIdeal() leaves boundary faces.";
                CPPUNIT_FAIL(msg.str());
            }

            // Make sure the invalid edges are left alone.
            unsigned oldInvEdges = 0, newInvEdges = 0;
            NTriangulation::EdgeIterator eit;
            for (eit = tri.getEdges().begin();
                    eit != tri.getEdges().end(); ++eit)
                if (! (*eit)->isValid())
                    ++oldInvEdges;
            for (eit = ideal.getEdges().begin();
                    eit != ideal.getEdges().end(); ++eit)
                if (! (*eit)->isValid())
                    ++newInvEdges;
            if (oldInvEdges != newInvEdges) {
                std::ostringstream msg;
                msg << triName << ": finiteToIdeal() changes "
                    "invalid edges .";
                CPPUNIT_FAIL(msg.str());
            }

            // In the case of a valid triangulation, ensure that the
            // boundary components are topologically unchanged, except
            // for sphere which must vanish.
            if (tri.isValid()) {
                typedef std::pair<long, bool> BCSpec;
                NTriangulation::BoundaryComponentIterator bcit;

                std::vector<BCSpec> bcOld;
                for (bcit = tri.getBoundaryComponents().begin();
                        bcit != tri.getBoundaryComponents().end(); ++bcit)
                    if ((*bcit)->getEulerCharacteristic() != 2)
                        bcOld.push_back(
                            BCSpec((*bcit)->getEulerCharacteristic(),
                            (*bcit)->isOrientable()));
                std::sort(bcOld.begin(), bcOld.end());

                std::vector<BCSpec> bcNew;
                for (bcit = ideal.getBoundaryComponents().begin();
                        bcit != ideal.getBoundaryComponents().end(); ++bcit)
                    bcNew.push_back(BCSpec((*bcit)->getEulerCharacteristic(),
                        (*bcit)->isOrientable()));
                std::sort(bcNew.begin(), bcNew.end());

                if (bcOld != bcNew) {
                    std::ostringstream msg;
                    msg << triName << ": finiteToIdeal() changes "
                        "the topology of one or more non-sphere "
                        "boundary components.";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void finiteToIdeal() {
            verifyFiniteToIdeal(empty, "Empty triangulation");
            verifyFiniteToIdeal(singleTet, "Single tetrahedron");
            verifyFiniteToIdeal(singleTet_bary,
                "Single tetrahedron subdivided");
            verifyFiniteToIdeal(s3, "S^3");
            verifyFiniteToIdeal(s2xs1, "S^2 x S^1");
            verifyFiniteToIdeal(rp3_1, "RP^3 (1 vtx)");
            verifyFiniteToIdeal(rp3_2, "RP^3 (2 vtx)");
            verifyFiniteToIdeal(lens3_1, "L(3,1)");
            verifyFiniteToIdeal(lens8_3, "L(8,3)");
            verifyFiniteToIdeal(lens8_3_large, "Large L(8,3)");
            verifyFiniteToIdeal(lens7_1_loop, "Layered loop L(7,1)");
            verifyFiniteToIdeal(rp3rp3, "RP^3 # RP^3");
            verifyFiniteToIdeal(q32xz3, "S^3 / Q_32 x Z_3");
            verifyFiniteToIdeal(q28, "S^3 / Q_28");
            verifyFiniteToIdeal(weberSeifert, "Weber-Seifert");
            verifyFiniteToIdeal(lens100_1, "L(100,1)");
            verifyFiniteToIdeal(ball_large, "4-tetrahedron ball");
            verifyFiniteToIdeal(ball_large_pillows,
                "4-tetrahedron pillow ball");
            verifyFiniteToIdeal(ball_large_snapped,
                "3-tetrahedron snapped ball");
            verifyFiniteToIdeal(lst3_4_7, "LST(3,4,7)");
            verifyFiniteToIdeal(figure8, "Figure eight knot complement");
            verifyFiniteToIdeal(fig8_bary,
                "Figure eight knot complement subdivided");
            verifyFiniteToIdeal(rp2xs1, "RP^2 x S^1");
            verifyFiniteToIdeal(solidKB, "Solid Klein bottle");
            verifyFiniteToIdeal(gieseking, "Gieseking manifold");
            verifyFiniteToIdeal(invalidEdges,
                "Triangulation with invalid edges");
            verifyFiniteToIdeal(twoProjPlaneCusps,
                "Triangulation with RP^2 cusps");
            verifyFiniteToIdeal(cuspedGenusTwoTorus,
                "Cusped solid genus 2 torus");
            verifyFiniteToIdeal(pinchedSolidTorus, "Pinched solid torus");
            verifyFiniteToIdeal(pinchedSolidKB, "Pinched solid Klein bottle");
        }

        void drillEdge() {
            // Start with the snapped 1-tetrahedron triangulation of the
            // 3-sphere.  Edges 0 and 2 make a Hopf link, and edge 1 is
            // just an interval.
            {
                NTriangulation snap;
                NTetrahedron* tet = snap.newTetrahedron();
                tet->joinTo(0, tet, NPerm4(0, 1));
                tet->joinTo(2, tet, NPerm4(2, 3));

                NTriangulation tmp0(snap);
                tmp0.drillEdge(tmp0.getEdge(0));
                if (! tmp0.isSolidTorus())
                    CPPUNIT_FAIL("Snapped 3-sphere: drilling edge 0 "
                        "does not give a solid torus.");

                NTriangulation tmp1(snap);
                tmp1.drillEdge(tmp1.getEdge(1));
                if (! tmp1.isBall())
                    CPPUNIT_FAIL("Snapped 3-sphere: drilling edge 1 "
                        "does not give a 3-ball.");

                NTriangulation tmp2(snap);
                tmp2.drillEdge(tmp2.getEdge(2));
                if (! tmp2.isSolidTorus())
                    CPPUNIT_FAIL("Snapped 3-sphere: drilling edge 2 "
                        "does not give a solid torus.");
            }

            // Move on to the layered 1-tetrahedron triangulation of the
            // 3-sphere.
            // Edge 0 forms a trefoil, and edge 1 is unknotted.
            {
                NTriangulation layer;
                NTetrahedron* tet = layer.newTetrahedron();
                tet->joinTo(0, tet, NPerm4(1, 2, 3, 0));
                tet->joinTo(2, tet, NPerm4(2, 3));

                NTriangulation tmp0(layer);
                tmp0.drillEdge(tmp0.getEdge(0));
                if (! (tmp0.isValid() && (! tmp0.isIdeal()) &&
                        (tmp0.getNumberOfBoundaryComponents() == 1) &&
                        (tmp0.getHomologyH1().isZ()) &&
                        (tmp0.getBoundaryComponent(0)->isOrientable()) &&
                        (tmp0.getBoundaryComponent(0)->
                            getEulerCharacteristic() == 0) &&
                        (! tmp0.isSolidTorus())))
                    CPPUNIT_FAIL("Layered 3-sphere: drilling edge 0 "
                        "does not give a non-trivial knot complement.");

                NTriangulation tmp1(layer);
                tmp1.drillEdge(tmp1.getEdge(1));
                if (! tmp1.isSolidTorus())
                    CPPUNIT_FAIL("Layered 3-sphere: drilling edge 1 "
                        "does not give a solid torus.");
            }

            // Next try a (1,2,3) layered solid torus.  Every edge is a
            // boundary edge, and so drilling it should leave us with
            // a solid torus again.
            {
                NTriangulation lst;
                lst.insertLayeredSolidTorus(1, 2);

                NTriangulation tmp0(lst);
                tmp0.drillEdge(tmp0.getEdge(0));
                if (! tmp0.isSolidTorus())
                    CPPUNIT_FAIL("LST(1,2,3): drilling edge 0 "
                        "does not give a solid torus.");

                NTriangulation tmp1(lst);
                tmp1.drillEdge(tmp1.getEdge(1));
                if (! tmp1.isSolidTorus())
                    CPPUNIT_FAIL("LST(1,2,3): drilling edge 1 "
                        "does not give a solid torus.");

                NTriangulation tmp2(lst);
                tmp2.drillEdge(tmp2.getEdge(2));
                if (! tmp2.isSolidTorus())
                    CPPUNIT_FAIL("LST(1,2,3): drilling edge 2 "
                        "does not give a solid torus.");
            }

            // Now try a 2-tetrahedron ball, where we drill the internal edge
            // between the two tetrahedra.  The result should be a solid torus.
            {
                NTriangulation ball;
                NTetrahedron* a = ball.newTetrahedron();
                NTetrahedron* b = ball.newTetrahedron();
                a->joinTo(0, b, NPerm4());
                a->joinTo(1, b, NPerm4());

                // The internal edge joins vertices 2-3.
                NTriangulation tmp(ball);
                tmp.drillEdge(tmp.getTetrahedron(0)->getEdge(5));
                if (! tmp.isSolidTorus())
                    CPPUNIT_FAIL("2-tetrahedron ball: drilling the "
                        "internal edge does not give a solid torus.");
            }
        }

        void verifyDehydration(const NTriangulation& tri, const char* name) {
            std::string dehydrate = tri.dehydrate();
            if (dehydrate.empty()) {
                std::ostringstream msg;
                msg << name << ": Cannot dehydrate.";
                CPPUNIT_FAIL(msg.str());
            }

            NTriangulation rehydrate;
            if (! rehydrate.insertRehydration(dehydrate)) {
                std::ostringstream msg;
                msg << name << ": Cannot rehydrate \"" << dehydrate << "\".";
                CPPUNIT_FAIL(msg.str());
            }

            if (! rehydrate.isIsomorphicTo(tri).get()) {
                std::ostringstream msg;
                msg << name << ": Rehydration of \"" << dehydrate
                    << "\" is not isomorphic to the original.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyNoDehydration(const NTriangulation& tri, const char* name) {
            std::string dehydrate = tri.dehydrate();
            if (! dehydrate.empty()) {
                std::ostringstream msg;
                msg << name
                    << ": Should not dehydrate, but instead dehydrates to \""
                    << dehydrate << "\".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void dehydration() {
            verifyDehydration(empty, "Empty triangulation");
            verifyNoDehydration(singleTet, "Single tetrahedron");
            verifyDehydration(s3, "S^3");
            verifyDehydration(s2xs1, "S^2 x S^1");
            verifyDehydration(rp3_1, "RP^3 (1 vtx)");
            verifyDehydration(rp3_2, "RP^3 (2 vtx)");
            verifyDehydration(lens3_1, "L(3,1)");
            verifyDehydration(lens8_3, "L(8,3)");
            verifyDehydration(lens8_3_large, "Large L(8,3)");
            verifyDehydration(lens7_1_loop, "Layered loop L(7,1)");
            verifyDehydration(rp3rp3, "RP^3 # RP^3");
            verifyDehydration(q32xz3, "S^3 / Q_32 x Z_3");
            verifyDehydration(q28, "S^3 / Q_28");
            verifyDehydration(weberSeifert, "Weber-Seifert");
            verifyNoDehydration(lens100_1, "L(100,1)");
            verifyNoDehydration(ball_large, "4-tetrahedron ball");
            verifyNoDehydration(ball_large_pillows,
                "4-tetrahedron pillow ball");
            verifyNoDehydration(ball_large_snapped,
                "4-tetrahedron snapped ball");
            verifyNoDehydration(lst3_4_7, "LST(3,4,7)");
            verifyDehydration(figure8, "Figure eight knot complement");
            verifyDehydration(rp2xs1, "RP^2 x S^1");
            verifyNoDehydration(solidKB, "Solid Klein bottle");
            verifyDehydration(gieseking, "Gieseking manifold");
            verifyDehydration(invalidEdges, "Triangulation with invalid edges");
            verifyDehydration(twoProjPlaneCusps,
                "Triangulation with RP^2 cusps");
            verifyDehydration(cuspedGenusTwoTorus,
                "Cusped solid genus 2 torus");
            verifyNoDehydration(pinchedSolidTorus, "Pinched solid torus");
            verifyNoDehydration(pinchedSolidKB, "Pinched solid Klein bottle");
        }

        void verifyIsoSig(const NTriangulation& tri, const char* name) {
            std::string sig = tri.isoSig();

            if (sig.empty()) {
                std::ostringstream msg;
                msg << name << ": Cannot create isomorphism signature.";
                CPPUNIT_FAIL(msg.str());
            }

            NTriangulation* rebuild = NTriangulation::fromIsoSig(sig);
            if (! rebuild) {
                std::ostringstream msg;
                msg << name << ": Cannot reconstruct from isomorphism "
                    "signature \"" << sig << "\".";
                CPPUNIT_FAIL(msg.str());
            }
            if (! rebuild->isIsomorphicTo(tri).get()) {
                std::ostringstream msg;
                msg << name << ": Reconstruction from \"" << sig
                    << "\" is not isomorphic to the original.";
                CPPUNIT_FAIL(msg.str());
            }
            delete rebuild;

            if (tri.getNumberOfTetrahedra() == 0)
                return;

            std::string otherSig;
            for (unsigned i = 0; i < 10; ++i) {
                NIsomorphism* iso = NIsomorphism::random(
                    tri.getNumberOfTetrahedra());
                NTriangulation* other = iso->apply(&tri);

                otherSig = other->isoSig();
                if (otherSig != sig) {
                    std::ostringstream msg;
                    msg << name << ": Random isomorphism gives different "
                        "signature: " << otherSig << " != " << sig << std::endl;
                    CPPUNIT_FAIL(msg.str());
                }

                delete other;
                delete iso;
            }
            for (unsigned i = 0; i < 10; ++i) {
                NIsomorphism* iso = NIsomorphism::random(
                    tri.getNumberOfTetrahedra());
                NTriangulation* other = new NTriangulation(tri);
                iso->applyInPlace(other);

                otherSig = other->isoSig();
                if (otherSig != sig) {
                    std::ostringstream msg;
                    msg << name << ": Random in-place isomorphism gives "
                        "different signature: "
                        << otherSig << " != " << sig << std::endl;
                    CPPUNIT_FAIL(msg.str());
                }

                delete other;
                delete iso;
            }
        }

        void isomorphismSignature() {
            verifyIsoSig(empty, "Empty triangulation");
            verifyIsoSig(singleTet, "Single tetrahedron");
            verifyIsoSig(s3, "S^3");
            verifyIsoSig(s2xs1, "S^2 x S^1");
            verifyIsoSig(rp3_1, "RP^3 (1 vtx)");
            verifyIsoSig(rp3_2, "RP^3 (2 vtx)");
            verifyIsoSig(lens3_1, "L(3,1)");
            verifyIsoSig(lens8_3, "L(8,3)");
            verifyIsoSig(lens8_3_large, "Large L(8,3)");
            verifyIsoSig(lens7_1_loop, "Layered loop L(7,1)");
            verifyIsoSig(rp3rp3, "RP^3 # RP^3");
            verifyIsoSig(q32xz3, "S^3 / Q_32 x Z_3");
            verifyIsoSig(q28, "S^3 / Q_28");
            verifyIsoSig(weberSeifert, "Weber-Seifert");
            verifyIsoSig(lens100_1, "L(100,1)");
            verifyIsoSig(ball_large, "4-tetrahedron ball");
            verifyIsoSig(ball_large_pillows, "4-tetrahedron pillow ball");
            verifyIsoSig(ball_large_snapped, "4-tetrahedron snapped ball");
            verifyIsoSig(lst3_4_7, "LST(3,4,7)");
            verifyIsoSig(figure8, "Figure eight knot complement");
            verifyIsoSig(rp2xs1, "RP^2 x S^1");
            verifyIsoSig(solidKB, "Solid Klein bottle");
            verifyIsoSig(gieseking, "Gieseking manifold");
            verifyIsoSig(invalidEdges, "Triangulation with invalid edges");
            verifyIsoSig(twoProjPlaneCusps, "Triangulation with RP^2 cusps");
            verifyIsoSig(cuspedGenusTwoTorus, "Cusped solid genus 2 torus");
            verifyIsoSig(pinchedSolidTorus, "Pinched solid torus");
            verifyIsoSig(pinchedSolidKB, "Pinched solid Klein bottle");

            NTriangulation t;
            t.insertTriangulation(lens8_3);
            t.insertTriangulation(ball_large_pillows);
            verifyIsoSig(t, "L(8,3) U B^3");
            t.insertTriangulation(cuspedGenusTwoTorus);
            verifyIsoSig(t, "L(8,3) U B^3 U (cusped genus 2 torus)");
        }

        void verifySimplification(const NTriangulation& tri,
                unsigned simpleSize, const char* simpleName) {
            NTriangulation t(tri);
            t.intelligentSimplify();

            if (t.getNumberOfTetrahedra() != simpleSize) {
                std::ostringstream msg;
                msg << "Large triangulation should simplify to " << simpleName
                    << ", but simplifies to " << t.getNumberOfTetrahedra()
                    << " tetrahedra instead of the expected "
                    << simpleSize << ".";
                CPPUNIT_FAIL(msg.str());
            }

            NStandardTriangulation* std =
                NStandardTriangulation::isStandardTriangulation(&t);
            if (std == 0) {
                std::ostringstream msg;
                msg << "Large triangulation should simplify to " << simpleName
                    << ", but instead simplifies to something unrecognisable.";
                CPPUNIT_FAIL(msg.str());
            }

            if (std->getName() != simpleName) {
                std::ostringstream msg;
                msg << "Large triangulation should simplify to " << simpleName
                    << ", but instead simplifies to " << std->getName() << ".";
                CPPUNIT_FAIL(msg.str());
            }

            // Make sure it does not simplify any further.
            NTriangulation t2(t);
            if (t2.intelligentSimplify()) {
                std::ostringstream msg;
                msg << "The simple triangulation " << std->getName()
                    << " should not simplify any further, but it does.";
                CPPUNIT_FAIL(msg.str());
            }
            if (t2.dumpConstruction() != t.dumpConstruction()) {
                std::ostringstream msg;
                msg << "The simple triangulation " << std->getName()
                    << " should not change at all when simplified again, "
                    "but it does.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void verifyNoSimplification(const NTriangulation& tri,
                const char* name) {
            NTriangulation t(tri);
            if (t.intelligentSimplify()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " simplifies but should not.";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void simplification() {
            verifySimplification(s3_large, 1, "L(1,0)");
            verifySimplification(rp3_large, 2, "L(2,1)");
            verifySimplification(lens8_3_large, 2, "L(8,3)");
            verifySimplification(q20_large, 5, "C~(5)");
            verifySimplification(ball_large, 1, "B3 (3-vtx)");
            verifySimplification(ball_large_pillows, 1, "B3 (3-vtx)");
            verifySimplification(ball_large_snapped, 1, "B3 (3-vtx)");
            verifySimplification(fig8_bary, 2, "SnapPea m004");
            verifySimplification(singleTet_bary, 1, "B3 (3-vtx)");

            // Some triangulations that should not simplify.
            NTriangulation* tri;
            NTetrahedron* tet[4];

            // A triangulation with two degree two projective plane cusps
            // (that should not be simplified away):
            tri = new NTriangulation();
            tri->insertRehydration("cabbbbxww");
            if (tri->getNumberOfTetrahedra() != 2)
                CPPUNIT_FAIL("Custom two-cusped triangulation failed "
                    "to rehydrate.");
            verifyNoSimplification(*tri, "Custom two-cusped triangluation");
            delete tri;

            // A triangulation with an invalid edge that simplifies
            // (where the invalid edge must not be simplified away):
            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[1] = tri->newTetrahedron();
            tet[2] = tri->newTetrahedron();
            tet[3] = tri->newTetrahedron();
            tet[0]->joinTo(3, tet[2], NPerm4());
            tet[0]->joinTo(2, tet[1], NPerm4(2, 3));
            tet[3]->joinTo(3, tet[2], NPerm4(2, 3));
            tet[3]->joinTo(2, tet[1], NPerm4(1, 0));
            if (tri->isValid())
                CPPUNIT_FAIL("Custom invalid triangulation was not built "
                    "properly.");

            tri->intelligentSimplify();
            if (tri->getNumberOfTetrahedra() != 1)
                CPPUNIT_FAIL("Custom invalid triangulation did not simplify "
                    "to 1 tetrahedron.");
            if (tri->isValid() || tri->getEdge(0)->isValid())
                CPPUNIT_FAIL("Custom invalid triangulation did not simplify "
                    "to an invalid triangulation with an invalid edge.");
            delete tri;

            // A solid torus that once upon a time was incorrectly simplified
            // away to a ball.
            tri = new NTriangulation();
            tet[0] = tri->newTetrahedron();
            tet[1] = tri->newTetrahedron();
            tet[2] = tri->newTetrahedron();
            tet[2]->joinTo(3, tet[2], NPerm4(2, 3));
            tet[2]->joinTo(1, tet[1], NPerm4(0, 2, 3, 1));
            tet[2]->joinTo(0, tet[0], NPerm4(3, 0, 1, 2));
            tet[1]->joinTo(3, tet[0], NPerm4(0, 3, 1, 2));
            tet[1]->joinTo(1, tet[0], NPerm4());
            if (tri->getHomologyH1().toString() != "Z")
                CPPUNIT_FAIL("Custom solid torus has incorrect H1.");
            tri->intelligentSimplify();
            if (tri->getHomologyH1().toString() != "Z")
                CPPUNIT_FAIL("Custom solid torus simplifies to "
                    "something different.");
            delete tri;
        }

        void testReordering(const NTriangulation& t, const char* name) {
            NTriangulation a(t);
            a.reorderTetrahedraBFS();

            NTriangulation b(t);
            b.reorderTetrahedraBFS(true);

            NIsomorphism* iso = NIsomorphism::random(t.getNumberOfTetrahedra());
            NTriangulation* c = iso->apply(&t);
            delete iso;

            NTriangulation d(*c);
            d.reorderTetrahedraBFS();

            NTriangulation e(*c);
            e.reorderTetrahedraBFS(true);

            if (! t.isIsomorphicTo(a).get()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " changes its isomorphism class when its tetrahedra "
                    "are reordered in the forward direction.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! t.isIsomorphicTo(b).get()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " changes its isomorphism class when its tetrahedra "
                    "are reordered in the reverse direction.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! t.isIsomorphicTo(*c).get()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " changes its isomorphism class when a random "
                    "isomorphism is applied.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! t.isIsomorphicTo(d).get()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " changes its isomorphism class when a random "
                    "isomorphism is applied and then its tetrahedra are "
                    "reordered in the forward direction.";
                CPPUNIT_FAIL(msg.str());
            }
            if (! t.isIsomorphicTo(e).get()) {
                std::ostringstream msg;
                msg << "Triangulation " << name
                    << " changes its isomorphism class when a random "
                    "isomorphism is applied and then its tetrahedra are "
                    "reordered in the reverse direction.";
                CPPUNIT_FAIL(msg.str());
            }

            delete c;
        }

        void reordering() {
            testReordering(empty, "Empty triangulation");
            testReordering(singleTet, "Single tetrahedron");
            testReordering(s3, "S^3");
            testReordering(s2xs1, "S^2 x S^1");
            testReordering(rp3_1, "RP^3 (1 vtx)");
            testReordering(rp3_2, "RP^3 (2 vtx)");
            testReordering(lens3_1, "L(3,1)");
            testReordering(lens8_3, "L(8,3)");
            testReordering(lens8_3_large, "Large L(8,3)");
            testReordering(lens7_1_loop, "Layered loop L(7,1)");
            testReordering(rp3rp3, "RP^3 # RP^3");
            testReordering(q32xz3, "S^3 / Q_32 x Z_3");
            testReordering(q28, "S^3 / Q_28");
            testReordering(weberSeifert, "Weber-Seifert");
            testReordering(lens100_1, "L(100,1)");
            testReordering(ball_large, "4-tetrahedron ball");
            testReordering(ball_large_pillows, "4-tetrahedron pillow ball");
            testReordering(ball_large_snapped, "3-tetrahedron snapped ball");
            testReordering(lst3_4_7, "LST(3,4,7)");
            testReordering(figure8, "Figure eight knot complement");
            testReordering(rp2xs1, "RP^2 x S^1");
            testReordering(solidKB, "Solid Klein bottle");
            testReordering(gieseking, "Gieseking manifold");
            testReordering(invalidEdges, "Triangulation with invalid edges");
            testReordering(twoProjPlaneCusps, "Triangulation with RP^2 cusps");
            testReordering(cuspedGenusTwoTorus, "Cusped solid genus 2 torus");
            testReordering(pinchedSolidTorus, "Pinched solid torus");
            testReordering(pinchedSolidKB, "Pinched solid Klein bottle");

            // Try this with some disconnected triangulations also.
            {
                NTriangulation t;
                t.insertTriangulation(s2xs1);
                t.insertTriangulation(singleTet);
                t.insertTriangulation(figure8);
                testReordering(t, "(S^2 x S^1) U (Single tetrahedron) U "
                    "(Figure eight knot complement)");
            }
        }

        void propertyUpdates() {
            // Begin with an empty triangulation and calculate various
            // properties.
            NTriangulation t;

            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not valid.",
                t.isValid());
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not orientable.",
                t.isOrientable());
            verifyGroup(t.getHomologyH1(),
                "H1(empty triangulation)", 0);
            verifyGroup(t.getHomologyH1Bdry(),
                "Boundary H1(empty triangulation)", 0);
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation is not "
                "0-efficient.", t.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation has a "
                "splitting surface.", ! t.hasSplittingSurface());
            CPPUNIT_ASSERT_MESSAGE("The empty triangulation has a "
                "two-sphere boundary component.",
                ! t.hasTwoSphereBoundaryComponents());

            // Add a single tetrahedron.
            t.newTetrahedron();

            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron is "
                "0-efficient.", ! t.isZeroEfficient());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron has no "
                "splitting surface.", t.hasSplittingSurface());
            CPPUNIT_ASSERT_MESSAGE("A single tetrahedron has no "
                "two-sphere boundary components.",
                t.hasTwoSphereBoundaryComponents());

            // Glue the tetrahedron to itself to form a solid torus.
            t.getTetrahedron(0)->joinTo(0, t.getTetrahedron(0),
                NPerm4(1, 2, 3, 0));

            verifyGroup(t.getHomologyH1(),
                "H1(LST(1,2,3))", 1);
            verifyGroup(t.getHomologyH1Bdry(),
                "Boundary H1(LST(1,2,3))", 2);

            // Glue the remaining two faces in a non-orientable fashion.
            t.getTetrahedron(0)->joinTo(2, t.getTetrahedron(0),
                NPerm4(1, 0, 3, 2));

            CPPUNIT_ASSERT_MESSAGE("A bad 1-tetrahedron triangulation "
                "is valid.", ! t.isValid());
            CPPUNIT_ASSERT_MESSAGE("A bad 1-tetrahedron triangulation "
                "is orientable.", ! t.isOrientable());
        }
};

void addNTriangulation(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(NTriangulationTest::suite());
}

