
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

#include "manifold/nsfs.h"
#include "subcomplex/nsatblocktypes.h"
#include "subcomplex/nlayeredsolidtorus.h"
#include "triangulation/nedge.h"
#include "triangulation/nfacepair.h"
#include "triangulation/ntetrahedron.h"
#include "triangulation/ntriangulation.h"
#include <algorithm>
#include <cstdlib> // For exit().
#include <iterator>
#include <list>

namespace regina {

bool NSatBlock::operator < (const NSatBlock& compare) const {
    const NSatTriPrism* tri1 = dynamic_cast<const NSatTriPrism*>(this);
    const NSatTriPrism* tri2 = dynamic_cast<const NSatTriPrism*>(&compare);
    if (tri1 && ! tri2)
        return true;
    if (tri2 && ! tri1)
        return false;
    if (tri1 && tri2) {
        // Major first, then minor.
        return (tri1->isMajor() && ! tri2->isMajor());
    }

    const NSatCube* cube1 = dynamic_cast<const NSatCube*>(this);
    const NSatCube* cube2 = dynamic_cast<const NSatCube*>(&compare);
    if (cube1 && ! cube2)
        return true;
    if (cube2 && ! cube1)
        return false;
    if (cube1 && cube2) {
        // All cubes are considered equal.
        return false;
    }

    const NSatReflectorStrip* ref1 =
        dynamic_cast<const NSatReflectorStrip*>(this);
    const NSatReflectorStrip* ref2 =
        dynamic_cast<const NSatReflectorStrip*>(&compare);
    if (ref1 && ! ref2)
        return true;
    if (ref2 && ! ref1)
        return false;
    if (ref1 && ref2) {
        // Always put untwisted before twisted.
        if (ref1->twistedBoundary() && ! ref2->twistedBoundary())
            return false;
        if (ref2->twistedBoundary() && ! ref1->twistedBoundary())
            return true;
        return (ref1->nAnnuli() < ref2->nAnnuli());
    }

    const NSatLST* lst1 = dynamic_cast<const NSatLST*>(this);
    const NSatLST* lst2 = dynamic_cast<const NSatLST*>(&compare);
    if (lst1 && ! lst2)
        return true;
    if (lst2 && ! lst1)
        return false;
    if (lst1 && lst2) {
        // Order first by LST parameters, then by roles.
        if (lst1->lst()->getMeridinalCuts(2) < lst2->lst()->getMeridinalCuts(2))
            return true;
        if (lst1->lst()->getMeridinalCuts(2) > lst2->lst()->getMeridinalCuts(2))
            return false;
        if (lst1->lst()->getMeridinalCuts(1) < lst2->lst()->getMeridinalCuts(1))
            return true;
        if (lst1->lst()->getMeridinalCuts(1) > lst2->lst()->getMeridinalCuts(1))
            return false;
        if (lst1->lst()->getMeridinalCuts(0) < lst2->lst()->getMeridinalCuts(0))
            return true;
        if (lst1->lst()->getMeridinalCuts(0) > lst2->lst()->getMeridinalCuts(0))
            return false;

        // Sorts by which edge group is joined to the vertical annulus
        // edges, then horizontal, then diagonal (though we won't bother
        // testing diagonal, since by that stage we will know the roles
        // permutations to be equal).
        if (lst1->roles()[0] < lst2->roles()[0])
            return true;
        if (lst1->roles()[0] > lst2->roles()[0])
            return false;
        if (lst1->roles()[1] < lst2->roles()[1])
            return true;
        if (lst1->roles()[1] > lst2->roles()[1])
            return false;

        // All equal.
        return false;
    }

    const NSatMobius* mob1 = dynamic_cast<const NSatMobius*>(this);
    const NSatMobius* mob2 = dynamic_cast<const NSatMobius*>(&compare);
    if (mob1 && ! mob2)
        return true;
    if (mob2 && ! mob1)
        return false;
    if (mob1 && mob2) {
        // Order by position in descending order (vertical first, then
        // horizontal, then finally diagonal).
        return (mob1->position() > mob2->position());
    }

    const NSatLayering* layer1 = dynamic_cast<const NSatLayering*>(this);
    const NSatLayering* layer2 = dynamic_cast<const NSatLayering*>(&compare);
    if (layer1 && ! layer2)
        return true;
    if (layer2 && ! layer1)
        return false;
    if (layer1 && layer2) {
        // Horizontal, then diagonal.
        return (layer1->overHorizontal() && ! layer2->overHorizontal());
    }

    return false;
}

NSatBlock* NSatBlock::isBlock(const NSatAnnulus& annulus, TetList& avoidTets) {
    NSatBlock* ans;

    // Run through the types of blocks that we know about.
    if ((ans = NSatMobius::isBlockMobius(annulus, avoidTets)))
        return ans;
    if ((ans = NSatLST::isBlockLST(annulus, avoidTets)))
        return ans;
    if ((ans = NSatTriPrism::isBlockTriPrism(annulus, avoidTets)))
        return ans;
    if ((ans = NSatCube::isBlockCube(annulus, avoidTets)))
        return ans;
    if ((ans = NSatReflectorStrip::isBlockReflectorStrip(annulus, avoidTets)))
        return ans;

    // As a last attempt, try a single layering.  We don't have to worry
    // about the degeneracy, since we'll never get a loop of these
    // things (since that would form a disconnected component, and we
    // never use one as a starting block).
    if ((ans = NSatLayering::isBlockLayering(annulus, avoidTets)))
        return ans;

    // Nothing was found.
    return 0;
}

void NSatMobius::adjustSFS(NSFSpace& sfs, bool reflect) const {
    if (position_ == 0) {
        // Diagonal:
        sfs.insertFibre(1, reflect ? 1 : -1);
    } else if (position_ == 1) {
        // Horizontal:
        sfs.insertFibre(1, reflect ? -2 : 2);
    } else {
        // Vertical:
        sfs.insertFibre(2, reflect ? -1 : 1);
    }
}

void NSatMobius::writeTextShort(std::ostream& out) const {
    out << "Saturated Mobius band, boundary on ";
    if (position_ == 0)
        out << "diagonal";
    else if (position_ == 1)
        out << "horizontal";
    else if (position_ == 2)
        out << "vertical";
    else
        out << "invalid";
    out << " edge";
}

void NSatMobius::writeAbbr(std::ostream& out, bool tex) const {
    out << (tex ? "M_" : "Mob(");
    if (position_ == 0)
        out << 'd';
    else if (position_ == 1)
        out << 'h';
    else if (position_ == 2)
        out << 'v';
    if (! tex)
        out << ')';
}

NSatMobius* NSatMobius::isBlockMobius(const NSatAnnulus& annulus, TetList&) {
    // The two tetrahedra must be joined together along the annulus faces.

    if (annulus.tet[0]->adjacentTetrahedron(annulus.roles[0][3]) !=
            annulus.tet[1])
        return 0;

    NPerm4 annulusGluing = annulus.roles[1].inverse() *
        annulus.tet[0]->adjacentGluing(annulus.roles[0][3]) *
        annulus.roles[0];

    if (annulusGluing[3] != 3)
        return 0;

    // The faces are glued together.  Is it one of the allowable
    // (orientable) permutations?

    int position = -1;
    if (annulusGluing == NPerm4(0, 1))
        position = 2; // Vertical
    else if (annulusGluing == NPerm4(0, 2))
        position = 1; // Horizontal
    else if (annulusGluing == NPerm4(1, 2))
        position = 0; // Diagonal

    if (position < 0) {
        // Nope.  It must be a non-orientable permutation.
        return 0;
    }

    // Got it!

    NSatMobius* ans = new NSatMobius(position);
    ans->annulus_[0] = annulus;
    return ans;
}

NSatLST::NSatLST(const NSatLST& cloneMe) : NSatBlock(cloneMe),
        lst_(cloneMe.lst_->clone()), roles_(cloneMe.roles_) {
}

NSatLST::~NSatLST() {
    delete lst_;
}

void NSatLST::adjustSFS(NSFSpace& sfs, bool reflect) const {
    long cutsVert = lst_->getMeridinalCuts(roles_[0]);
    long cutsHoriz = lst_->getMeridinalCuts(roles_[1]);
    if (roles_[2] == 2) {
        // Most cuts are on the diagonal, which means the meridinal
        // curve is negative.
        cutsHoriz = -cutsHoriz;
    }

    sfs.insertFibre(cutsVert, reflect ? -cutsHoriz : cutsHoriz);
}

void NSatLST::writeTextShort(std::ostream& out) const {
    out << "Saturated ("
        << lst_->getMeridinalCuts(0) << ", "
        << lst_->getMeridinalCuts(1) << ", "
        << lst_->getMeridinalCuts(2) << ") layered solid torus";
}

void NSatLST::writeAbbr(std::ostream& out, bool tex) const {
    out << (tex ? "\\mathrm{LST}_{" : "LST(")
        << lst_->getMeridinalCuts(0) << ", "
        << lst_->getMeridinalCuts(1) << ", "
        << lst_->getMeridinalCuts(2) << (tex ? '}' : ')');
}

void NSatLST::transform(const NTriangulation* originalTri,
        const NIsomorphism* iso, NTriangulation* newTri) {
    // Start with the parent implementation.
    NSatBlock::transform(originalTri, iso, newTri);

    // Transform the layered solid torus also.
    lst_->transform(originalTri, iso, newTri);
}

NSatLST* NSatLST::isBlockLST(const NSatAnnulus& annulus, TetList& avoidTets) {
    // Do we move to a common usable tetrahedron?

    if (annulus.tet[0] != annulus.tet[1])
        return 0;
    if (isBad(annulus.tet[0], avoidTets))
        return 0;

    // Is it a layering?

    // Here we find the endpoints of the edge from which the two layered
    // faces fold out.
    NFacePair centralEdge =
        NFacePair(annulus.roles[0][3], annulus.roles[1][3]).complement();

    if (annulus.roles[1] !=
            NPerm4(annulus.roles[0][3], annulus.roles[1][3]) *
            NPerm4(centralEdge.upper(), centralEdge.lower()) *
            annulus.roles[0])
        return 0;

    // Find the layered solid torus.
    NLayeredSolidTorus* lst = NLayeredSolidTorus::formsLayeredSolidTorusTop(
        annulus.tet[0], annulus.roles[0][3], annulus.roles[1][3]);
    if (! lst)
        return 0;

    // Make sure we're not about to create a (0,k) curve.
    NPerm4 lstRoles(
        lst->getTopEdgeGroup(
            NEdge::edgeNumber[annulus.roles[0][0]][annulus.roles[0][1]]),
        lst->getTopEdgeGroup(
            NEdge::edgeNumber[annulus.roles[0][0]][annulus.roles[0][2]]),
        lst->getTopEdgeGroup(
            NEdge::edgeNumber[annulus.roles[0][1]][annulus.roles[0][2]]),
        3);

    if (lst->getMeridinalCuts(lstRoles[0]) == 0)
        return 0;

    // Make two runs through the full set of tetrahedra.
    // The first run verifies that each tetrahedron is usable.
    // The second run inserts the tetrahedra into avoidTets.
    NTetrahedron* current = annulus.tet[0];
    NFacePair currPair = centralEdge;
    NFacePair nextPair;
    while (current != lst->getBase()) {
        // INV: The current tetrahedron is usable.
        // INV: The next two faces to push through are in currPair.

        // Push through to the next tetrahedron.
        nextPair = NFacePair(
            current->adjacentFace(currPair.upper()),
            current->adjacentFace(currPair.lower())
            ).complement();
        current = current->adjacentTetrahedron(currPair.upper());
        currPair = nextPair;

        // Make sure this next tetrahedron is usable.
        if (isBad(current, avoidTets))
            return 0;
    }

    // All good!
    current = annulus.tet[0];
    currPair = centralEdge;
    avoidTets.insert(current);
    while (current != lst->getBase()) {
        // INV: All tetrahedra up to and including current have been added.
        // INV: The next two faces to push through are in currPair.

        // Push through to the next tetrahedron.
        nextPair = NFacePair(
            current->adjacentFace(currPair.upper()),
            current->adjacentFace(currPair.lower())
            ).complement();
        current = current->adjacentTetrahedron(currPair.upper());
        currPair = nextPair;

        // Add this next tetrahedron to the list.
        avoidTets.insert(current);
    }

    NSatLST* ans = new NSatLST(lst, lstRoles);
    ans->annulus_[0] = annulus;
    return ans;
}

void NSatTriPrism::adjustSFS(NSFSpace& sfs, bool reflect) const {
    if (major_)
        sfs.insertFibre(1, reflect ? -1 : 1);
    else
        sfs.insertFibre(1, reflect ? -2 : 2);
}

NSatTriPrism* NSatTriPrism::isBlockTriPrism(const NSatAnnulus& annulus,
        TetList& avoidTets) {
    NSatTriPrism* ans;

    // First try for one of major type.
    if ((ans = isBlockTriPrismMajor(annulus, avoidTets)))
        return ans;

    // Now try the reflected version.
    NSatAnnulus altAnnulus = annulus.verticalReflection();
    if ((ans = isBlockTriPrismMajor(altAnnulus, avoidTets))) {
        // Reflect it back again but mark it as a minor variant.
        ans->major_ = false;

        ans->annulus_[0].reflectVertical();
        ans->annulus_[1].reflectVertical();
        ans->annulus_[2].reflectVertical();

        return ans;
    }

    // Neither variant was found.
    return 0;
}

NSatTriPrism* NSatTriPrism::isBlockTriPrismMajor(const NSatAnnulus& annulus,
        TetList& avoidTets) {
    if (annulus.tet[0] == annulus.tet[1])
        return 0;
    if (isBad(annulus.tet[0], avoidTets) || isBad(annulus.tet[1], avoidTets))
        return 0;
    if (annulus.tet[0]->adjacentTetrahedron(annulus.roles[0][0]) !=
            annulus.tet[1])
        return 0;
    if (annulus.tet[0]->adjacentGluing(annulus.roles[0][0]) *
            annulus.roles[0] * NPerm4(1, 2) != annulus.roles[1])
        return 0;

    // The two tetrahedra forming the annulus are joined together as
    // expected.  Look for the third tetrahedron.

    NTetrahedron* adj = annulus.tet[0]->adjacentTetrahedron(
        annulus.roles[0][1]);
    if (adj == 0 || adj == annulus.tet[0] || adj == annulus.tet[1])
        return 0;
    if (isBad(adj, avoidTets))
        return 0;

    NPerm4 adjRoles =
        annulus.tet[0]->adjacentGluing(annulus.roles[0][1]) *
        annulus.roles[0] * NPerm4(0, 3);

    if (annulus.tet[1]->adjacentTetrahedron(annulus.roles[1][1]) != adj)
        return 0;
    if (annulus.tet[1]->adjacentGluing(annulus.roles[1][1]) *
            annulus.roles[1] * NPerm4(1, 3, 0, 2) != adjRoles)
        return 0;

    // All three tetrahedra are joined together as expected!

    NSatTriPrism* ans = new NSatTriPrism(true);

    const NPerm4 pairSwap(1, 0, 3, 2);
    ans->annulus_[0] = annulus;
    ans->annulus_[1].tet[0] = annulus.tet[1];
    ans->annulus_[1].tet[1] = adj;
    ans->annulus_[1].roles[0] = annulus.roles[1] * pairSwap;
    ans->annulus_[1].roles[1] = adjRoles;
    ans->annulus_[2].tet[0] = adj;
    ans->annulus_[2].tet[1] = annulus.tet[0];
    ans->annulus_[2].roles[0] = adjRoles * pairSwap;
    ans->annulus_[2].roles[1] = annulus.roles[0] * pairSwap;

    avoidTets.insert(annulus.tet[0]);
    avoidTets.insert(annulus.tet[1]);
    avoidTets.insert(adj);

    return ans;
}

NSatTriPrism* NSatTriPrism::insertBlock(NTriangulation& tri, bool major) {
    NTetrahedron* a = tri.newTetrahedron();
    NTetrahedron* b = tri.newTetrahedron();
    NTetrahedron* c = tri.newTetrahedron();
    a->joinTo(1, c, NPerm4(2, 0, 3, 1));
    b->joinTo(1, a, NPerm4(2, 0, 3, 1));
    c->joinTo(1, b, NPerm4(2, 0, 3, 1));

    NSatTriPrism* ans = new NSatTriPrism(major);

    const NPerm4 id;
    const NPerm4 pairSwap(1, 0, 3, 2);
    ans->annulus_[0].tet[0] = a;
    ans->annulus_[0].tet[1] = b;
    ans->annulus_[0].roles[0] = id;
    ans->annulus_[0].roles[1] = pairSwap;
    ans->annulus_[1].tet[0] = b;
    ans->annulus_[1].tet[1] = c;
    ans->annulus_[1].roles[0] = id;
    ans->annulus_[1].roles[1] = pairSwap;
    ans->annulus_[2].tet[0] = c;
    ans->annulus_[2].tet[1] = a;
    ans->annulus_[2].roles[0] = id;
    ans->annulus_[2].roles[1] = pairSwap;

    if (! major) {
        ans->annulus_[0].reflectVertical();
        ans->annulus_[1].reflectVertical();
        ans->annulus_[2].reflectVertical();
    }

    return ans;
}

void NSatCube::adjustSFS(NSFSpace& sfs, bool reflect) const {
    sfs.insertFibre(1, reflect ? -2 : 2);
}

NSatCube* NSatCube::isBlockCube(const NSatAnnulus& annulus,
        TetList& avoidTets) {
    if (annulus.tet[0] == annulus.tet[1])
        return 0;
    if (isBad(annulus.tet[0], avoidTets) || isBad(annulus.tet[1], avoidTets))
        return 0;

    NTetrahedron* central0 = annulus.tet[0]->adjacentTetrahedron(
        annulus.roles[0][0]);
    NTetrahedron* central1 = annulus.tet[0]->adjacentTetrahedron(
        annulus.roles[0][1]);

    if (central0 == 0 || central0 == annulus.tet[0] ||
            central0 == annulus.tet[1] || isBad(central0, avoidTets))
        return 0;
    if (central1 == 0 || central1 == annulus.tet[0] ||
            central1 == annulus.tet[1] || central1 == central0 ||
            isBad(central0, avoidTets))
        return 0;

    NPerm4 roles0 = annulus.tet[0]->adjacentGluing(
        annulus.roles[0][0]) * annulus.roles[0];
    NPerm4 roles1 = annulus.tet[0]->adjacentGluing(
        annulus.roles[0][1]) * annulus.roles[0];

    // We've got the two central tetrahedra.  Now look for the remaining
    // three boundary tetrahedra.

    if (annulus.tet[1]->adjacentTetrahedron(annulus.roles[1][0]) !=
            central0)
        return 0;
    if (annulus.tet[1]->adjacentTetrahedron(annulus.roles[1][1]) !=
            central1)
        return 0;
    if (annulus.tet[1]->adjacentGluing(annulus.roles[1][0]) *
            annulus.roles[1] * NPerm4(3, 2, 1, 0) != roles0)
        return 0;
    if (annulus.tet[1]->adjacentGluing(annulus.roles[1][1]) *
            annulus.roles[1] * NPerm4(2, 3, 0, 1) != roles1)
        return 0;

    // We've got the two tetrahedra from the annulus boundary completely
    // sorted out.  Just the two new boundary tetrahedra to go.

    NTetrahedron* bdry2 = central0->adjacentTetrahedron(roles0[1]);
    NPerm4 roles2 = central0->adjacentGluing(roles0[1]) * roles0;

    NTetrahedron* bdry3 = central0->adjacentTetrahedron(roles0[2]);
    NPerm4 roles3 = central0->adjacentGluing(roles0[2]) * roles0;

    if (bdry2 == 0 || bdry2 == annulus.tet[0] || bdry2 == annulus.tet[1] ||
            bdry2 == central0 || bdry2 == central1 ||
            isBad(bdry2, avoidTets))
        return 0;
    if (bdry3 == 0 || bdry3 == annulus.tet[0] || bdry3 == annulus.tet[1] ||
            bdry3 == central0 || bdry3 == central1 || bdry3 == bdry2 ||
            isBad(bdry3, avoidTets))
        return 0;
    if (central1->adjacentTetrahedron(roles1[0]) != bdry2)
        return 0;
    if (central1->adjacentTetrahedron(roles1[2]) != bdry3)
        return 0;
    if (central1->adjacentGluing(roles1[0]) * roles1 != roles2)
        return 0;
    if (central1->adjacentGluing(roles1[2]) * roles1 *
            NPerm4(1, 0, 3, 2) != roles3)
        return 0;

    // All looking good!

    NSatCube* ans = new NSatCube();

    const NPerm4 id;
    ans->annulus_[0] = annulus;
    ans->annulus_[1].tet[0] = annulus.tet[1];
    ans->annulus_[1].tet[1] = bdry2;
    ans->annulus_[1].roles[0] = annulus.roles[1] * NPerm4(1, 0, 3, 2);
    ans->annulus_[1].roles[1] = roles2;
    ans->annulus_[2].tet[0] = bdry2;
    ans->annulus_[2].tet[1] = bdry3;
    ans->annulus_[2].roles[0] = roles2 * NPerm4(1, 0, 3, 2);
    ans->annulus_[2].roles[1] = roles3 * NPerm4(2, 3, 0, 1);
    ans->annulus_[3].tet[0] = bdry3;
    ans->annulus_[3].tet[1] = annulus.tet[0];
    ans->annulus_[3].roles[0] = roles3 * NPerm4(3, 2, 1, 0);
    ans->annulus_[3].roles[1] = annulus.roles[0] * NPerm4(1, 0, 3, 2);

    avoidTets.insert(annulus.tet[0]);
    avoidTets.insert(annulus.tet[1]);
    avoidTets.insert(bdry2);
    avoidTets.insert(bdry3);
    avoidTets.insert(central0);
    avoidTets.insert(central1);

    return ans;
}

NSatCube* NSatCube::insertBlock(NTriangulation& tri) {
    NTetrahedron* bdry0 = tri.newTetrahedron();
    NTetrahedron* bdry1 = tri.newTetrahedron();
    NTetrahedron* bdry2 = tri.newTetrahedron();
    NTetrahedron* bdry3 = tri.newTetrahedron();
    NTetrahedron* central0 = tri.newTetrahedron();
    NTetrahedron* central1 = tri.newTetrahedron();

    const NPerm4 id;
    bdry0->joinTo(1, central0, id);
    bdry0->joinTo(0, central1, NPerm4(0, 1));
    bdry1->joinTo(2, central0, NPerm4(2, 1, 3, 0));
    bdry1->joinTo(0, central1, NPerm4(0, 3));
    bdry2->joinTo(0, central0, id);
    bdry2->joinTo(1, central1, NPerm4(0, 1));
    bdry3->joinTo(3, central0, NPerm4(0, 3, 1, 2));
    bdry3->joinTo(1, central1, NPerm4(1, 2));

    NSatCube* ans = new NSatCube();

    ans->annulus_[0].tet[0] = bdry0;
    ans->annulus_[0].tet[1] = bdry1;
    ans->annulus_[1].tet[0] = bdry1;
    ans->annulus_[1].tet[1] = bdry2;
    ans->annulus_[2].tet[0] = bdry2;
    ans->annulus_[2].tet[1] = bdry3;
    ans->annulus_[3].tet[0] = bdry3;
    ans->annulus_[3].tet[1] = bdry0;

    ans->annulus_[0].roles[0] = NPerm4(0, 1);
    ans->annulus_[0].roles[1] = NPerm4(2, 0, 3, 1);
    ans->annulus_[1].roles[0] = NPerm4(1, 2);
    ans->annulus_[1].roles[1] = NPerm4(0, 1);
    ans->annulus_[2].roles[0] = NPerm4(2, 3);
    ans->annulus_[2].roles[1] = NPerm4(0, 3);
    ans->annulus_[3].roles[0] = NPerm4(1, 3, 0, 2);
    ans->annulus_[3].roles[1] = NPerm4(2, 3);

    return ans;
}

void NSatReflectorStrip::adjustSFS(NSFSpace& sfs, bool) const {
    if (! twistedBoundary_)
        sfs.addReflector(false);
}

NSatReflectorStrip* NSatReflectorStrip::isBlockReflectorStrip(
        const NSatAnnulus& annulus, TetList& avoidTets) {
    // Hunt for the initial segment of the reflector strip that lies
    // behind the given annulus.

    if (annulus.tet[0] == annulus.tet[1])
        return 0;
    if (isBad(annulus.tet[0], avoidTets) || isBad(annulus.tet[1], avoidTets))
        return 0;

    NTetrahedron* middle = annulus.tet[0]->adjacentTetrahedron(
        annulus.roles[0][0]);
    NPerm4 middleRoles = annulus.tet[0]->adjacentGluing(
        annulus.roles[0][0]) * annulus.roles[0] * NPerm4(3, 1, 0, 2);

    if (notUnique(middle, annulus.tet[0], annulus.tet[1]) ||
            isBad(middle, avoidTets))

    if (middle != annulus.tet[0]->adjacentTetrahedron(
            annulus.roles[0][1]))
        return 0;
    if (middle != annulus.tet[1]->adjacentTetrahedron(
            annulus.roles[1][0]))
        return 0;
    if (middle != annulus.tet[1]->adjacentTetrahedron(
            annulus.roles[1][1]))
        return 0;
    if (middleRoles != annulus.tet[0]->adjacentGluing(
            annulus.roles[0][1]) * annulus.roles[0] * NPerm4(1, 3))
        return 0;
    if (middleRoles != annulus.tet[1]->adjacentGluing(
            annulus.roles[1][0]) * annulus.roles[1] * NPerm4(0, 2, 3, 1))
        return 0;
    if (middleRoles != annulus.tet[1]->adjacentGluing(
            annulus.roles[1][1]) * annulus.roles[1] * NPerm4(0, 2))
        return 0;

    // We've found the initial segment.
    // Do we just have a segment of length one?
    if (annulus.tet[0]->adjacentTetrahedron(annulus.roles[0][2]) ==
            annulus.tet[1]) {
        // It's either length one or nothing.
        if (annulus.roles[1] == annulus.tet[0]->adjacentGluing(
                annulus.roles[0][2]) * annulus.roles[0] * NPerm4(0, 1)) {
            // Got one that's untwisted.
            NSatReflectorStrip* ans = new NSatReflectorStrip(1, false);
            ans->annulus_[0] = annulus;

            avoidTets.insert(annulus.tet[0]);
            avoidTets.insert(middle);
            avoidTets.insert(annulus.tet[1]);

            return ans;
        }

        if (annulus.roles[1] == annulus.tet[0]->adjacentGluing(
                annulus.roles[0][2]) * annulus.roles[0]) {
            // Got one that's twisted.
            NSatReflectorStrip* ans = new NSatReflectorStrip(1, true);
            ans->annulus_[0] = annulus;

            avoidTets.insert(annulus.tet[0]);
            avoidTets.insert(middle);
            avoidTets.insert(annulus.tet[1]);

            return ans;
        }
        // Nup.  Nothing.
        return 0;
    }

    // If anything, we have a segment of length >= 2.  Start following
    // it around.

    // Make a list storing the tetrahedra from left to right around the
    // boundary ring.  We must use a list and not a set, since we will
    // rely on the tetrahedra being stored in a particular order.
    std::list<NTetrahedron*> foundSoFar;
    foundSoFar.push_back(annulus.tet[0]);
    foundSoFar.push_back(middle);
    foundSoFar.push_back(annulus.tet[1]);

    // Also make a list of tetrahedron vertex roles for the two
    // tetrahedra in each segment that meet the boundary annuli.
    std::list<NPerm4> rolesSoFar;
    rolesSoFar.push_back(annulus.roles[0]);
    rolesSoFar.push_back(annulus.roles[1]);

    unsigned length = 1;
    bool twisted = false;

    NTetrahedron *nextLeft, *nextMiddle, *nextRight;
    NPerm4 nextLeftRoles, nextMiddleRoles, nextRightRoles;

    while (1) {
        // Run off the right hand side looking for the next tetrahedron.
        nextLeft = foundSoFar.back()->adjacentTetrahedron(
            rolesSoFar.back()[2]);
        nextLeftRoles = foundSoFar.back()->adjacentGluing(
            rolesSoFar.back()[2]) * rolesSoFar.back() * NPerm4(0, 1);

        if (nextLeft == annulus.tet[0]) {
            // The right _might_ have completed!
            if (nextLeftRoles == annulus.roles[0]) {
                // All good!  An untwisted strip.
            } else if (nextLeftRoles == annulus.roles[0] * NPerm4(0, 1)) {
                // A complete twisted strip.
                twisted = true;
            } else {
                // Nothing.
                return 0;
            }

            NSatReflectorStrip* ans = new NSatReflectorStrip(length, twisted);

            std::copy(foundSoFar.begin(), foundSoFar.end(),
                std::inserter(avoidTets, avoidTets.begin()));

            std::list<NTetrahedron*>::const_iterator tit = foundSoFar.begin();
            std::list<NPerm4>::const_iterator pit = rolesSoFar.begin();
            for (unsigned i = 0; i < length; i++) {
                ans->annulus_[i].tet[0] = *tit++;
                tit++; // Skip the middle tetrahedron from each block.
                ans->annulus_[i].tet[1] = *tit++;

                ans->annulus_[i].roles[0] = *pit++;
                ans->annulus_[i].roles[1] = *pit++;
            }

            return ans;
        }

        // Look for a new adjacent block.
        if (notUnique(nextLeft) ||
                isBad(nextLeft, avoidTets) || isBad(nextLeft, foundSoFar))
            return 0;

        nextMiddle = nextLeft->adjacentTetrahedron(nextLeftRoles[0]);
        nextMiddleRoles = nextLeft->adjacentGluing(
            nextLeftRoles[0]) * nextLeftRoles * NPerm4(3, 1, 0, 2);

        if (notUnique(nextMiddle, nextLeft) ||
                isBad(nextMiddle, avoidTets) || isBad(nextMiddle, foundSoFar))
            return 0;

        if (nextMiddle != nextLeft->adjacentTetrahedron(nextLeftRoles[1]))
            return 0;
        if (nextMiddleRoles != nextLeft->adjacentGluing(
                nextLeftRoles[1]) * nextLeftRoles * NPerm4(1, 3))
            return 0;

        nextRight = nextMiddle->adjacentTetrahedron(nextMiddleRoles[0]);
        nextRightRoles = nextMiddle->adjacentGluing(
            nextMiddleRoles[0]) * nextMiddleRoles * NPerm4(0, 3, 1, 2);

        if (notUnique(nextRight, nextLeft, nextMiddle) ||
                isBad(nextRight, avoidTets) || isBad(nextRight, foundSoFar))
            return 0;

        if (nextRight != nextMiddle->adjacentTetrahedron(nextMiddleRoles[1]))
            return 0;
        if (nextRightRoles != nextMiddle->adjacentGluing(
                nextMiddleRoles[1]) * nextMiddleRoles * NPerm4(0, 2))
            return 0;

        // Yup, we have a new block.
        foundSoFar.push_back(nextLeft);
        foundSoFar.push_back(nextMiddle);
        foundSoFar.push_back(nextRight);

        rolesSoFar.push_back(nextLeftRoles);
        rolesSoFar.push_back(nextRightRoles);

        length++;
    }

    // We should never get out of the loop this way.
    return 0;
}

NSatReflectorStrip* NSatReflectorStrip::insertBlock(NTriangulation& tri,
        unsigned length, bool twisted) {
    NSatReflectorStrip* ans = new NSatReflectorStrip(length, twisted);

    const NPerm4 id;
    NTetrahedron *upper, *lower, *middle;
    NTetrahedron *prevRight = 0, *firstLeft = 0;
    for (unsigned i = 0; i < length; i++) {
        // Create the three tetrahedra behind boundary annulus #i.
        upper = tri.newTetrahedron();
        lower = tri.newTetrahedron();
        middle = tri.newTetrahedron();

        upper->joinTo(0, middle, NPerm4(2, 1, 3, 0));
        lower->joinTo(0, middle, NPerm4(0, 3, 1, 2));
        upper->joinTo(1, middle, NPerm4(1, 3));
        lower->joinTo(1, middle, NPerm4(0, 2));

        if (i == 0)
            firstLeft = upper;
        else
            upper->joinTo(2, prevRight, NPerm4(0, 1));

        prevRight = lower;

        ans->annulus_[i].tet[0] = upper;
        ans->annulus_[i].tet[1] = lower;
        ans->annulus_[i].roles[0] = id;
        ans->annulus_[i].roles[1] = id;
    }

    if (twisted)
        firstLeft->joinTo(2, prevRight, id);
    else
        firstLeft->joinTo(2, prevRight, NPerm4(0, 1));

    return ans;
}

void NSatLayering::adjustSFS(NSFSpace& sfs, bool reflect) const {
    if (overHorizontal_)
        sfs.insertFibre(1, reflect ? -2 : 2);

    // Over the diagonal, there is no change at all.
}

NSatLayering* NSatLayering::isBlockLayering(const NSatAnnulus& annulus,
        TetList& avoidTets) {
    // Must be a common usable tetrahedron.
    if (annulus.tet[0] != annulus.tet[1])
        return 0;
    if (isBad(annulus.tet[0], avoidTets))
        return 0;

    // Is it a layering over the horizontal edge?
    if (annulus.roles[0][0] == annulus.roles[1][2] &&
            annulus.roles[0][2] == annulus.roles[1][0]) {
        avoidTets.insert(annulus.tet[0]);

        NSatLayering* ans = new NSatLayering(true);
        ans->annulus_[0] = annulus;
        ans->annulus_[1].tet[0] = ans->annulus_[1].tet[1] = annulus.tet[0];
        ans->annulus_[1].roles[0] = annulus.roles[1] * NPerm4(1, 0, 3, 2);
        ans->annulus_[1].roles[1] = annulus.roles[0] * NPerm4(1, 0, 3, 2);

        return ans;
    }

    // Is it a layering over the diagonal edge?
    if (annulus.roles[0][1] == annulus.roles[1][2] &&
            annulus.roles[0][2] == annulus.roles[1][1]) {
        avoidTets.insert(annulus.tet[0]);

        NSatLayering* ans = new NSatLayering(false);
        ans->annulus_[0] = annulus;
        ans->annulus_[1].tet[0] = ans->annulus_[1].tet[1] = annulus.tet[0];
        ans->annulus_[1].roles[0] = annulus.roles[1] * NPerm4(1, 0, 3, 2);
        ans->annulus_[1].roles[1] = annulus.roles[0] * NPerm4(1, 0, 3, 2);

        return ans;
    }

    // No layering.
    return 0;
}

} // namespace regina

