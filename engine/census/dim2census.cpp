
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

#include <algorithm>
#include <sstream>
#include "census/dim2census.h"
#include "census/dim2gluingpermsearcher.h"
#include "dim2/dim2triangulation.h"
#include "utilities/memutils.h"

namespace regina {

unsigned long Dim2Census::formCensus(NPacket* parent, unsigned nTriangles,
        NBoolSet orientability, NBoolSet boundary,
        int nBdryEdges, AcceptTriangulation sieve,
        void* sieveArgs) {
    // If obviously nothing is going to happen but we won't realise
    // it until we've actually generated the facet pairings, change
    // nTriangles to 0 so we'll realise it immediately once the new
    // thread starts.
    if (orientability == NBoolSet::sNone)
        nTriangles = 0;

    // Start the census!
    Dim2Census* census = new Dim2Census(parent, orientability,
        sieve, sieveArgs);

    Dim2EdgePairing::findAllPairings(nTriangles, boundary, nBdryEdges,
        Dim2Census::foundEdgePairing, census);
    unsigned long ans = census->whichSoln_ - 1;
    delete census;
    return ans;
}

unsigned long Dim2Census::formPartialCensus(const Dim2EdgePairing* pairing,
        NPacket* parent, NBoolSet orientability,
        AcceptTriangulation sieve, void* sieveArgs) {
    // Is it obvious that nothing will happen?
    if (orientability == NBoolSet::sNone)
        return 0;

    // Make a list of automorphisms.
    Dim2EdgePairing::IsoList autos;
    pairing->findAutomorphisms(autos);

    // Select the individual gluing permutations.
    Dim2Census census(parent, orientability, sieve, sieveArgs);
    Dim2GluingPermSearcher::findAllPerms(pairing, &autos,
        ! census.orientability_.hasFalse(),
        Dim2Census::foundGluingPerms, &census);

    // Clean up.
    std::for_each(autos.begin(), autos.end(), FuncDelete<Dim2Isomorphism>());
    return census.whichSoln_ - 1;
}

Dim2Census::Dim2Census(NPacket* parent,
        const NBoolSet& orientability, AcceptTriangulation sieve,
        void* sieveArgs) :
        parent_(parent), orientability_(orientability),
        sieve_(sieve), sieveArgs_(sieveArgs),
        whichSoln_(1) {
}

void Dim2Census::foundEdgePairing(const Dim2EdgePairing* pairing,
        const Dim2EdgePairing::IsoList* autos, void* census) {
    Dim2Census* realCensus = static_cast<Dim2Census*>(census);
    if (pairing) {
        // We've found another edge pairing.
        // Select the individual gluing permutations.
        Dim2GluingPermSearcher::findAllPerms(pairing, autos,
            ! realCensus->orientability_.hasFalse(),
            Dim2Census::foundGluingPerms, census);
    } else {
        // Census generation has finished.
    }
}

void Dim2Census::foundGluingPerms(const Dim2GluingPermSearcher* perms,
        void* census) {
    if (perms) {
        // We've found another permutation set.
        // Triangulate and see what we've got.
        Dim2Triangulation* tri = perms->triangulate();
        Dim2Census* realCensus = static_cast<Dim2Census*>(census);

        bool ok = true;
        if ((! realCensus->orientability_.hasTrue()) && tri->isOrientable())
            ok = false;
        else if (realCensus->sieve_ &&
                ! realCensus->sieve_(tri, realCensus->sieveArgs_))
            ok = false;

        if (ok) {
            // Put it in the census!
            // Make sure it has a charming label.
            // Don't insist on unique labels, since this requirement
            // will soon be dropped and it multiplies the running time
            // by a factor of #triangulations.
            std::ostringstream out;
            out << "Item " << realCensus->whichSoln_;
            tri->setLabel(out.str());
            realCensus->parent_->insertChildLast(tri);
            realCensus->whichSoln_++;
        } else {
            // Bad triangulation.
            delete tri;
        }
    }
}

} // namespace regina

