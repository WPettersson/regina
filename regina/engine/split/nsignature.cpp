
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
#include <cctype>
#include "split/nsignature.h"
#include "triangulation/ntriangulation.h"
#include "utilities/memutils.h"

namespace regina {

namespace {
    NPerm4 exitFace(bool firstOccurrence, bool lowerCase) {
        if (firstOccurrence) {
            if (lowerCase)
                return NPerm4(2,3,1,0);
            else
                return NPerm4(2,3,0,1);
        } else {
            if (lowerCase)
                return NPerm4(0,1,3,2);
            else
                return NPerm4(0,1,2,3);
        }
    }
}

NSignature::NSignature(const NSignature& sig) : ShareableObject(),
        order(sig.order), label(new unsigned[2 * sig.order]),
        labelInv(new bool[2 * sig.order]), nCycles(sig.nCycles),
        cycleStart(new unsigned[sig.nCycles + 1]),
        nCycleGroups(sig.nCycleGroups),
        cycleGroupStart(new unsigned[sig.nCycleGroups + 1]) {
    std::copy(sig.label, sig.label + 2 * sig.order, label);
    std::copy(sig.labelInv, sig.labelInv + 2 * sig.order, labelInv);
    std::copy(sig.cycleStart, sig.cycleStart + sig.nCycles + 1, cycleStart);
    std::copy(sig.cycleGroupStart, sig.cycleGroupStart + sig.nCycleGroups + 1,
        cycleGroupStart);
}

NSignature* NSignature::parse(const std::string& str) {
    // See if the string looks correctly formed.
    // Note that we're not yet counting the individual frequency of each
    // letter, just the overall number of letters.
    // Cycles are assumed to be separated by any non-space
    // non-alphabetic characters.
    unsigned nAlpha = 0;
    int largestLetter = -1;

    size_t len = str.length();
    size_t pos;
    for (pos = 0; pos < len; pos++)
        // Avoid isalpha(), etc. and be explicit, in case the signature
        // string contains international characters.
        if (str[pos] >= 'A' && str[pos] <= 'Z') {
            nAlpha++;
            if (largestLetter < str[pos] - 'A')
                largestLetter = str[pos] - 'A';
        } else if (str[pos] >= 'a' && str[pos] <= 'z') {
            nAlpha++;
            if (largestLetter < str[pos] - 'a')
                largestLetter = str[pos] - 'a';
        }

    if (static_cast<int>(nAlpha) != 2 * (largestLetter + 1))
        return 0;
    if (nAlpha == 0)
        return 0;

    // Looks fine so far.
    // Build the signature and cycle structure (but not cycle groups yet).
    unsigned order = largestLetter + 1;
    unsigned* label = new unsigned[nAlpha];
    bool* labelInv = new bool[nAlpha];
    unsigned nCycles = 0;
    unsigned* cycleStart = new unsigned[nAlpha + 1];
    cycleStart[0] = 0;

    unsigned* freq = new unsigned[order];
    std::fill(freq, freq + order, 0);

    unsigned whichPos = 0;
        /* Position in the signature, as opposed to position in the string. */
    unsigned letterIndex;
    for (pos = 0; pos < len; pos++) {
        if (isspace(str[pos]))
            continue;
        if (! ((str[pos] >= 'A' && str[pos] <= 'Z') ||
                (str[pos] >= 'a' && str[pos] <= 'z'))) {
            if (cycleStart[nCycles] < whichPos) {
                // We've just ended a cycle.
                nCycles++;
                cycleStart[nCycles] = whichPos;
            }
        } else {
            if (str[pos] >= 'A' && str[pos] <= 'Z')
                letterIndex = str[pos] - 'A';
            else
                letterIndex = str[pos] - 'a';

            freq[letterIndex]++;
            if (freq[letterIndex] > 2) {
                // We've seen this letter a third time!
                delete[] label;
                delete[] labelInv;
                delete[] cycleStart;
                delete[] freq;
                return 0;
            }
            label[whichPos] = letterIndex;
            labelInv[whichPos] = (str[pos] >= 'A' && str[pos] <= 'Z');
            whichPos++;
        }
    }

    delete[] freq;

    if (cycleStart[nCycles] < whichPos) {
        nCycles++;
        cycleStart[nCycles] = whichPos;
    }

    // We now have a valid signature!
    NSignature* sig = new NSignature();
    sig->order = order;
    sig->label = label;
    sig->labelInv = labelInv;
    sig->nCycles = nCycles;
    sig->cycleStart = cycleStart;

    // Fill in the rest of the data members.
    sig->nCycleGroups = 0;
    sig->cycleGroupStart = new unsigned[nCycles];
    for (pos = 0; pos < nCycles; pos++)
        if (pos == 0 || sig->cycleStart[pos + 1] - sig->cycleStart[pos] !=
                sig->cycleStart[pos] - sig->cycleStart[pos - 1]) {
            // New cycle group.
            sig->cycleGroupStart[sig->nCycleGroups] =
                static_cast<unsigned>(pos);
            sig->nCycleGroups++;
        }

    return sig;
}

NTriangulation* NSignature::triangulate() const {
    unsigned sigLen = 2 * order;
    NTriangulation* tri = new NTriangulation();

    // Create a new set of tetrahedra.
    // Tetrahedron vertices will be:
    //   bottom left -> top right: 0 -> 1
    //   bottom right -> top left: 2 -> 3
    NTetrahedron** tet = new NTetrahedron*[order];
    unsigned pos;
    for (pos = 0; pos < order; pos++)
        tet[pos] = tri->newTetrahedron();

    // Store the first occurrence of each symbol.
    unsigned* first = new unsigned[order];
    std::fill(first, first + order, sigLen);

    for (pos = 0; pos < sigLen; pos++)
        if (first[label[pos]] == sigLen)
            first[label[pos]] = pos;

    // Make the face gluings.
    unsigned currCycle = 0;
    unsigned adjPos;
    NPerm4 myFacePerm, yourFacePerm;
    for (pos = 0; pos < sigLen; pos++) {
        if (cycleStart[currCycle + 1] == pos + 1) {
            adjPos = cycleStart[currCycle];
            currCycle++;
        } else
            adjPos = pos + 1;

        myFacePerm = exitFace(first[label[pos]] == pos, ! labelInv[pos]);
        yourFacePerm = exitFace(first[label[adjPos]] == adjPos,
            labelInv[adjPos]);
        tet[label[pos]]->joinTo(myFacePerm[3], tet[label[adjPos]],
            yourFacePerm * myFacePerm.inverse());
    }

    // Clean up.
    delete[] first;
    delete[] tet;
    return tri;
}

int NSignature::cycleCmp(const NSignature& sig1, unsigned cycle1,
        unsigned start1, int dir1, unsigned* relabel1, const NSignature& sig2,
        unsigned cycle2, unsigned start2, int dir2, unsigned* relabel2) {
    unsigned len = sig1.cycleStart[cycle1 + 1] - sig1.cycleStart[cycle1];
    unsigned* arr1 = sig1.label + sig1.cycleStart[cycle1];
    unsigned* arr2 = sig2.label + sig2.cycleStart[cycle2];
    unsigned pos1 = start1;
    unsigned pos2 = start2;
    for (unsigned i = 0; i < len; i++) {
        if ((relabel1 ? relabel1[arr1[pos1]] : arr1[pos1]) <
                (relabel2 ? relabel2[arr2[pos2]] : arr2[pos2]))
            return -1;
        if ((relabel1 ? relabel1[arr1[pos1]] : arr1[pos1]) >
                (relabel2 ? relabel2[arr2[pos2]] : arr2[pos2]))
            return 1;

        if (dir1 > 0) {
            pos1++;
            if (pos1 == len)
                pos1 = 0;
        } else {
            if (pos1 == 0)
                pos1 = len - 1;
            else
                pos1--;
        }
        if (dir2 > 0) {
            pos2++;
            if (pos2 == len)
                pos2 = 0;
        } else {
            if (pos2 == 0)
                pos2 = len - 1;
            else
                pos2--;
        }
    }
    return 0;
}

void NSignature::writeCycles(std::ostream& out, const std::string& cycleOpen,
        const std::string& cycleClose, const std::string& cycleJoin) const {
    out << cycleOpen;

    unsigned cycle = 0;
    for (unsigned pos = 0; pos < 2 * order; pos++) {
        if (cycleStart[cycle] == pos) {
            if (cycle > 0)
                out << cycleClose << cycleJoin << cycleOpen;
            cycle++;
        }
        out << char((labelInv[pos] ? 'A' : 'a') + label[pos]);
    }

    out << cycleClose;
}

} // namespace regina
