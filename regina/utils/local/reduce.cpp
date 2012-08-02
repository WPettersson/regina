
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Attempt to find triangulations related by few elementary moves        *
 *                                                                        *
 *  Copyright (c) 2005-2011, Ben Burton                                   *
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

/**
 * The given topology data file is read, and an attempt is made to
 * reduce the size of each triangulation found within it (specifically,
 * to reduce the number of tetrahedra used).
 *
 * In trying to reduce triangulations, the following steps are performed:
 *
 * (i) making up to a given number of 2-3 moves to expand the triangulation;
 * (ii) making up to a given number of 4-4 moves to modify the triangulation;
 * (iii) making up to a given number of reduction moves to simplify the
 * triangulation;
 * (iv) doing a greedy simplification to further reduce the size of the
 * triangulation if possible.
 *
 * Steps (i), (ii) and (iii) will be tried in all possible ways.  This
 * means that the program will be infeasibly slow unless the maximum numbers
 * of 2-3, 4-4 and reduction moves are very small.
 *
 * If a smaller triangulation is found as a result, it will be stored
 * beneath the original in the packet tree (and no further attempts will
 * be made to reduce it).  The original triangulations will remain
 * unchanged.
 *
 * Progress reporting and final statistics will be written to standard output.
 * If the option -o is passed, the resulting packet tree (including original
 * and reduced triangulations) will be saved to the given output data file.
 */

#include <file/nxmlfile.h>
#include <packet/ncontainer.h>
#include <triangulation/ntriangulation.h>

#include <cstdio>
#include <popt.h>

using namespace regina;

// Maximum number of moves of each type to make.
int argUp = 1;
int argAcross = 1;
int argDown = 1;

// The output file to write if required.
const char* outFile = 0;

// The input packet tree.
NPacket* tree = 0;

// The original triangulation currently being processed.
NTriangulation* orig;

// Do we know that the original triangulation is non-minimal?
bool nonMin;

// Statistics.
unsigned long nTris = 0;
unsigned long nReduced = 0;

/**
 * Do the final greedy simplify and process.
 */
void processAlt(NTriangulation* t) {
    t->intelligentSimplify();

    if (t->getNumberOfTetrahedra() >= orig->getNumberOfTetrahedra())
        return;

    // Got one.
    nonMin = true;

    NTriangulation* final = new NTriangulation(*t);
    final->setPacketLabel(tree->makeUniqueLabel(orig->getPacketLabel() +
        " - reduced"));

    // Insert as the first child so it doesn't get processed.
    orig->insertChildFirst(final);
}

/**
 * Perform reduction moves.  The given triangulation may be changed.
 */
void tryMovesDown(NTriangulation* t, int maxLevels) {
    if (maxLevels == 0) {
        processAlt(t);
        return;
    }

    NTriangulation* alt;
    unsigned i, j;
    bool found = false;

    for (i = 0; i < t->getNumberOfEdges(); i++)
        if (t->twoZeroMove(t->getEdge(i), true, false)) {
            alt = new NTriangulation(*t);
            alt->twoZeroMove(alt->getEdge(i));
            tryMovesDown(alt, maxLevels - 1);
            found = true;
            delete alt;

            if (nonMin)
                return;
        }

    for (i = 0; i < t->getNumberOfEdges(); i++)
        for (j = 0; j < 2; j++)
            if (t->twoOneMove(t->getEdge(i), j, true, false)) {
                alt = new NTriangulation(*t);
                alt->twoOneMove(alt->getEdge(i), j);
                tryMovesDown(alt, maxLevels - 1);
                found = true;
                delete alt;

                if (nonMin)
                    return;
            }

    // Only try 3-2 moves if nothing better has worked so far.
    if (! found)
        for (i = 0; i < t->getNumberOfEdges(); i++)
            if (t->threeTwoMove(t->getEdge(i), true, false)) {
                alt = new NTriangulation(*t);
                alt->threeTwoMove(alt->getEdge(i));
                tryMovesDown(alt, maxLevels - 1);
                found = true;
                delete alt;

                if (nonMin)
                    return;
            }

    // Only try 4-4 moves if nothing else has worked.
    if (! found)
        for (i = 0; i < t->getNumberOfEdges(); i++)
            for (j = 0; j < 2; j++)
                if (t->fourFourMove(t->getEdge(i), j, true, false)) {
                    alt = new NTriangulation(*t);
                    alt->fourFourMove(alt->getEdge(i), j);
                    tryMovesDown(alt, maxLevels - 1);
                    found = true;
                    delete alt;

                    if (nonMin)
                        return;
                }

    // No more moves to make.  Run straight to intelligentSimplify().
    if (! found)
        processAlt(t);
}

/**
 * Perform 4-4 moves.  The given triangulation may be changed.
 * Moves that revert to prev, prev2 or prev3 will not be considered.
 */
void tryMovesAcross(NTriangulation* t, int maxLevels,
        NTriangulation* prev = 0, NTriangulation* prev2 = 0,
        NTriangulation* prev3 = 0) {
    unsigned i, j;
    NTriangulation* alt;

    if (maxLevels > 0)
        for (i = 0; i < t->getNumberOfEdges(); i++)
            for (j = 0; j < 2; j++)
                if (t->fourFourMove(t->getEdge(i), j, true, false)) {
                    alt = new NTriangulation(*t);
                    alt->fourFourMove(alt->getEdge(i), j);
                    if (prev && alt->isIsomorphicTo(*prev).get()) {
                        // Ignore, reversion.
                    } else if (prev2 && alt->isIsomorphicTo(*prev2).get()) {
                        // Ignore, reversion.
                    } else if (prev3 && alt->isIsomorphicTo(*prev3).get()) {
                        // Ignore, reversion.
                    } else
                        tryMovesAcross(alt, maxLevels - 1, t, prev, prev2);
                    delete alt;

                    if (nonMin)
                        return;
                }

    // Try just going for the simplify.
    tryMovesDown(t, argDown);
}

/**
 * Perform 2-3 moves.  The given triangulation will not be changed.
 */
void tryMovesUp(NTriangulation* t, int levelsRemaining) {
    NTriangulation* alt;

    if (levelsRemaining == 0) {
        // We're not allowed to change the original, so clone it.
        alt = new NTriangulation(*t);
        tryMovesAcross(alt, argAcross);
        delete alt;
    } else {
        for (unsigned i = 0; i < t->getNumberOfFaces(); i++) {
            alt = new NTriangulation(*t);
            if (alt->twoThreeMove(alt->getFace(i))) {
                if (levelsRemaining > 1)
                    tryMovesUp(alt, levelsRemaining - 1);
                else
                    tryMovesAcross(alt, argAcross);
            }
            delete alt;

            if (nonMin)
                return;
        }
    }
}

void processTree() {
    NPacket* p = tree;
    NPacket* next;

    while (p) {
        next = p->nextTreePacket();

        if (p->getPacketType() == NTriangulation::packetType) {
            // A triangulation to process.
            orig = static_cast<NTriangulation*>(p);
            printf("Processing %s ...\n", orig->getPacketLabel().c_str());
            nTris++;

            nonMin = false;
            tryMovesUp(orig, argUp);
            if (nonMin)
                nReduced++;
        }

        p = next;
    }

    // Finished progress reporting.
    printf("\n");

    printf("Final statistics:\n");
    printf("    Triangulations read:    %ld\n", nTris);
    printf("    Triangulations reduced: %ld\n", nReduced);
}

int main(int argc, const char* argv[]) {
    // Set up the command-line arguments.
    poptOption opts[] = {
        { "up", 'u', POPT_ARG_INT, &argUp, 0,
            "Number of initial 2-3 moves (default is 1).", "<moves>" },
        { "across", 'a', POPT_ARG_INT, &argAcross, 0,
            "Number of subsequent 4-4 moves (default is 1).", "<moves>" },
        { "down", 'd', POPT_ARG_INT, &argDown, 0,
            "Number of final reduction moves before the greedy simplification (default is 1).", "<moves>" },
        { "output", 'o', POPT_ARG_STRING, &outFile, 0,
            "Output original and reduced triangulations to the given file.",
            "<output-file>" },
        POPT_AUTOHELP
        { 0, 0, 0, 0, 0, 0, 0 }
    };

    poptContext optCon = poptGetContext(0, argc, argv, opts, 0);
    poptSetOtherOptionHelp(optCon, "<file.rga>");

    // Parse the command-line arguments.
    int rc = poptGetNextOpt(optCon);
    if (rc != -1) {
        fprintf(stderr, "%s: %s\n\n",
            poptBadOption(optCon, POPT_BADOPTION_NOALIAS), poptStrerror(rc));
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return 1;
    }

    std::string filename;
    const char** otherOpts = poptGetArgs(optCon);
    if (otherOpts && otherOpts[0]) {
        filename = otherOpts[0];
        if (otherOpts[1]) {
            fprintf(stderr, "Only one filename may be supplied.\n\n");
            poptPrintHelp(optCon, stderr, 0);
            poptFreeContext(optCon);
            return 1;
        }
    } else {
        fprintf(stderr, "No filename was supplied.\n\n");
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return 1;
    }

    // Run a sanity check on the command-line arguments.
    bool broken = false;
    if (argUp < 0) {
        fprintf(stderr, "The number of 2-3 moves up may not be negative.\n");
        broken = true;
    }
    if (argAcross < 0) {
        fprintf(stderr, "The number of 4-4 moves across may not be "
            "negative.\n");
        broken = true;
    }
    if (argDown < 0) {
        fprintf(stderr, "The number of reduction moves down may not be "
            "negative.\n");
        broken = true;
    }

    if (broken) {
        fprintf(stderr, "\n");
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return 1;
    }

    // Done parsing the command line.
    poptFreeContext(optCon);

    // Read the data file.
    if (! (tree = readXMLFile(filename.c_str()))) {
        fprintf(stderr, "ERROR: Could not read data from %s.\n",
            filename.c_str());
        return 1;
    }

    // Off we go.
    processTree();

    // Are we saving results?
    if (outFile) {
        printf("\nSaving results to %s...\n", outFile);
        writeXMLFile(outFile, tree);
    } else
        printf("\nNot saving results.\n");

    // Clean up.
    delete tree;
    return 0;
}

