
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
 * It is assumed that the user is interested in finding minimal
 * triangulations of 3-manifolds.  With this in mind, the given topology
 * data file is read and processed as follows.
 *
 * Each triangulation will be modified by:
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
 * As a result:
 *
 * - some triangulations will be found to be non-minimal;
 * - some triangulations will be found to be homeomorphic to others in
 *   the file;
 * - some triangulations will be found to be homeomorphic to others not
 *   in the file.
 *
 * Any triangulation found to be non-minimal will be described as such
 * in the output.  The remaining triangulations (including any new
 * triangulations that were not in the data file) will be grouped into
 * equivalence classes according to which have been found to be homeomorphic
 * to which others.
 *
 * All of these results will be written to standard output.  Progress
 * reporting and other diagnostic messages will be written to standard
 * error.
 *
 * The results can also be optionally saved to an output data file.  In
 * this output file, each equivalence class of triangulations will be
 * stored beneath a separate container packet.  Triangulations that were
 * found to be non-minimal will not be saved at all.
 */

#include <file/nxmlfile.h>
#include <packet/ncontainer.h>
#include <triangulation/ntriangulation.h>

#include <cstdio>
#include <list>
#include <map>
#include <popt.h>
#include <set>
#include <sstream>

using namespace regina;

// Maximum number of moves of each type to make.
int argUp = 1;
int argAcross = 1;
int argDown = 1;

// The output file to write if required.
const char* outFile = 0;

// The input and output packet trees.
NPacket* tree = 0;
NPacket* newTree = 0;

// The original triangulation currently being processed.
NTriangulation* orig;

// Triangulations that have been found to be homeomorphic to orig.
// This list includes orig itself.
typedef std::set<NTriangulation*> TriSet;
TriSet equivs;

// Do we know that the original triangulation is non-minimal?
bool nonMin;

// All triangulations that were found to be non-minimal.
std::list<NTriangulation*> allNonMin;

// A mapping from triangulations to equivalence classes.
typedef std::map<NTriangulation*, int> ClassMap;
ClassMap eClass;
int nextClass = 0;

// Statistics.
unsigned long nTris = 0;
unsigned long nClasses = 0;
unsigned long nNonMin = 0;
unsigned long nNew = 0;

/**
 * We have a homeomorphic triangulation of the same size as the original.
 */
void sameSize(NTriangulation* t) {
    // Hunt for it in the packet tree.
    NTriangulation* found = 0;
    for (NPacket* p = tree; p; p = p->nextTreePacket())
        if (p->getPacketType() == NTriangulation::packetType)
            if (static_cast<NTriangulation*>(p)->isIsomorphicTo(*t).get()) {
                found = static_cast<NTriangulation*>(p);
                break;
            }

    if (found == orig)
        return;

    if (! found) {
        found = new NTriangulation(*t);
        found->setPacketLabel(orig->makeUniqueLabel(
            orig->getPacketLabel() + " - New"));
        orig->insertChildLast(found);
        nNew++;
    }

    // Yes, it's a set and so the insert should just fail if it's
    // already present.  But this is leading to a small memory leak on
    // my machine which, when multiplied by the _very_ large number of
    // times this is called, has rather unpleasant effects (such as
    // exhausting all swap space).  So let's keep the test here for now.
    if (! equivs.count(found))
        equivs.insert(found);
}

/**
 * Do the final greedy simplify and process.
 */
void processAlt(NTriangulation* t) {
    t->intelligentSimplify();

    if (t->getNumberOfTetrahedra() < orig->getNumberOfTetrahedra())
        nonMin = true;
    else if (t->getNumberOfTetrahedra() == orig->getNumberOfTetrahedra())
        sameSize(t);
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
    TriSet::iterator tit;
    ClassMap::iterator cit, cit2;
    int c, cOld;
    NTriangulation* t;

    for (NPacket* p = tree; p; p = p->nextTreePacket())
        if (p->getPacketType() == NTriangulation::packetType) {
            // A triangulation to process.
            t = static_cast<NTriangulation*>(p);
            fprintf(stderr, "Processing %s ...\n", t->getPacketLabel().c_str());
            nTris++;

            nonMin = false;
            orig = static_cast<NTriangulation*>(p);
            equivs.clear();
            equivs.insert(orig);

            tryMovesUp(orig, argUp);

            if (nonMin) {
                allNonMin.push_back(orig);
                nNonMin++;
                continue;
            }

            // In equivs we now have a list of all triangulations
            // equivalent to orig.

            // Is this an equivalence class we're already seen?
            for (tit = equivs.begin(); tit != equivs.end(); tit++) {
                cit = eClass.find(*tit);
                if (cit != eClass.end())
                    break;
            }
            if (tit != equivs.end()) {
                // We found an equivalence class.  Insert everything we
                // haven't seen yet, and merge the classes of everything
                // we have.
                c = cit->second;
                for (tit = equivs.begin(); tit != equivs.end(); tit++) {
                    cit = eClass.find(*tit);
                    if (cit == eClass.end())
                        eClass.insert(std::make_pair(*tit, c));
                    else if (cit->second != c) {
                        // Merge the two equivalence classes.
                        cOld = cit->second;
                        for (cit = eClass.begin(); cit != eClass.end(); cit++)
                            if (cit->second == cOld)
                                cit->second = c;
                        nClasses--;
                    }
                }
            } else {
                // No such equivalence class.  Insert everything.
                c = nextClass++;
                for (tit = equivs.begin(); tit != equivs.end(); tit++)
                    eClass.insert(std::make_pair(*tit, c));
                nClasses++;
            }
        }

    // Finished progress reporting.
    fprintf(stderr, "\n");

    // Write the summary of results.
    if (! allNonMin.empty()) {
        printf("NON-MINIMAL TRIANGULATIONS:\n\n");
        for (std::list<NTriangulation*>::const_iterator it = allNonMin.begin();
                it != allNonMin.end(); it++)
            printf("    %s\n", (*it)->getPacketLabel().c_str());
        printf("\n");
    }

    if (nClasses) {
        printf("EQUIVALENCE CLASSES:\n\n");

        if (outFile) {
            newTree = new NContainer();
            newTree->setPacketLabel("Equivalence Classes");
        }

        int classNum = 1;
        std::string className;
        NContainer* classCnt = 0;

        for (cit = eClass.begin(); cit != eClass.end(); cit++)
            if (cit->second >= 0) {
                // The first triangulation of a new equivalence class.
                c = cit->second;

                std::ostringstream s;
                s << "Class " << classNum << " : " <<
                    cit->first->getHomologyH1().toString();
                className = s.str();
                classNum++;

                printf("%s\n\n", className.c_str());
                if (outFile) {
                    classCnt = new NContainer();
                    classCnt->setPacketLabel(className);
                    newTree->insertChildLast(classCnt);
                }

                // Find the triangulations in this class, and erase the
                // class as we go.
                for (cit2 = cit; cit2 != eClass.end(); cit2++)
                    if (cit2->second == c) {
                        printf("    %s\n",
                            cit2->first->getPacketLabel().c_str());
                        if (outFile) {
                            t = new NTriangulation(*(cit2->first));
                            t->setPacketLabel(cit2->first->getPacketLabel());
                            classCnt->insertChildLast(t);
                        }

                        cit2->second = -1;
                    }

                printf("\n");
            }
    }

    printf("Final statistics:\n");
    printf("    Triangulations read:        %ld\n", nTris);
    printf("    Equivalence classes:        %ld\n", nClasses);
    printf("    New triangulations:         %ld\n", nNew);
    printf("    Non-minimal triangulations: %ld\n", nNonMin);
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
            "Output equivalence classes (plus new triangulations) to the given file.", "<output-file>" },
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
    if (outFile && newTree) {
        fprintf(stderr, "\nSaving results to %s...\n", outFile);
        writeXMLFile(outFile, newTree);
    } else
        fprintf(stderr, "\nNot saving results.\n");

    // Clean up.
    delete tree;
    if (newTree)
        delete newTree;

    return 0;
}

