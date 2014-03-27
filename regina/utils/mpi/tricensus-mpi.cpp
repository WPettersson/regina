
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Manage a distributed census of triangulations using MPI               *
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

#define SUPPORT_DIM4 0

#include "mpi.h"

#include <cctype>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <popt.h>
#include "census/dim2census.h"
#include "census/dim2gluingpermsearcher.h"
#if SUPPORT_DIM4
#include "census/dim4census.h"
#include "census/dim4gluingpermsearcher.h"
#endif
#include "census/ncensus.h"
#include "census/ngluingpermsearcher.h"
#include "dim2/dim2triangulation.h"
#if SUPPORT_DIM4
#include "dim4/dim4triangulation.h"
#endif
#include "file/nxmlfile.h"
#include "packet/ncontainer.h"
#include "packet/ntext.h"
#include "triangulation/ntriangulation.h"
#include "../cycledecomp.h"
// Write messages tailored to the working dimension.
#define WORD_face (dim4 ? "facet" : dim2 ? "edge" : "face")
#define WORD_Face (dim4 ? "Facet" : dim2 ? "Edge" : "Face")

// MPI constants.
#define TAG_REQUEST_TASK 10
#define TAG_REQUEST_PAIRING 11
#define TAG_REQUEST_SUBSEARCH 12
#define TAG_RESULT 13

// Time constants:
#define MIN_SEC 60
#define HOUR_SEC (60 * MIN_SEC)
#define DAY_SEC (24 * HOUR_SEC)

// Forward declarations required for the templates below:
template <class CensusType>
void slaveFoundGluingPerms(const typename CensusType::GluingPermSearcher*,
    void*);

template <class CensusType>
void ctrlFarmPartialSearch(const typename CensusType::GluingPermSearcher*,
    void*);

// Differences between censuses of 2, 3 and 4-manifolds:
struct Dim2Params {
    typedef regina::Dim2EdgePairing Pairing;
    typedef regina::Dim2GluingPermSearcher GluingPermSearcher;
    typedef regina::Dim2Triangulation Triangulation;

    inline static GluingPermSearcher* bestSearcher(Pairing* p,
            bool orientableOnly, bool /* finiteOnly */, int /* whichPurge */) {
        return GluingPermSearcher::bestSearcher(p, 0 /* autos */,
            orientableOnly, ctrlFarmPartialSearch<Dim2Params>);
    }

    inline static void findAllPerms(Pairing* p, bool orientableOnly,
            bool /* finiteOnly */, int /* whichPurge */,
            regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, 0,
            orientableOnly, slaveFoundGluingPerms<Dim2Params>, dest);
    }

    inline static bool mightBeMinimal(Triangulation* tri) {
        return tri->isMinimal();
    }

    inline static const Pairing* pairingFor(const GluingPermSearcher* s) {
        return s->getFacetPairing();
    }
};

struct Dim3Params {
    typedef regina::NFacePairing Pairing;
    typedef regina::NGluingPermSearcher GluingPermSearcher;
    typedef regina::NTriangulation Triangulation;

    inline static GluingPermSearcher* bestSearcher(Pairing* p,
            bool orientableOnly, bool finiteOnly, int whichPurge) {
        return GluingPermSearcher::bestSearcher(p, 0 /* autos */,
            orientableOnly, finiteOnly, whichPurge,
            ctrlFarmPartialSearch<Dim3Params>);
    }

    inline static void findAllPerms(Pairing* p, bool orientableOnly,
            bool finiteOnly, int whichPurge, regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, 0,
            orientableOnly, finiteOnly, whichPurge,
            slaveFoundGluingPerms<Dim3Params>, dest);
    }

    inline static bool mightBeMinimal(Triangulation* tri) {
        return regina::NCensus::mightBeMinimal(tri, 0);
    }

    inline static const Pairing* pairingFor(const GluingPermSearcher* s) {
        return s->getFacePairing();
    }
};

struct Dim3DecompParams {
    typedef regina::NFacePairing Pairing;
    typedef CycleDecompSearcher GluingPermSearcher;
    typedef regina::NTriangulation Triangulation;

    inline static GluingPermSearcher* bestSearcher(Pairing* p,
            bool orientableOnly, bool finiteOnly, int whichPurge) {
        return  new CycleDecompSearcher(p, 0,
                orientableOnly, ctrlFarmPartialSearch<Dim3DecompParams>);
    }

    inline static void findAllPerms(Pairing* p, bool orientableOnly,
            bool finiteOnly, int whichPurge, regina::NPacket* dest) {
            CycleDecompSearcher *searcher = new CycleDecompSearcher(p, 0,
                orientableOnly, slaveFoundGluingPerms<Dim3DecompParams>, dest);
            searcher->runSearch();
            delete searcher;
    }

    inline static bool mightBeMinimal(Triangulation* tri) {
        return regina::NCensus::mightBeMinimal(tri, 0);
    }

    inline static const Pairing* pairingFor(const GluingPermSearcher* s) {
        return s->getFacetPairing();
    }
};

#if SUPPORT_DIM4
struct Dim4Params {
    typedef regina::Dim4FacetPairing Pairing;
    typedef regina::Dim4GluingPermSearcher GluingPermSearcher;
    typedef regina::Dim4Triangulation Triangulation;

    inline static GluingPermSearcher* bestSearcher(Pairing* p,
            bool orientableOnly, bool finiteOnly, int /* whichPurge */) {
        return GluingPermSearcher::bestSearcher(p, 0 /* autos */,
            orientableOnly, finiteOnly,
            ctrlFarmPartialSearch<Dim4Params>);
    }

    inline static void findAllPerms(Pairing* p, bool orientableOnly,
            bool finiteOnly, int /* whichPurge */, regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, 0,
            orientableOnly, finiteOnly,
            slaveFoundGluingPerms<Dim4Params>, dest);
    }

    inline static bool mightBeMinimal(Triangulation*) {
        return true;
    }

    inline static const Pairing* pairingFor(const GluingPermSearcher* s) {
        return s->getFacetPairing();
    }
};
#endif

// Census parameters.
regina::NBoolSet
    finiteness(true, true),
    orientability(true, true);
int minimal = 0;
int minimalPrime = 0;
int minimalPrimeP2 = 0;
int whichPurge = 0;
int dim2 = 0;
int dim4 = 0;
long depth = 0;
int dryRun = 0;
int sigs = 0;
int cycleDecomp = 0;

// Filenames read from the command line.
std::string outputStub;
std::string pairsFile;

// Task identification globals.
// For processing a face pairing: ( pairing ID, -1, text data length )
// For processing a partial search: ( pairing ID, search ID, text data length )
// For end of processing: ( -1, -1, -1 )
long taskID[3] = { 0, 0, 0 };

// Controller-specific globals.
std::ofstream logger;
bool controllerError = false;

struct Task {
    long pairing, subtask;
    time_t start;
};
Task* slaveTask;
int nSlaves, nRunningSlaves;
long totTri;

// Slave-specific globals.
long nSolns;
std::string sigFile;
std::ofstream sigStream;
bool sigStreamErr;

/**
 * Generic helper routine.
 *
 * Parse command-line arguments for census and filename options.
 */
int parseCmdLine(int argc, const char* argv[], bool isController) {
    // Set up the command-line arguments.
    int argOr = 0;
    int argNor = 0;
    int argFinite = 0;
    int argIdeal = 0;
    poptOption opts[] = {
        { "orientable", 'o', POPT_ARG_NONE, &argOr, 0,
            "Must be orientable.", 0 },
        { "nonorientable", 'n', POPT_ARG_NONE, &argNor, 0,
            "Must be non-orientable.", 0 },
        { "finite", 'f', POPT_ARG_NONE, &argFinite, 0,
            "Must be finite (no ideal vertices).", 0 },
        { "ideal", 'd', POPT_ARG_NONE, &argIdeal, 0,
            "Must have at least one ideal vertex.", 0 },
        { "minimal", 'm', POPT_ARG_NONE, &minimal, 0,
            "Ignore obviously non-minimal triangulations.", 0 },
        { "minprime", 'M', POPT_ARG_NONE, &minimalPrime, 0,
            "Ignore obviously non-minimal, non-prime and/or disc-reducible "
            "triangulations.", 0 },
        { "minprimep2", 'N', POPT_ARG_NONE, &minimalPrimeP2, 0,
            "Ignore obviously non-minimal, non-prime, disc-reducible and/or "
            "P2-reducible triangulations.", 0 },
        { "cycledecomp" , 'C', POPT_ARG_NONE, &cycleDecomp, 0,
            "Find permutations by cycle decompositions of the face pairing graph.", 0 },
        { "dim2", '2', POPT_ARG_NONE, &dim2, 0,
            "Run a census of 2-manifold triangulations, "
            "not 3-manifold triangulations.", 0 },
#if SUPPORT_DIM4
        { "dim4", '4', POPT_ARG_NONE, &dim4, 0,
            "Run a census of 4-manifold triangulations, "
            "not 3-manifold triangulations.", 0 },
#endif
        { "sigs", 's', POPT_ARG_NONE, &sigs, 0,
            "Write isomorphism signatures only, not full Regina data files.",
            0 },
        { "depth", 'D', POPT_ARG_LONG, &depth, 0,
            "Split each face pairing into subsearches at the given depth.",
            "<depth>" },
        { "dryrun", 'x', POPT_ARG_NONE, &dryRun, 0,
            "Have slaves ignore each task and simply report zero "
            "triangulations instead.  This allows a quick overview of the "
            "search space.", 0 },
        POPT_AUTOHELP
        { 0, 0, 0, 0, 0, 0, 0 }
    };

    poptContext optCon = poptGetContext(0, argc, argv, opts, 0);
    poptSetOtherOptionHelp(optCon, "<pairs-file> <output-stub>");

    // Parse the command-line arguments.
    int rc = poptGetNextOpt(optCon);
    if (rc != -1) {
        if (isController) {
            std::cerr << poptBadOption(optCon, POPT_BADOPTION_NOALIAS)
                << ": " << poptStrerror(rc) << "\n\n";
            poptPrintHelp(optCon, stderr, 0);
        }
        poptFreeContext(optCon);
        return 1;
    }

    const char** otherOpts = poptGetArgs(optCon);
    if (otherOpts && otherOpts[0] && otherOpts[1]) {
        pairsFile = otherOpts[0];
        outputStub = otherOpts[1];
        if (otherOpts[2]) {
            if (isController) {
                std::cerr << "Too many arguments.\n\n";
                poptPrintHelp(optCon, stderr, 0);
            }
            poptFreeContext(optCon);
            return 1;
        }
    } else {
        if (isController) {
            std::cerr << "Not enough arguments.\n\n";
            poptPrintHelp(optCon, stderr, 0);
        }
        poptFreeContext(optCon);
        return 1;
    }

    // Run a sanity check on the command-line arguments.
    bool broken = false;
    if (pairsFile.empty()) {
        if (isController)
            std::cerr << "A pairs file must be specified.\n";
        broken = true;
    } else if (outputStub.empty()) {
        if (isController)
            std::cerr << "An output filename stub must be specified.\n";
        broken = true;
    } else if (argOr && argNor) {
        if (isController)
            std::cerr << "Options -o/--orientable and -n/--nonorientable "
                << "cannot be used together.\n";
        broken = true;
    } else if (dim2 && dim4) {
        if (isController)
            std::cerr << "Options -2/--dim2 and -4/--dim4 "
                "cannot be used together.\n";
        broken = true;
    } else if (dim2 && (argFinite || argIdeal)) {
        if (isController)
            std::cerr << "Finiteness options cannot be used with -2/--dim2.\n";
        broken = true;
    } else if (dim2 && (minimalPrime || minimalPrimeP2)) {
        if (isController)
            std::cerr << "Primeness options cannot be used with -2/--dim2 "
                "(the weaker -m/--minimal can).\n";
        broken = true;
    } else if (dim4 && (minimal || minimalPrime || minimalPrimeP2)) {
        if (isController)
            std::cerr << "Minimality options cannot be used with -4/--dim4.\n";
        broken = true;
    } else if (cycleDecomp && (dim2 || dim4)) {
        std::cerr << "Cycle decompositions are only available in 3 "
            << "dimensions.\n";
        broken = true;
    } else if (depth < 0) {
        if (isController)
            std::cerr << "The branching depth must be non-negative.\n";
        broken = true;
    }

    if (broken) {
        if (isController) {
            std::cerr << '\n';
            poptPrintHelp(optCon, stderr, 0);
        }
        poptFreeContext(optCon);
        return 1;
    }

    // Done parsing the command line.
    poptFreeContext(optCon);

    // Finalise the census parameters.
    finiteness = regina::NBoolSet(! argIdeal, ! argFinite);
    orientability = regina::NBoolSet(! argNor, ! argOr);

    if (minimalPrimeP2)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL_PRIME |
            regina::NCensus::PURGE_P2_REDUCIBLE;
    else if (minimalPrime)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL_PRIME;
    else if (minimal)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL;
    else
        whichPurge = 0;

    return 0;
}

/**
 * Generic helper routine.
 *
 * Write the given time in human-readable form to the given output stream.
 */
void writeTime(std::ostream& out, long seconds) {
    bool started = false;
    if (seconds >= DAY_SEC) {
        out << (seconds / DAY_SEC) << " days ";
        seconds = seconds % DAY_SEC;
        started = true;
    }
    if (started || seconds >= HOUR_SEC) {
        out << (seconds / HOUR_SEC) << " hrs ";
        seconds = seconds % HOUR_SEC;
        started = true;
    }
    if (started || seconds >= MIN_SEC) {
        out << (seconds / MIN_SEC) << " min ";
        seconds = seconds % MIN_SEC;
        started = true;
    }
    out << seconds << " sec";
}

/**
 * Controller helper routine.
 *
 * Write the current date and time plus whitespace to the log stream and
 * return a reference to the log stream.
 */
std::ostream& ctrlLogStamp() {
    time_t t = time(0);
    std::string time = asctime(localtime(&t));
    // Remove the trailing newline.
    if (time[time.length() - 1] == '\n')
        time.resize(time.length() - 1);
    return logger << time << "  ";
}

/**
 * Controller helper routine.
 *
 * Return the next face pairing representation from the given input stream,
 * or an empty string if end-of-stream is reached.
 */
std::string ctrlNextPairing(std::istream& input) {
    std::string pairing;
    while (true) {
        std::getline(input, pairing);

        if (pairing.length() > 0) {
            // We have a real pairing.
            return pairing;
        }

        // Is the input stream exhausted?
        if (input.eof())
            return std::string();
    }
}

/**
 * Controller helper routine.
 *
 * Return the next slave available for processing work.
 */
int ctrlWaitForSlave(bool runningSlavesOnly = false) {
    if ((! runningSlavesOnly) && (nRunningSlaves < nSlaves)) {
        // We can use a slave that hasn't been started yet.
        // It's probably number (nRunningSlaves + 1).
        if (slaveTask[nRunningSlaves + 1].pairing < 0)
            return nRunningSlaves + 1;

        // Sigh.  Better go look.
        for (int slave = 1; slave <= nSlaves; slave++)
            if (slaveTask[slave].pairing < 0)
                return slave;

        // Strange.
        // Shouldn't get here ever.
        std::cerr << "ERROR: Something is seriously wrong." << std::endl;
        std::cerr << "A free slave could not be found." << std::endl;
        controllerError = true;
        return 0;
    }

    // All slaves are currently working.  Wait for the next one to finish.
    long results[3];
    MPI_Status status;
    MPI_Recv(results, 3, MPI_LONG, MPI_ANY_SOURCE, TAG_RESULT, MPI_COMM_WORLD,
        &status);

    int slave = status.MPI_SOURCE;
    time_t totSec = time(0) - slaveTask[slave].start;

    if (results[0] != slaveTask[slave].pairing ||
            results[1] != slaveTask[slave].subtask) {
        ctrlLogStamp() << "ERROR: Mismatched task data for slave "
            << slave << "." << std::endl;
        controllerError = true;
    }

    if (results[2] < 0) {
        ctrlLogStamp() << "ERROR: Unable to run task (see next line)."
            << std::endl;
        controllerError = true;
    }

    ctrlLogStamp() << "Task [" << results[0];
    if (results[1] >= 0)
        logger << '-' << results[1];
    logger << " @ slave " << slave << "]:  ";

    if (results[2] >= 0)
        logger << results[2] << " found";
    else
        logger << "UNABLE TO RUN TASK";

    logger << ", time " << totSec
        << "s (";
    writeTime(logger, totSec);
    logger << ")" << std::endl;

    if (results[2] >= 0)
        totTri += results[2];

    // Tidy up and return our answer.
    slaveTask[slave].pairing = slaveTask[slave].subtask = -1;
    nRunningSlaves--;
    return slave;
}

/**
 * Controller helper routine.
 *
 * Send the given face pairing to the next available slave for processing.
 */
void ctrlFarmPairing(const std::string& pairingRep) {
    taskID[2] = pairingRep.length();

    int slave = ctrlWaitForSlave();

    nRunningSlaves++;
    ctrlLogStamp() << "Farmed pairing " << taskID[0]
        << " --> slave " << slave << " ..." << std::endl;

    slaveTask[slave].pairing = taskID[0];
    slaveTask[slave].subtask = -1;
    slaveTask[slave].start = time(0);

    MPI_Send(taskID, 3, MPI_LONG, slave, TAG_REQUEST_TASK, MPI_COMM_WORLD);
    MPI_Send(const_cast<char*>(pairingRep.c_str()), taskID[2] + 1, MPI_CHAR,
        slave, TAG_REQUEST_PAIRING, MPI_COMM_WORLD);
}

/**
 * Controller helper routine.
 *
 * Send the given partial search to the next available slave for processing.
 */
template <class CensusType>
void ctrlFarmPartialSearch(
        const typename CensusType::GluingPermSearcher* search, void*) {
    if (! search) {
        // That's it for this face pairing.
        ctrlLogStamp() << "Pairing " << taskID[0] << ": Farmed "
            << taskID[1] << " subsearch(es) in total." << std::endl;
        return;
    }

    // We have a real subsearch.
    taskID[1]++;

    std::ostringstream searchRep;
    search->dumpTaggedData(searchRep);

    taskID[2] = searchRep.str().length();

    int slave = ctrlWaitForSlave();

    nRunningSlaves++;
    ctrlLogStamp() << "Farmed subsearch " << taskID[0] << '-' << taskID[1]
        << " --> slave " << slave << " ..." << std::endl;

    slaveTask[slave].pairing = taskID[0];
    slaveTask[slave].subtask = taskID[1];
    slaveTask[slave].start = time(0);

    MPI_Send(taskID, 3, MPI_LONG, slave, TAG_REQUEST_TASK, MPI_COMM_WORLD);
    MPI_Send(const_cast<char*>(searchRep.str().c_str()), taskID[2] + 1,
        MPI_CHAR, slave, TAG_REQUEST_SUBSEARCH, MPI_COMM_WORLD);
}

/**
 * Controller helper routine.
 *
 * Inform the given slave that it has no more tasks to process, and
 * that it may therefore exit.
 */
void ctrlStopSlave(int slave) {
    long msg[3] = { -1, -1, -1 };
    MPI_Send(&msg, 3, MPI_LONG, slave, TAG_REQUEST_TASK, MPI_COMM_WORLD);

    ctrlLogStamp() << "Slave " << slave << " stopped." << std::endl;
}

/**
 * Main routine for the controller.
 */
template <class CensusType>
int mainController() {
    // Prepare to read in the face pairings.
    std::ifstream input(pairsFile.c_str());
    if (! input) {
        std::cerr << "Could not open pairs file " << pairsFile
            << " for reading.\n";
        return 1;
    }

    // Start logging.
    std::string logFile = outputStub + ".log";
    logger.open(logFile.c_str());
    if (! logger) {
        std::cerr << "Could not open log file " << logFile << " for writing.\n";
        return 1;
    }

    // Go do it.
    slaveTask = new Task[nSlaves + 1];
    int slave;
    for (slave = 1; slave <= nSlaves; slave++)
        slaveTask[slave].pairing = slaveTask[slave].subtask = -1;
    nRunningSlaves = 0;
    totTri = 0;

    std::string pairingRep;
    if (depth > 0) {
        // Generate the face pairings and prepare subsearches.
        typename CensusType::Pairing* pairing;
        typename CensusType::GluingPermSearcher* searcher;
        while (! (pairingRep = ctrlNextPairing(input)).empty()) {
            taskID[0]++;
            taskID[1] = 0;

            pairing = CensusType::Pairing::fromTextRep(pairingRep);
            if (! pairing) {
                ctrlLogStamp() << "ERROR: Pairing " << taskID[0]
                    << " is invalid: " << pairingRep << std::endl;
                controllerError = true;
                continue;
            }
            if (! pairing->isCanonical()) {
                ctrlLogStamp() << "ERROR: Pairing " << taskID[0]
                    << " is not canonical: " << pairingRep << std::endl;
                controllerError = true;
                continue;
            }

            searcher = CensusType::bestSearcher(pairing,
                ! orientability.hasFalse(), ! finiteness.hasFalse(),
                whichPurge);
            searcher->runSearch(depth);
            delete searcher;
        }
    } else {
        // Just farm out the face pairing strings.
        while (! (pairingRep = ctrlNextPairing(input)).empty()) {
            taskID[0]++;
            taskID[1] = -1;

            ctrlFarmPairing(pairingRep);
        }
    }

    // Wait for everyone else to finish their tasks.
    while (nRunningSlaves > 0)
        ctrlWaitForSlave(true);

    // Stop all the slaves and finish!
    for (slave = 1; slave <= nSlaves; slave++)
        ctrlStopSlave(slave);

    delete[] slaveTask;

    // All done!
    input.close();
    ctrlLogStamp() << "Done: " << taskID[0] << " pairing(s) read" << std::endl;
    ctrlLogStamp() << "Done: " << totTri
        << " triangulation(s) found" << std::endl;
    if (controllerError)
        ctrlLogStamp() << "Done: ONE OR MORE ERRORS OCCURRED!" << std::endl;
    logger.close();
    return 0;
}

/**
 * Slave helper routine.
 *
 * Called each time the slave finds a complete triangulation.
 */
template <class CensusType>
void slaveFoundGluingPerms(const typename CensusType::GluingPermSearcher* perms,
        void* container) {
    if (perms) {
        typename CensusType::Triangulation* tri = perms->triangulate();

        bool ok = true;
        if (! tri->isValid())
            ok = false;
        else if ((! finiteness.hasFalse()) && tri->isIdeal())
            ok = false;
        else if ((! finiteness.hasTrue()) && (! tri->isIdeal()))
            ok = false;
        else if ((! orientability.hasTrue()) && tri->isOrientable())
            ok = false;
        else if ((minimal || minimalPrime || minimalPrimeP2) &&
                ! CensusType::mightBeMinimal(tri))
            ok = false;

        if (ok) {
            // Put it in the census!
            if (sigs) {
                if (! nSolns) {
                    sigStream.open(sigFile.c_str());
                    if (! sigStream)
                        sigStreamErr = true;
                }
                if (! sigStreamErr)
                    sigStream << tri->isoSig() << std::endl;
                delete tri;
            } else {
                regina::NPacket* dest =
                    static_cast<regina::NPacket*>(container);

                std::ostringstream out;
                out << "Item " << (nSolns + 1);
                tri->setPacketLabel(out.str());

                dest->insertChildLast(tri);
            }
            nSolns++;
        } else {
            // The fish that John West reject.
            delete tri;
        }
    }
}

/**
 * Slave helper routine.
 *
 * Write global census parameters in human-readable form to the given
 * outupt stream.
 */
void slaveDescribeCensusParameters(std::ostream& out) {
    if (dim4)
        out << "Searching for 4-manifold triangulations\n";
    else if (dim2)
        out << "Searching for 2-manifold triangulations\n";
    else
        out << "Searching for 3-manifold triangulations\n";

    if (cycleDecomp)
        out << "Searching for manifold decompositions\n";

    if (finiteness == regina::NBoolSet::sTrue)
        out << "Finite only\n";
    else if (finiteness == regina::NBoolSet::sFalse)
        out << "Ideal only\n";
    else
        out << "Finite and ideal\n";

    if (orientability == regina::NBoolSet::sTrue)
        out << "Orientable only\n";
    else if (orientability == regina::NBoolSet::sFalse)
        out << "Non-orientable only\n";
    else
        out << "Orientable and non-orientable\n";

    if (minimalPrimeP2)
        out << "Ignored obviously non-minimal, non-prime, "
            << "disc-reducible and/or P2-reducible triangulations\n";
    else if (minimalPrime)
        out << "Ignored obviously non-minimal, non-prime and/or "
            << "disc-reducible triangulations\n";
    else if (minimal)
        out << "Ignored obviously non-minimal triangulations\n";
}

/**
 * Slave helper routine.
 *
 * Construct a skeleton packet tree for the census data file, for use
 * when processing a face pairing subsearch.
 *
 * The census container will not be included, and should be inserted
 * as the last child of the parent.
 */
template <class CensusType>
regina::NPacket* slaveSkeletonTree(
        const typename CensusType::GluingPermSearcher* search,
        const char* searchRep) {
    // Create the overall parent packet.
    regina::NContainer* parent = new regina::NContainer();
    parent->setPacketLabel("Partial MPI census");

    // Create a child packet that describes the census parameters.
    regina::NText* desc = new regina::NText();
    desc->setPacketLabel("Parameters");
    std::ostringstream descStream;

    descStream << "Processed a " << WORD_face << " pairing subsearch.\n\n"
        << WORD_Face << " pairing:\n"
        << CensusType::pairingFor(search)->toString() << "\n\nSubsearch:\n"
        << searchRep << "\n\n";

    slaveDescribeCensusParameters(descStream);

    desc->setText(descStream.str());
    parent->insertChildLast(desc);

    // All done!
    return parent;
}

/**
 * Slave helper routine.
 *
 * Construct a skeleton packet tree for the census data file, for use
 * when processing a single face pairing.
 *
 * The census container will not be included, and should be inserted
 * as the last child of the parent.
 */
template <class CensusType>
regina::NPacket* slaveSkeletonTree(
        const typename CensusType::Pairing* pairing) {
    // Create the overall parent packet.
    regina::NContainer* parent = new regina::NContainer();
    parent->setPacketLabel("Partial MPI census");

    // Create a child packet that describes the census parameters.
    regina::NText* desc = new regina::NText();
    desc->setPacketLabel("Parameters");
    std::ostringstream descStream;

    descStream << "Processed a single " << WORD_face << " pairing:\n"
        << pairing->toString() << "\n\n";

    slaveDescribeCensusParameters(descStream);

    desc->setText(descStream.str());
    parent->insertChildLast(desc);

    // All done!
    return parent;
}

/**
 * Slave helper routine.
 *
 * Creates a filename whose name is unique to the current task being
 * processed.  This name will include the face pairing ID, as well as the
 * subsearch ID if appropriate.
 *
 * The given suffix will be appended to the filename; generally this
 * will include an extension (such as ".rga").
 */
void slaveMakeTaskFilename(std::string& result, const char* suffix) {
    std::ostringstream s;
    s << outputStub << '_' << taskID[0];
    if (taskID[1] >= 0)
        s << '-' << taskID[1];
    s << suffix;

    result = s.str();
}

/**
 * Slave helper routine.
 *
 * Inform the controller that the current task was successfully completed.
 */
void slaveSendResult(long nTriangulations) {
    long data[3] = { taskID[0], taskID[1], nTriangulations };
    MPI_Send(data, 3, MPI_LONG, 0, TAG_RESULT, MPI_COMM_WORLD);
}

/**
 * Slave helper routine.
 *
 * Signal that a fatal error occurred whilst working on the current task.
 *
 * A message is written to a newly created error file (specific to this
 * task), and the controller is informed also.
 */
void slaveBail(const std::string& error) {
    // Write the error to file.
    std::string errFile;
    slaveMakeTaskFilename(errFile, ".err");

    std::ofstream out(errFile.c_str());
    if (out) {
        if (taskID[1] < 0)
            out << "ERROR (pairing " << taskID[0] << "):" << std::endl;
        else
            out << "ERROR (pairing " << taskID[0] << ", subsearch "
                << taskID[1] << "):" << std::endl;
        out << error << std::endl;
    }
    out.close();

    // Inform the controller (by returning a negative number of
    // triangulations).
    slaveSendResult(-1);
}

/**
 * Slave helper routine.
 *
 * Request full data for a face pairing subsearch and process that subsearch.
 *
 * A data file is written if and only if at least one triangulation is
 * found (to reduce clutter for large censuses).
 *
 * The controller is informed of the final number of triangulations.
 */
template <class CensusType>
void slaveProcessPartialSearch() {
    char* searchRep = new char[taskID[2] + 1];

    // Fetch the subsearch data.
    MPI_Status status;
    MPI_Recv(searchRep, taskID[2] + 1, MPI_CHAR, 0,
        TAG_REQUEST_SUBSEARCH, MPI_COMM_WORLD, &status);

    // Construct the subsearch.
    regina::NPacket* parent = 0;
    regina::NPacket* dest = 0;
    if (! sigs) {
        dest = new regina::NContainer();
        dest->setPacketLabel("Triangulations");
    }

    typename CensusType::GluingPermSearcher* search;
    {
        std::istringstream s(searchRep);
        search = CensusType::GluingPermSearcher::readTaggedData(s,
            slaveFoundGluingPerms<CensusType>, dest);
    }

    if (search == 0 || search->inputError()) {
        slaveBail(std::string("Invalid search data:\n") + searchRep);
        delete dest;
        return;
    }

    // Prepare a packet tree (or output file) to wrap around the search.
    if (sigs) {
        slaveMakeTaskFilename(sigFile, ".sig");
        sigStreamErr = false;
    } else {
        parent = slaveSkeletonTree<CensusType>(search, searchRep);
        parent->insertChildLast(dest);
    }

    // Run the partial census.
    nSolns = 0;
    if (! dryRun)
        search->runSearch();

    if (nSolns > 0) {
        // Write the completed census to file.
        // Use a unique filename for each task.
        if (sigs) {
            sigStream.close();
            if (sigStreamErr)
                slaveBail("Signature file could not be written.");
            else
                slaveSendResult(nSolns);
        } else {
            std::string outFile;
            slaveMakeTaskFilename(outFile, ".rga");

            if (regina::writeXMLFile(outFile.c_str(), parent))
                slaveSendResult(nSolns);
            else
                slaveBail("Output file could not be written.");
        }
    } else {
        // No triangulations.  Just inform the controller.
        slaveSendResult(0);
    }

    delete parent;
    delete const_cast<typename CensusType::Pairing*>(
        CensusType::pairingFor(search));
    delete search;
    delete[] searchRep;
}

/**
 * Slave helper routine.
 *
 * Request full data for a single face pairing and process that face pairing.
 *
 * A data file is written if and only if at least one triangulation is
 * found (to reduce clutter for large censuses).
 *
 * The controller is informed of the final number of triangulations.
 */
template <class CensusType>
void slaveProcessPairing() {
    char* pairingRep = new char[taskID[2] + 1];

    // Fetch the pairing itself.
    MPI_Status status;
    MPI_Recv(pairingRep, taskID[2] + 1, MPI_CHAR, 0,
        TAG_REQUEST_PAIRING, MPI_COMM_WORLD, &status);

    // Parse the face pairing.
    typename CensusType::Pairing* pairing =
        CensusType::Pairing::fromTextRep(pairingRep);
    if (! pairing) {
        if (dim4)
            slaveBail(std::string("Invalid facet pairing: ") + pairingRep);
        else if (dim2)
            slaveBail(std::string("Invalid edge pairing: ") + pairingRep);
        else
            slaveBail(std::string("Invalid face pairing: ") + pairingRep);
        return;
    }
    if (! pairing->isCanonical()) {
        if (dim4)
            slaveBail(std::string("Non-canonical facet pairing: ") +
                pairingRep);
        else if (dim2)
            slaveBail(std::string("Non-canonical edge pairing: ") +
                pairingRep);
        else
            slaveBail(std::string("Non-canonical face pairing: ") +
                pairingRep);
        return;
    }

    // Run the partial census.
    regina::NPacket* parent = 0;
    regina::NPacket* dest = 0;
    if (sigs) {
        slaveMakeTaskFilename(sigFile, ".sig");
        sigStreamErr = false;
    } else {
        parent = slaveSkeletonTree<CensusType>(pairing);
        dest = new regina::NContainer();
        dest->setPacketLabel("Triangulations");
        parent->insertChildLast(dest);
    }

    nSolns = 0;
    if (! dryRun)
        CensusType::findAllPerms(pairing, ! orientability.hasFalse(),
            ! finiteness.hasFalse(), whichPurge, dest);

    if (nSolns > 0) {
        // Write the completed census to file.
        // Use a unique filename for each task.
        if (sigs) {
            sigStream.close();
            if (sigStreamErr)
                slaveBail("Signature file could not be written.");
            else
                slaveSendResult(nSolns);
        } else {
            std::string outFile;
            slaveMakeTaskFilename(outFile, ".rga");

            if (regina::writeXMLFile(outFile.c_str(), parent))
                slaveSendResult(nSolns);
            else
                slaveBail("Output file could not be written.");
        }
    } else {
        // No triangulations.  Just inform the controller.
        slaveSendResult(0);
    }

    delete parent;
    delete pairing;
    delete[] pairingRep;
}

/**
 * Main routine for a slave (ranks 1..size).
 */
template <class CensusType>
int mainSlave() {
    // Keep fetching and processing tasks until there are no more.
    MPI_Status status;
    while (true) {
        // Get identification data for the new processing task.
        // See the taskID declaration at the top of this file for what
        // each array entry means.
        MPI_Recv(taskID, 3, MPI_LONG, 0, TAG_REQUEST_TASK,
            MPI_COMM_WORLD, &status);

        if (taskID[0] < 0)
            break;

        if (taskID[1] < 0)
            slaveProcessPairing<CensusType>();
        else
            slaveProcessPartialSearch<CensusType>();
    }

    return 0;
}

/**
 * Main routine for all processors.
 *
 * Parse the command-line arguments for census options, then determine
 * whether we are controller or slave and run a specialised main routine
 * accordingly.
 */
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Which processor are we?
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Extract census options.
    int retVal = parseCmdLine(argc, (const char**)argv, (rank == 0));

    if (retVal == 0) {
        // No breakage yet.  Controller or slave?
        if (rank == 0) {
            // We're the controller.
            // How many processors in total?
            int size;
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            if (size <= 1) {
                std::cerr << "ERROR: At least two processors are required "
                    "(one controller and one slave).\n";
                retVal = 1;
            } else {
                nSlaves = size - 1;
                if (dim2)
                    retVal = mainController<Dim2Params>();
#if SUPPORT_DIM4
                else if (dim4)
                    retVal = mainController<Dim4Params>();
#endif
                else if (cycleDecomp)
                    retVal = mainController<Dim3DecompParams>();
                else
                    retVal = mainController<Dim3Params>();
            }
        } else {
            // We're one of many slaves.
            if (dim2)
                retVal = mainSlave<Dim2Params>();
#if SUPPORT_DIM4
            else if (dim4)
                retVal = mainSlave<Dim4Params>();
#endif
            else if (cycleDecomp)
                retVal = mainSlave<Dim3DecompParams>();
            else
                retVal = mainSlave<Dim3Params>();
        }
    }

    MPI_Finalize();
    return retVal;
}

