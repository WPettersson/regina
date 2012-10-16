
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Form a census of triangulations that satisfy given properties         *
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

#include <cctype>
#include <fstream>
#include <memory>
#include <sstream>
#include <popt.h>
#include <unistd.h>
#include "census/dim2census.h"
#include "census/dim2edgepairing.h"
#include "census/dim2gluingpermsearcher.h"
#if SUPPORT_DIM4
#include "census/dim4census.h"
#include "census/dim4facetpairing.h"
#include "census/dim4gluingpermsearcher.h"
#endif
#include "census/ncensus.h"
#include "census/nfacepairing.h"
#include "census/ngluingpermsearcher.h"
#include "dim2/dim2triangulation.h"
#if SUPPORT_DIM4
#include "dim4/dim4triangulation.h"
#endif
#include "file/nxmlfile.h"
#include "packet/ncontainer.h"
#include "packet/ntext.h"
#include "triangulation/ntriangulation.h"

#include "cycledecomp.h"

#define WORD_face (dim4 ? "facet" : dim2 ? "edge" : "face")
#define WORD_Face (dim4 ? "Facet" : dim2 ? "Edge" : "Face")
#define WORD_faces (dim4 ? "facets" : dim2 ? "edges" : "faces")
#define WORD_tetrahedra (dim4 ? "pentachora" : dim2 ? "triangles" : "tetrahedra")
#define WORD_tetrahedron (dim4 ? "pentachoron" : dim2 ? "triangle" : "tetrahedron")
#define WORD_Tetrahedron (dim4 ? "Pentachoron" : dim2 ? "Triangle" : "Tetrahedron")

// Constants.
const long MAXTET = 15;
const unsigned SLEEP_SECONDS = 1;

// Census parameters.
long nTet = 0, nBdryFaces = -1;
regina::NBoolSet
    finiteness(true, true),
    orientability(true, true),
    boundary(true, true);
int minimal = 0;
int minimalPrime = 0;
int minimalPrimeP2 = 0;
int cycleDecomp = 0;
int dim2 = 0;
int dim4 = 0;
int usePairs = 0;
int sigs = 0;
int whichPurge = 0;
int genPairs = 0;
std::string outFile;

// Variables used for a dump of face pairings.
std::auto_ptr<std::ostream> dumpStream;
unsigned long totPairings = 0;

// Variables used for output.
long long nSolns;
std::ofstream sigStream;

// Forward declarations:
template <class CensusType>
void foundGluingPerms(const typename CensusType::GluingPermSearcher*, void*);

template <class CensusType>
void foundCycleDecomp(const typename CensusType::GluingPermSearcher*, void*); 

template <class CensusType>
int runCensus();

// Differences between censuses of 2, 3 and 4-manifolds:
struct Dim2Params {
    typedef regina::Dim2EdgePairing Pairing;
    typedef regina::Dim2GluingPermSearcher GluingPermSearcher;
    typedef regina::Dim2Triangulation Triangulation;

    inline static void findAllPerms(const Pairing* p,
            const Pairing::IsoList* autos, bool orientableOnly,
            bool finiteOnly, int whichPurge, regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, autos,
            orientableOnly,
            foundGluingPerms<Dim2Params>, dest);
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

    inline static void findAllPerms(const Pairing* p,
            const Pairing::IsoList* autos, bool orientableOnly,
            bool finiteOnly, int whichPurge, regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, autos,
            orientableOnly, finiteOnly, whichPurge,
            foundGluingPerms<Dim3Params>, dest);
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
    inline static void findAllPerms(const Pairing* p,
            const Pairing::IsoList* autos, bool orientableOnly,
            bool finiteOnly, int whichPurge, regina::NPacket* dest) {
            CycleDecompSearcher *searcher = new CycleDecompSearcher(p, autos,
            orientableOnly, foundCycleDecomp<Dim3DecompParams>, dest, true);
            searcher->runSearch();
            delete searcher;
    }

    inline static bool mightBeMinimal(Triangulation* tri) {
        return true;
    }
};

#if SUPPORT_DIM4
struct Dim4Params {
    typedef regina::Dim4FacetPairing Pairing;
    typedef regina::Dim4GluingPermSearcher GluingPermSearcher;
    typedef regina::Dim4Triangulation Triangulation;

    inline static void findAllPerms(const Pairing* p,
            const Pairing::IsoList* autos, bool orientableOnly,
            bool finiteOnly, int /* whichPurge */, regina::NPacket* dest) {
        GluingPermSearcher::findAllPerms(p, autos,
            orientableOnly, finiteOnly,
            foundGluingPerms<Dim4Params>, dest);
    }

    inline static bool mightBeMinimal(Triangulation*) {
        return true;
    }

    inline static const Pairing* pairingFor(const GluingPermSearcher* s) {
        return s->getFacetPairing();
    }
};
#endif

/**
 * What to do with each complete triangulation that is generated.
 */
template <class CensusType>
void foundGluingPerms(const typename CensusType::GluingPermSearcher* perms,
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
                sigStream << tri->isoSig() << std::endl;
                //std::cout << tri->isoSig() << std::endl;;
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

template <class CensusType>
void foundCycleDecomp(const typename CensusType::GluingPermSearcher* perms,
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
                sigStream << tri->isoSig() << std::endl;
                //std::cout << tri->isoSig();
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
        //std::cout << std::endl;
    }
}

/**
 * What to do with each face/facet pairing that is generated.
 */
template <class CensusType>
void foundFacePairing(const typename CensusType::Pairing* pairing,
        const typename CensusType::Pairing::IsoList* autos, void* container) {
    long long before = nSolns;
    if (pairing) {
        std::cout << pairing->toString() << std::endl;

        CensusType::findAllPerms(pairing, autos,
            ! orientability.hasFalse(), ! finiteness.hasFalse(),
            whichPurge, static_cast<regina::NPacket*>(container));
        std::cout << "Found " << (nSolns - before) << " triangulations." << std::endl;
    }

}

/**
 * Dump the given face/facet pairing to dumpStream.
 */
template <class CensusType>
void dumpPairing(const typename CensusType::Pairing* pair,
        const typename CensusType::Pairing::IsoList*, void*) {
    if (pair) {
        if (dumpStream.get())
            (*dumpStream) << (*pair).toTextRep() << std::endl;
        else
            std::cout << (*pair).toTextRep() << std::endl;
        totPairings++;
    }
}

/**
 * Return a new text packet storing the census parameters.
 */
regina::NText* parameterPacket() {
    regina::NText* desc = new regina::NText();
    desc->setPacketLabel("Parameters");
    std::ostringstream descStream;

    if (dim4)
        descStream << "Searching for 4-manifold triangulations\n";
    else if (dim2)
        descStream << "Searching for 2-manifold triangulations\n";
    else
        descStream << "Searching for 3-manifold triangulations\n";

    if (usePairs)
        descStream << "Only used a subset of all available "
            << WORD_face << " pairings\n";
    else {
        descStream << nTet << ' ' <<
            (nTet == 1 ? WORD_tetrahedron : WORD_tetrahedra) << "\n";

        if (boundary == regina::NBoolSet::sTrue)
            descStream << "Boundary " << WORD_faces << " only\n";
        else if (boundary == regina::NBoolSet::sFalse)
            descStream << "No boundary " << WORD_faces << " only\n";
        else
            descStream << "With and without boundary " << WORD_faces << "\n";

        if (nBdryFaces >= 0)
            descStream << "Requires precisely " << nBdryFaces <<
                " boundary " << WORD_faces << "\n";
    }

    if (finiteness == regina::NBoolSet::sTrue)
        descStream << "Finite only\n";
    else if (finiteness == regina::NBoolSet::sFalse)
        descStream << "Ideal only\n";
    else
        descStream << "Finite and ideal\n";

    if (orientability == regina::NBoolSet::sTrue)
        descStream << "Orientable only\n";
    else if (orientability == regina::NBoolSet::sFalse)
        descStream << "Non-orientable only\n";
    else
        descStream << "Orientable and non-orientable\n";

    if (minimalPrimeP2)
        descStream << "Ignored obviously non-minimal, non-prime, "
            << "disc-reducible and/or P2-reducible triangulations\n";
    else if (minimalPrime)
        descStream << "Ignored obviously non-minimal, non-prime and/or "
            << "disc-reducible triangulations\n";
    else if (minimal)
        descStream << "Ignored obviously non-minimal triangulations\n";

    desc->setText(descStream.str());
    return desc;
}

/**
 * Parse the command-line arguments and then farm out to runCensus(),
 * which does the real work.
 */
int main(int argc, const char* argv[]) {
    // Set up the command-line arguments.
    int argBdry = 0;
    int argNoBdry = 0;
    int argOr = 0;
    int argNor = 0;
    int argFinite = 0;
    int argIdeal = 0;
    poptOption opts[] = {
        { "tetrahedra", 't', POPT_ARG_LONG, &nTet, 0,
            "Number of tetrahedra.", "<tetrahedra>" },
        { "boundary", 'b', POPT_ARG_NONE, &argBdry, 0,
            "Must have at least one boundary face.", 0 },
        { "internal", 'i', POPT_ARG_NONE, &argNoBdry, 0,
            "Must have all faces internal (no boundary faces).", 0 },
        { "bdryfaces", 'B', POPT_ARG_LONG, &nBdryFaces, 0,
            "Must have fixed number (>= 1) of boundary faces.", "<faces>" },
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
            "Ignore obviously non-minimal, non-prime and/or disc-reducible triangulations.", 0 },
        { "minprimep2", 'N', POPT_ARG_NONE, &minimalPrimeP2, 0,
            "Ignore obviously non-minimal, non-prime, disc-reducible and/or P2-reducible triangulations.", 0 },
        { "cycledecomp" , 'c', POPT_ARG_NONE, &cycleDecomp, 0,
            "Find permutations by cycle decompositions of the face pairing graph.", 0 },
        { "dim2", '2', POPT_ARG_NONE, &dim2, 0,
            "Run a census of 2-manifold triangulations, "
            "not 3-manifold triangulations.  Here --tetrahedra counts "
            "triangles, and --bdryfaces counts boundary edges.", 0 },
#if SUPPORT_DIM4
        { "dim4", '4', POPT_ARG_NONE, &dim4, 0,
            "Run a census of 4-manifold triangulations, "
            "not 3-manifold triangulations.  Here --tetrahedra counts "
            "pentachora, and --bdryfaces counts boundary facets.", 0 },
#endif
        { "sigs", 's', POPT_ARG_NONE, &sigs, 0,
            "Write isomorphism signatures only, not full Regina data files.",
            0 },
        { "genpairs", 'p', POPT_ARG_NONE, &genPairs, 0,
            "Only generate face pairings, not triangulations.", 0 },
        { "usepairs", 'P', POPT_ARG_NONE, &usePairs, 0,
            "Only use face pairings read from standard input.", 0 },
        POPT_AUTOHELP
        { 0, 0, 0, 0, 0, 0, 0 }
    };

    poptContext optCon = poptGetContext(0, argc, argv, opts, 0);
    poptSetOtherOptionHelp(optCon, "<output-file>");

    // Parse the command-line arguments.
    int rc = poptGetNextOpt(optCon);
    if (rc != -1) {
        std::cerr << poptBadOption(optCon, POPT_BADOPTION_NOALIAS)
            << ": " << poptStrerror(rc) << "\n\n";
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return 1;
    }

    const char** otherOpts = poptGetArgs(optCon);
    if (otherOpts && otherOpts[0]) {
        outFile = otherOpts[0];
        if (otherOpts[1]) {
            std::cerr << "Too many arguments.\n\n";
            poptPrintHelp(optCon, stderr, 0);
            poptFreeContext(optCon);
            return 1;
        }
    }

    // Run a sanity check on the command-line arguments.
    bool broken = false;
    if ((! genPairs) && outFile.empty()) {
        std::cerr << "An output file must be specified.\n";
        broken = true;
    } else if (genPairs && (argOr || argNor)) {
        std::cerr << "Orientability options cannot be used with "
            << "-p/--genpairs.\n";
        broken = true;
    } else if (genPairs && (argFinite || argIdeal)) {
        std::cerr << "Finiteness options cannot be used with -p/--genpairs.\n";
        broken = true;
    } else if (genPairs && (minimal || minimalPrime || minimalPrimeP2)) {
        std::cerr << "Minimality options cannot be used with -p/--genpairs.\n";
        broken = true;
    } else if (genPairs && sigs) {
        std::cerr << "Signature output cannot be used with -p/--genpairs.\n";
        broken = true;
    } else if (dim2 && dim4) {
        std::cerr << "Options -2/--dim2 and -4/--dim4 cannot be used together.\n";
        broken = true;
    } else if (dim2 && (argFinite || argIdeal)) {
        std::cerr << "Finiteness options cannot be used with -2/--dim2.\n";
        broken = true;
    } else if (dim2 && (minimalPrime || minimalPrimeP2)) {
        std::cerr << "Primeness options cannot be used with -2/--dim2 "
            "(the weaker -m/--minimal can).\n";
        broken = true;
    } else if (dim4 && (minimal || minimalPrime || minimalPrimeP2)) {
        std::cerr << "Minimality options cannot be used with -4/--dim4.\n";
        broken = true;
    } else if (usePairs && nTet) {
        std::cerr << WORD_Tetrahedron
            << " options cannot be used with -P/--usepairs.\n";
        broken = true;
    } else if (usePairs && (argBdry || argNoBdry || (nBdryFaces != -1))) {
        std::cerr << "Boundary options cannot be used with -P/--usepairs.\n";
        broken = true;
    } else if ((! usePairs) && nTet == 0) {
        std::cerr << "The number of " << WORD_tetrahedra
            << " must be specified using option -t/--tetrahedra.\n";
        broken = true;
    } else if ((! usePairs) && (nTet < 1 || nTet > MAXTET)) {
        std::cerr << "The number of " << WORD_tetrahedra
            << " must be between 1 and " << MAXTET << " inclusive.\n";
        broken = true;
    } else if (argBdry && argNoBdry) {
        std::cerr << "Options -b/--boundary and -i/--internal "
            << "cannot be used together.\n";
        broken = true;
    } else if (argOr && argNor) {
        std::cerr << "Options -o/--orientable and -n/--nonorientable "
            << "cannot be used together.\n";
        broken = true;
    } else if (argFinite && argIdeal) {
        std::cerr << "Options -f/--finite and -d/--ideal "
            << "cannot be used together.\n";
        broken = true;
    } else if (genPairs && usePairs) {
        std::cerr << "Options -p/--genpairs and -P/--usepairs "
            << "cannot be used together.\n";
        broken = true;
    } else if (cycleDecomp && (dim2 || dim4)) {
        std::cerr << "Cycle decompositions are only available in 3 "
            << "dimensions.\n";
        broken = true;
    } else if (cycleDecomp && (argBdry)) {
        std::cerr << "Cycle decompositions cannot be used with "
            << "boundary facets (-b/--boundary).\n";
        broken = true;
    }

    if ((! broken) && (nBdryFaces != -1)) {
        if (nBdryFaces < 0) {
            std::cerr << "Number of boundary "
                << WORD_faces << " cannot be negative.\n";
            broken = true;
        } else if (nBdryFaces == 0) {
            // Asking for no boundary faces.
            if (argBdry) {
                std::cerr << "Option -b/--boundary cannot be used with "
                    << "-B/--bdryfaces=0.\n";
                broken = true;
            } else
                argNoBdry = 1;
        } else if ((! dim4) && (! dim2) && (nBdryFaces % 2 != 0)) {
            std::cerr << "Number of boundary faces must be even.\n";
            broken = true;
        } else if (dim2 && ((nTet + nBdryFaces) % 2 != 0)) {
            std::cerr << "Number of boundary edges must have the "
                "same parity as the number of triangles.\n";
            broken = true;
        } else if (dim4 && ((nTet + nBdryFaces) % 2 != 0)) {
            std::cerr << "Number of boundary facets must have the "
                "same parity as the number of pentachora.\n";
            broken = true;
        } else if ((! dim4) && (! dim2) && (nBdryFaces > 2 * nTet + 2)) {
            std::cerr << "Number of boundary faces for " << nTet
                << (nTet == 1 ? " tetrahedron" : " tetrahedra")
                << " can be at most " << (2 * nTet + 2) << ".\n";
            broken = true;
        } else if (dim2 && (3 * nTet - nBdryFaces < 2 * (nTet - 1))) {
            std::cerr << "Number of boundary edges for " << nTet
                << (nTet == 1 ? " triangle" : " triangles")
                << " can be at most " << (nTet + 2) << ".\n";
            broken = true;
        } else if (dim4 && (5 * nTet - nBdryFaces < 2 * (nTet - 1))) {
            std::cerr << "Number of boundary facets for " << nTet
                << (nTet == 1 ? " pentachoron" : " pentachora")
                << " can be at most " << (3 * nTet + 2) << ".\n";
            broken = true;
        } else {
            // Asking for a valid positive number of boundary faces.
            if (argNoBdry) {
                std::cerr << "Option -i/--internal cannot be used with "
                    << "-B/--bdryfaces=<non-zero>.\n";
                broken = true;
            } else
                argBdry = 1;
        }
    }

    if (broken) {
        std::cerr << '\n';
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return 1;
    }

    // Done parsing the command line.
    poptFreeContext(optCon);

    // Finalise the census parameters.
    finiteness = regina::NBoolSet(! argIdeal, ! argFinite);
    orientability = regina::NBoolSet(! argNor, ! argOr);
    boundary = regina::NBoolSet(! argNoBdry, ! argBdry);

    // And off we go!
    if (dim2)
        return runCensus<Dim2Params>();
#if SUPPORT_DIM4
    else if (dim4)
        return runCensus<Dim4Params>();
#endif
    else if (cycleDecomp)
        return runCensus<Dim3DecompParams>();
    else
        return runCensus<Dim3Params>();
}

/**
 * The main census generation routine.
 */
template <class CensusType>
int runCensus() {
    // Are we only dumping face pairings?
    if (genPairs) {
        if (! outFile.empty()) {
            dumpStream.reset(new std::ofstream(outFile.c_str()));
            if ((! dumpStream.get()) || (! *dumpStream)) {
                std::cerr << "Could not write to file " << outFile << ".\n";
                return 1;
            }
        }
        CensusType::Pairing::findAllPairings(nTet, boundary,
            nBdryFaces, dumpPairing<CensusType>, 0, false);
        std::cerr << "Total " << WORD_face << " pairings: "
            << totPairings << std::endl;
        return 0;
    }

    // We're actually generating triangulations.
    nSolns = 0;

    // Prepare the packet tree (or signature file) for output.
    regina::NPacket* parent = 0;
    regina::NPacket* census = 0;
    regina::NPacket* desc = 0;
    if (sigs) {
        sigStream.open(outFile.c_str());
        if (! sigStream) {
            std::cerr << "Signature file " << outFile
                << " could not be written.\n";
            return 1;
        }
    } else {
        parent = new regina::NContainer();
        if (usePairs)
            parent->setPacketLabel("Partial command-line census");
        else
            parent->setPacketLabel("Command-line census");

        desc = parameterPacket();

        census = new regina::NContainer();
        census->setPacketLabel("Triangulations");

        parent->insertChildLast(desc);
        parent->insertChildLast(census);
    }

    // Start the census running.
    if (minimalPrimeP2)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL_PRIME |
            regina::NCensus::PURGE_P2_REDUCIBLE;
    else if (minimalPrime)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL_PRIME;
    else if (minimal)
        whichPurge = regina::NCensus::PURGE_NON_MINIMAL;

    if (usePairs) {
        // Only use the face pairings read from standard input.
        std::cout << "Trying " << WORD_face << " pairings..." << std::endl;
        std::string pairingList = WORD_Face;
        pairingList += " pairings:\n\n";

        std::string pairingRep;
        while (true) {
            std::getline(std::cin, pairingRep);

            if (pairingRep.length() > 0) {
                typename CensusType::Pairing* pairing =
                    CensusType::Pairing::fromTextRep(pairingRep);
                if (! pairing) {
                    std::cerr << "Invalid " << WORD_face << " pairing: "
                        << pairingRep << std::endl;
                    pairingList += "INVALID: ";
                    pairingList += pairingRep;
                    pairingList += '\n';
                } else if (! pairing->isCanonical()) {
                    std::cerr << "Non-canonical " << WORD_face
                        << " pairing: " << pairingRep << std::endl;
                    pairingList += "NON-CANONICAL: ";
                    pairingList += pairingRep;
                    pairingList += '\n';
                } else {
                    std::cout << pairing->toString() << std::endl;
                    foundFacePairing<CensusType>(pairing, 0, census);
                    pairingList += pairing->toString();
                    pairingList += '\n';
                }
                delete pairing;
            }

            if (std::cin.eof())
                break;
        }

        // Store the face pairings used with the census.
        if (! sigs) {
            regina::NText* pairingPacket = new regina::NText(pairingList);
            pairingPacket->setPacketLabel(
                dim4 ? "Facet Pairings" : dim2 ? "Edge Pairings" :
                "Face Pairings");
            parent->insertChildAfter(pairingPacket, desc);
        }
    } else {
        // An ordinary all-face-pairings census.
        std::cout << "Starting census generation..." << std::endl;
        CensusType::Pairing::findAllPairings(nTet, boundary, nBdryFaces,
            foundFacePairing<CensusType>, census /* dest */);

        std::cout << "Finished." << std::endl;
    }

    // Write the completed census to file.
    if (sigs) {
        sigStream.close();
    } else {
        if (! regina::writeXMLFile(outFile.c_str(), parent)) {
            std::cerr << "Output file " << outFile
                << " could not be written.\n";
            return 1;
        }
    }

    std::cout << "Total triangulations: " << nSolns << std::endl;
    return 0;
}

