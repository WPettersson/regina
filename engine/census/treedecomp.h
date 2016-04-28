
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

/*! \file census/treedecomp.h
 *  \brief Supports searching through all possible sets of tetrahedron
 *  gluing permutations for a given tetrahedron face pairing.
 */

#ifndef __TREEDECOMP_H
#ifndef __DOXYGEN
#define __TREEDECOMP_H
#endif

#include "regina-core.h"
#include "census/ngluingpermsearcher.h"

namespace regina {

/**
 * \weakgroup census
 * @{
 */

/**
 * A gluing permutation search class that
 *
 * The search algorithm uses
 * For details see
 *
 * \ifacespython Not present.
 */


// TODO Inheriting from NGluingPermSearcher so we can interface with regina
// census code easily (both library based and tricensus) but we barely use any
// of the NGluingPermSearcher-specific members.
class REGINA_API TreeDecompSearcher : public NGluingPermSearcher {

    // 2000 is extreme, but this comment exists as a warning that the TFE/TVE
    // objects will overrun at around 2730 tetrahedra.
    static const int maxTetrahedra = 2000;

    // Represents an edge on a face of a tetrahedron. Used to give each edge on
    // a boundary face a designation. These are used to pair up such edge-face
    // pairs, itself useful for tracking the links of the edges of the
    // triangulation.
    typedef uint16 TFE;
    TFE inline TFE_(int t, int f, int e) { return (t*4 + f)*6 + e; } ;
    // 4 faces per tet, 6 edges per tet. TODO work with edges per face?

    // Represents an edge which lies on the given tetrahedron but not incident to
    // the given vertex. Used to give the boundary of the link of each vertex a
    // designation, where the boundary edge of a link has a given TVE as a name
    // if it belongs to the link of the given vertex on the given tetrahedron,
    // and is parallel to that edge.
    typedef uint16 TVE;
    TVE inline TVE_(int t, int v, int e) { return (t*4 + v)*6 + e; } ;
    // 4 vertices per tet, 6 edges per tet. TODO work with edges per vertex?

    // Represents a vertex on a tetrahedron. Used to track which such pairs are
    // part of the same triangulation-vertex.
    typedef uint16 TV;
    TV inline TV_(int t, int v) { return t*4 + v; } ;

    // An arc in the face pairing graph represents a pair of "faces" being
    // identified. We actually store each "face" as a pair of face/tet ints.
    struct Arc {
        int t1,t2;
        int f1,f2;
    };

    // The configuration of the boundary of a partial triangulation. This can
    // be constructed either from an actual triangulation, or just from
    // building upon a configuration of a "smaller" configuration. This second
    // fact is the basis of the FPT algorithm for admissibility testing.
    class Config {
        // Bare constructor
        Config();
        // Copy constructor
        Config(const Config c);

        // For each TV, we want to easily be able to find the set of equivalent
        // TVs. That is, which TVs are actually part of the same
        // triangulation-vertex.
        std::map<TV, std::set<TV>*> equivMap;


        // Combine another configuration with this one. The two configurations
        // must have distinct tetrahedra.
        void mergeWith(const Config &other);

        // Undo such a merge. Since the two configurations had distinct
        // tetrahedra, we just need to delete any pieces of information
        // relating to the "new" tetrahedra
        void undoMerge(std::set<int> tets);

        // Is the given face on the given tetrahedra on the boundary?
        bool onBoundary(int t, int f);

        // Glue two faces together. The two faces are given by t1,f1 and t2,f2,
        // and the exact gluing is given by gluing (an entry into
        // NEdge::some_table TODO
        bool glue(int gluing, int t1, int f1, int t2, int f2);

        // Unglue two faces. TODO: Might need more info passed in
        void unGlue(int t1, int f1);

        // Mark a and b as identified pairs. If orientation is true, then the
        // pair are matched such that lowest vertex of a meets lowest vertex of
        // b.
        addTFEPair(TFE a, TFE b, bool orientation);
    };

    // Represents a triangulation we are building. Note that we need a bit more
    // than an NTriangulation object. Still not set on NTriangulation member or
    // inherited?
    class Triangulation {
        // Maybe inherit from NTriangulation instead of having member?
        // Means we don't need to "extract" the NTriangulation object out
        NTriangulation t_;

        // Basic constructor
        Triangulation();
        // Copy constructor
        Triangulation(const Triangulation t);

        // Combine another triangulation with this one. The two triangulations
        // must have distinct tetrahedra.
        void mergeWith(const Triangulation other);

        // Undo such a merge. Since the two triangulations had distinct
        // tetrahedra, we just need to delete any pieces of information
        // relating to the "new" tetrahedra
        void undoMerge(std::set<int> tets);


        // Glue two faces together. The two faces are given by t1,f1 and t2,f2,
        // and the exact gluing is given by gluing (an entry into
        // NEdge::some_table TODO
        bool glue(int gluing, int t1, int f1, int t2, int f2);

        // Extract the NTriangulation object (if it's a member)
        const NTriangulation * object() const;
    };

    // Represents a bag of the tree decomposition
    class Bag {
        // The integer contents of this particular bag
        int[] contents_;
        // Which tetrahedra are we adding in this bag
        int[] toAdd_;

        // The list of arcs which will be "added" (where adding means
        // identifying related faces) in this bag
        Arc[] arcs;
        // How many arcs will be added at this bag
        int numArcs;

        // How many children do we have
        int numChildren;
        // The children of this bag
        Bag*[] children;

        // The possible configurations achievable at this bag
        std::list<Config*> possibleConfigs_;
        // The possible triangulations achievable at this bag
        std::list<Triangulation*> possibleTriangulations_;

        std::list<Triangulation *>::iterator triangulations();

        bool configsFound;
        bool trisFound;

        void findConfigs();
        void findTriangulations();


        Config::iterator childConfigs();
        bool hasNoConfigs();
        // When asked for next Config, start again at start
        bool resetConfigCount();
        // Does this bag have more Configs?
        bool hasNextConfig();
        // Get the next Config
        Config getNextConfig();

        bool addArcs(Config& c);

        void storeConfig(Config& c);

        Triangulation::iterator childTriangulations();
        bool hasNoTriangulations();
        bool resetTriangulationCount();
        bool hasNextTriangulation();

        bool hasValidConfig(Triangulation& t);

        void storeTri(Triangulation *t);
    };




    static inline int combine(int tet, int vert) { return 4*tet+vert; }
};

// Inline functions for TreeDecompSearcher


} // namespace regina

#endif

