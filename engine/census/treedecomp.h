
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
// of the NGluingPermSearcher-specific members. Should probably be changed to
// NCensus (which needs to be implemented, see
// github.com/WPettersson/regina/issues/8)
class REGINA_API TreeDecompSearcher : public NGluingPermSearcher {

    typedef const typename std::iterator<Triangulation*> TriIterator;

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
    class Arc {
        public:
            Arc(FacetSpec<3> one, FacetSpec<3> two);
            FacetSpec<3> one, two;
    };


    // For a given TFE, states which TFE (as a value) come "before" and "after"
    // this TFE. For each of those, there's also an "orientation" boolean.
    // If next() == b, and nextO() == false, then to continue around the link
    // one would next look at b->prev() as orientation is not consistent
    // between these two. However, if nextO() == true, then one should continue
    // with b->next()

    class LinkEdge {
        public:
            LinkEdge(TFE next, bool nextO, TFE prev, bool prevO);
            inline TFE next(bool o) { if (o) return next_; else return prev_; }
            inline TFE next() { return next_; }
            inline TFE prev() { return prev_; }
            inline bool nextO(bool o) { if (o) return nextO_; else return prevO_; }
            inline bool nextO() { return nextO_; }
            inline bool prevO() { return prevO_; }
        private:
            TFE next_;
            TFE prev_;
            bool nextO_;
            bool prevO_;
    };

    // A pair of TVEs along with the orientation indicating how they are
    // identified.
    class Pair {
        public:
            inline const TVE getA() { return a_;}
            inline const TVE getB() { return b_;}
            inline const bool o() { return orientation_;}
            inline const opp(TVE thing) { return (thing == a_) ? b_ : a_; }
            inline const bool degreeOk() { return degree_ >= 4; }
            inline void addToDegree(int i) { degree_ += i; }
        private:
            TVE a_,b_;
            bool orientation_; // true means "lowest vertex on edge a meets
                    // lowest vertex on edge b"
            int degree_;
    };

    // Union-find structure for tracking equivalent vertices
    class Equiv {
        public:
            Equiv(TV a);
            inline root() {if (parent_ != NULL) return parent_->root(); return this;}
            inline Equiv* parent() { return parent_; }
            inline TV value() { return val_; }
            inline void setParent(Equiv *a) {parent_ = a;}
            inline int height() { return height_;}
            inline void increaseHeight() { height_++; }
        private:
            Equiv* parent_;
            int height_;
            TV val_;
    };

    // The configuration of the boundary of a partial triangulation. This can
    // be constructed either from an actual triangulation, or just from
    // building upon a configuration of a "smaller" configuration. This second
    // fact is the basis of the FPT algorithm for admissibility testing.
    class Config {
        public:
            // Constructor
            Config(int nChildren);
            // Copy constructor
            Config(const Config& c);

            // The three faces about a face, such that FACE_EDGES[f] give the
            // three edges about face "f"
            static const int[4][3] FACE_EDGES;
            // The three vertices about a face, as above.
            static const int[4][3] FACE_VERTICES;
            // The symmetries used in this gluing are numbered as follows.
            // 0 = 012->012
            // 1 = 012->021
            // 2 = 012->102
            // 3 = 012->120
            // 4 = 012->201
            // 5 = 012->210
            // In the above, we are mapping vertices to vertex, where 0,1,2
            // represent the three vertices, in natural ordering (so triangle
            // 230 is of type "120" since 0<2<3).
            // Given face f, and symmetry s, the vertex v is glued to VERT_SYM_MAP[s][v]
            static const int[4][3] VERT_SYM_MAP;

            // Given face f, and symmetry s, the edge e is glued to
            // EDGE_SYM_MAP[s][f], where this is a value in 0,1,2. To get an
            // actual edge, use FACE_EDGES[ EDGE_SYM_MAP ...
            static const int[4][3] EDGE_SYM_MAP;
            // Given the above edge gluing, EDGE_ORIENT_MAP is true or false
            // based on whether the orientations will agree (true) or disagree
            // (false)
            static const bool[4][3] EDGE_ORIENT_MAP;
            // Given a face f, the edge opposite vertex v is OPP_EDGE[f][v]
            static const int[4][4] OPP_EDGE;

        private:

            // The child configs that created this Config. We use this to
            // lookup details of the configuration, so we don't have to keep
            // multiple copies of the same data.
            Config* children;

            int numChildren_;

            // For each TV, we want to easily be able to find the set of equivalent
            // TVs. That is, which TVs are actually part of the same
            // triangulation-vertex.
            std::map<TV, std::set<TV>*> equivMap;

            // Maps two pairs of TVEs together, along with a boolean indicating
            // whether the mapping preserves orientability (true) or not
            // (false).
            std::map<TVE, Pair> pairs;

            // Is this config useful, aka does this config lead to some
            // actually interesting triangulation.
            bool useful;


            // For the following, we will temporarily store new data structures
            // here before adding them to the actual ... data structures.
            // This makes it much easier to "undo" any partial gluings.
            // Actually, it means we won't have to do anything to undo them.
            // However, we will need to make these class variables, so that the
            // various lookup functions can find them.

            // TODO Instead of just arrays, we could make them associative
            // immediately. Would make lookups faster, but inserts slower.

            // We will store any new pairs here temporarily. At the end of
            // glue(), if nothing has gone wrong, we assign them to the Config.
            Pair* newPair[3];
            TVE newTVE[6]; // newTVE[2*i] and newTVE[2*i+1] are the two TVEs in
                           // pair[i]
            int newPairCount;

            // We temporarily store new link edges here, and if nothing has gone wrong
            // by the end of glue() we assign them to the config. Note that for
            // each of the 3 pairs of link edges we merge, we need to update
            // the prev & next of both in the pair, which means we need (up to)
            // 4 new LinkEdges per merging.
            LinkEdge* newLinks[12];
            TFE newTFE[12]; //newTFE[i] is the TFE associated with newLinks[i]
            int newLinkCount;

        public:
            // Combine another configuration with this one. The two configurations
            // must have distinct tetrahedra. Due to the data structures in
            // use, we only need to add another child to merge configs.
            void addChild(Config *other);

            // Undo such a merge, which just involves removing a child
            void removeChild(int child);

            // Is the given face on the given tetrahedra on the boundary?
            bool onBoundary(int t, int f);

            // Glue two faces together. The two faces are given by Arc a,
            // and the exact gluing is given by gluing (an entry into
            // NEdge::some_table TODO
            bool glue(int gluing, Arc& a);

            // Unglue two faces. TODO: Might need more info passed in
            void unGlue(int t1, int f1);

            // Mark a and b as identified pairs. If orientation is true, then the
            // pair are matched such that lowest vertex of a meets lowest vertex of
            // b.
            void addTVEPair(TVE a, TVE b, bool orientation);

            // Given a TVE, return the Pair object that contains it.
            Pair* getTVEPair(TVE a);
    };

    // Represents a triangulation we are building. Note that we need a bit more
    // than an NTriangulation object. Still not set on NTriangulation member or
    // inherited?
    class Triangulation : public NTriangulation {
        public:
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


            // Glue two faces together. The two faces are given by a,
            // and the exact gluing is given by gluing (an entry into
            // NEdge::some_table TODO
            bool glue(int gluing, Arc& a);
    };

    // Represents a bag of the tree decomposition
    class Bag {
        private:
            // The integer contents of this particular bag
            int* contents_;
            // Which tetrahedra are we adding in this bag
            int[] toAdd_;

            // The list of arcs which will be "added" (where adding means
            // identifying related faces) in this bag
            Arc*[] arcs;
            // How many arcs will be added at this bag
            int numArcs;

            // The children of this bag
            std::list<Bag*> children;

            // Parent bag
            Bag * parent;

            // The possible configurations achievable at this bag
            std::list<Config*> configs_;
            // The possible triangulations achievable at this bag
            std::list<Triangulation*> triangulations_;
            // Have all viable configurations been found
            bool configsFound;
            // Have all viable triangulations been found
            bool trisFound;
        public:
            // Constructor
            Bag(NTreeBag *bag);

            // iterator over triangulations found here.
            std::list<Triangulation *>::iterator triangulations();

            // Find all configurations viable for this bag
            void findConfigs();
            // Find all triangulations viable for this bag
            void findTriangulations();


            // Each child has a possibly different number of triangulations.
            // childTriangulations() returns every possible combination of taking
            // one triangulation from each child.
            Config::iterator childConfigs();

            // Does this bag have no configs?
            bool hasNoConfigs();
            // When asked for next Config, start again at start
            bool resetConfigCount();
            // Does this bag have more Configs?
            bool hasNextConfig();
            // Get the next Config
            Config getNextConfig();

            // "Add" arcs to a config. This function "adds" all arcs, and arcs are
            // added to a config by trying each possible gluing for each arc. Any
            // Configs found after adding all arcs are stored at this bag.
            bool addArcs(Config& c);

            // Store a configuration
            void storeConfig(Config& c);

            // Each child has a possibly different number of triangulations.
            // childTriangulations() returns every possible combination of taking
            // one triangulation from each child.
            Triangulation::iterator childTriangulations();

            // Does this bag have no triangulations
            bool hasNoTriangulations();
            // When asked for the next Triangulation, start again at the start
            bool resetTriangulationCount();
            // Does this bag have more triangulations
            bool hasNextTriangulation();

            // Does the given triangulation have a valid boundary configuration
            bool hasValidConfig(const Triangulation& t);

            // Store a triangulation.
            void storeTri(Triangulation *t);
    };


    // The root bag for the tree decomposition
    Bag *root;

};

// Inline functions for TreeDecompSearcher

// Inline functions for TreeDecompSearcher::Arc
inline TreeDecompSearcher::Arc(FacetSpec<3> one_, FacetSpec<3> two_) one(one_), two(two_) { }


// Inline functions for TreeDecompSearcher::Bag

inline void TreeDecompSearcher::Bag::addChild(Bag *c) {
    children.push_back(c);
    c->parent = this;
}

// Inline functions for TreeDecompSearcher::Config

inline TreeDecompSearcher::Config::Config(int numChildren) : useful(false),
    numChildren_(numChildren) {
    children = new Config*[numChildren_];
}

inline TreeDecompSearcher::Config::Config(const Config& other) : useful(false),
    numChildren_(other.coundChildren()) {
    children = new Config*[numChildren_];
    std::copy(other.children, children, numChildren_*sizeof(Config*));
}

inline TreeDecompSearcher::Config::~Config() {
    delete[] children;
}

inline int TreeDecompSearcher::Config::countChildren() {
    return numChildren_;
}

inline void TreeDecompSearcher::Config::addChild(Config *c, int pos) {
    children[pos] = c;
}

inline void TreeDecompSearcher::Config::removeChild(int pos) {
    // NOP. Will get replaced when needed, else it never gets used anyway.
}

// Inline functions for TreeDecompSearcher::LinkEdge

inline TreeDecompSearcher::LinkEdge::LinkEdge(TFE next, bool nextO, TFE prev,
        bool prevO) : { next_(next), nextO_(nextO), prev_(prev), prevO_(prevO)
    { }

} // namespace regina

#endif

