
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
class REGINA_API TreeDecompSearcher : public NGluingPermSearcher {

    static const int maxTetrahedra = 20;

    class Config;
    class Triangulation;

    class Bag {
        int[] contents_;
        int[] toAdd_;
        int[] arcs;
        Bag*[] children;
        std::list<Config*> possibleConfigs_;
        int numChildren;
        int numArcs;

        bool configsFound;
        bool trisFound;

        std::list<Bag> children;
        std::list<Config *> childConfigs;
        std::list<Triangulation *> childTriangulations;


        void getChildConfigs();
        bool hasNoConfigs();
        // When asked for next Config, start again at start
        bool resetConfigCount();
        // Does this bag have more Configs?
        bool hasNextConfig();
        // Get the next Config
        Config getNextConfig();

        bool addArcs(Config& c);

        void storeConfig(Config& c);

        void getChildTriangulations();
        bool hasNoTriangulations();
        bool resetTriangulationCount();
        bool hasNextTriangulation();

        bool hasValidConfig(Triangulation& t);
    };

    typedef uint16 TFE;
    TFE inline TFE_(int t, int f, int e) { return (t*4 + f)*6 + e; } ;
    // 4 faces per tet, 6 edges per tet. TODO work with edges per face?

    typedef uint16 TVE;
    TVE inline TVE_(int t, int v, int e) { return (t*4 + v)*6 + e; } ;
    // 4 vertices per tet, 6 edges per tet. TODO work with edges per vertex?

    typedef uint16 TV;
    TV inline TV_(int t, int v) { return t*4 + v; } ;

    class Config {
        std::map<int, std::set<int>*> equivMap;

        Config();
        Config(Config c);

        void mergeWith(Config other);
        void undoMerge(std::set<int> tets);
        bool onBoundary(int t, int f);
        bool glue(int gluing, int t1, int f1, int t2, int f2);
        void unGlue(int t1, int f1);

        addPair(TFE a, TFE b);
    };

    class Triangulation {
        NTriangulation t_;

        Triangulation();
        Triangulation(Triangulation t);

        void mergeWith(Triangulation other);
        void undoMerge(std::set<int> tets);
        bool onBoundary(int t, int f);
        bool glue(int gluing, int t1, int f1, int t2, int f2);
    };



    static inline int combine(int tet, int vert) { return 4*tet+vert; }
};

// Inline functions for TreeDecompSearcher


} // namespace regina

#endif

