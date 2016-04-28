
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


namespace regina {


TreeDecompSearcher::runSearch() {

    // get tree decomp
    // assign Bags to each bag
    for(auto tri: rootBag.triangulations()) {
        use_(tri);
    }
}

TriIterator TreeDecompSearcher::Bag::triangulations() {
    if (! configsFound) {
        findConfigs();
    }
    if (! trisFound) {
        findTriangulations();
    }
    return storedTriIterator();
}

void TreeDecompSearcher::Bag::findConfigs() {
    for(auto config : childConfigs()) {
        for(auto tet: tetsAddedHere) {
            config.addTet(tet);
        }
        addArcs(config);
    }
    configsFound = true;
}

void TreeDecompSearcher::Bag::findTriangulations() {
    for(auto tri: childTriangulations()) {
        for(auto tet: tetsAddedHere) {
            tri.addTet(tet);
        }
        identifyFaces(tri);
    }
    trisFound = true;
}


void TreeDecompSearcher::Bag::getChildConfigs() {
    if (numChildren == 0) {
        childConfigs.append(new Config);
        return;
    }
    int childCount = 0;
    if (children[childCount].hasNoConfigs())
        return;
    children[childCount].resetCount();
    Config c = new Config;
    while (true) {
        if (!children[childCount].hasNextConfig()) {
            // Nothing more on this child. Go back one child if possible, else
            // we are done.
            if (childCount == 0) {
                return;
            } else {
                // Undo the last merge of new data
                // Only involves removing parts of the configuration which
                // relate to the last child joined on. By construction these
                // "last parts" must only contain things related to a
                // tetrahedron contained in this last child.
                c.undoMerge(children[childCount].contents());
                --childCount;
            }
        }
        // Get next config from child number childCount
        c.mergeWith(children[childCount].getNextConfig());
        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            childConfigs.append(new Config(c));
            c.undoMerge(children[childCount].contents());
        } else {
            ++childCount;
            children[childCount].resetCount();
        }
    }
}


void TreeDecompSearcher::Bag::addArcs(Config& c) {
    int arcNo = 0;
    int gluings[numArcs] = {0};
    if (! arcFacesKnown)
        findArcFaces(c);
    while (arcNo >= 0) {
        if (arcNo == numArcs) {
            storeConfig(new Config(c)); // Not with duplicates
        }
        if ((arcNo == numArcs) || (gluings[arcNo] == 6)) {
            --arcNo;
            c.unGlue(arcs[arcNo].t1, arcs[arcNo].f1); // TODO Do we need more info passed in to undo?
            ++gluings[arcNo];
        }
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        int f1 = arcs[arcNo].f1;
        int f2 = arcs[arcNo].f2;
        if (c.glue(gluings[arcNo], t1,f1, t2,f2)) {
            ++arcNo;
            gluings[arcNo] = 0;
        } else {
            ++gluings[arcNo];
        }
    }
}

void TreeDecompSearcher::Bag::findArcFaces(Config& c) {
    bool f[4*nTets_] = {false}; // 4 = num faces per tet
    for(int arcNo = 0; arcNo < numArcs; ++arcNo) {
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        while ( ( !c.onBoundary(t1,f1)) || (used[4*t1+f1] == true))
            ++f1;
        f[4*t1+f1] = true;
        while ( ( !c.onBoundary(t2,f2)) || (used[4*t2+f2] == true))
            ++f2;
        f[4*t2+f2] = true;
        arcs[arcNo].f1 = f1;
        arcs[arcNo].f2 = f2;
    }
    arcFacesKnown = true;
}

// Config has:
// circular list of edges (TetVertexEdge) around each vertex link
// map? from each TVE to its corresponding entry in circular list
// map from TetVertex to the set of equivalent TetVertex es. (pointer, multiple
// TetVertex should point to the same set.)
// list of pairs of edges (TetFaceEdge) on each face, along with degree of
// underlying edge and a bool tracking whether a tetrahedron has been repeated
// the edge from TetFaceEdge is signed to store orientation.

bool TreeDecompSearcher::Config::mergeWith(Config) {
    // copy ownership of circular list, merge maps from TVE to circ. list
    // merge maps from TetVertex to equivalent TetVertex
    // merge list of pairs of edges
}

void TreeDecompSearcher::Config::undoMerge(std::set<int> tets) {
    // delete any circular list nodes relating to tetrahedra in tets
    // delete same from equiv. map, and list of pairs of edges
}

bool TreeDecompSearcher::Config::glue(int gluing, int t1, int f1, int t2, int f2) {

    for(int i=0; i< 3; ++i) {
        // Get entries from circular list for t1,f1,t2,f2 and gluing using
        // edge i of f1. For now, TVE[0] and TVE[1]

        // Work out if they are part of same puncture (walk through circular
        // list)

        // Find out if orientations match up.

        // Check the map of equivalent TetVertex to see if TVE[0] and TVE[1]
        // are part of the same link

        // If same link
        //    If (different puncture || orientations don't match)
        //       fail.

        // Get TFE from list of pairs of edges corresponding to
        // t1,f1,t2,f2,gluing,i

        // If TFE[0] matches with TFE[1] in a pair in the list
        //    If orientations don't match || degree in pair < 3
        //       fail.
        //    If degree in pair == 3 && ! repeated_tetrahedron
        //       fail.

    }

    for(int i=0; i < 3; ++i) {
        // Joining things now.
        // For pairs, either close up (+delete?) pair or join two distinct
        // pairs, increasing degree, tracking repeated_tetrahedron and watching
        // orientability

        // For edges around vertex links, merge the two together but watch for
        // orientation of the two.
        // For equivalence map, merge sets and point them at the right place?
    }
}

void TreeDecompSearcher::Bag::getChildTriangulations() {
    if (numChildren == 0) {
        childTriangulations.append(new Triangulation);
        return;
    }
    int childCount = 0;
    if (children[childCount].hasNoTriangulations())
        return;
    children[childCount].resetTriangulationCount();
    Triangulation t = new Triangulation;
    while (true) {
        if (!children[childCount].hasNextTriangulation()) {
            // Nothing more on this child. Go back one child if possible, else
            // we are done.
            if (childCount == 0) {
                return;
            } else {
                // Undo the last merge of new data
                // Only involves removing parts of the configuration which
                // relate to the last child joined on. By construction these
                // "last parts" must only contain things related to a
                // tetrahedron contained in this last child (or subtree below).
                t.undoMerge(children[childCount].triangulationContains());
                --childCount;
            }
        }
        // Get next config from child number childCount
        t.mergeWith(children[childCount].getNextConfig());

        // TODO mergeWith() has to somehow maintain tetrahedron labels/indices
        // Possibly just use a Triangulation*[] array to access tetrahedron i.

        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            childTriangulations.append(new Triangulation(t));
            t.undoMerge(children[childCount].contents());
        } else {
            ++childCount;
            child[childrenCount].resetCount();
        }
    }
}

void TreeDecompSearcher::Triangulation::mergeWith(Triangulation other) {
    // copy tetrahedra over to this
    // be careful with tetrahedra indexing!
    // Maybe just glue based on other?
}

void TreeDecompSearcher::Triangulation::undoMerge(std::set<int> tets) {
    // delete any tetrahedra in tets
    // Or just unglue any glued faces on the tets?
}

void TreeDecompSearcher::Bag::identifyFaces(Triangulation& t) {
    int arcNo = 0;
    int gluings[numArcs] = {0};
    while (true) {
        if (arcNo == -1)
            return;
        if (arcNo == numArcs) {
            if (hasValidBoundaryConfig(t)
                storeTri(new Triangulation(t)); // Not with duplicates
        }
        if ((arcNo == numArcs) || (gluings[arcNo] == 6)) {
            --arcNo;
            t.unGlue(arcs[arcNo].t1, arcs[arcNo].f1); // What do we need to store to do this
            // Maybe store t1,f1,t2,f2 as an array
            ++gluings[arcNo];
        }
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        int f1 = arcs[arcNo].f1;
        int f2 = arcs[arcNo].f2;
        auto f = std::bind(&hasValidConfig, this, _1);
        // Don't modify triangulation if this fails
        if (t.glue(gluings[arcNo], t1,f1, t2,f2, f)) {
            ++arcNo;
            gluings[arcNo] = 0;
        } else {
            ++gluings[arcNo];
        }
    }
}

bool TreeDecompSearcher::Triangulation::glue(int gluing, int t1, int f1, int
        t2, int f2, bool (valid*)(Triangulation&) ) {
    // Glue t1,f1 to t2,f2 according to gluing
    // Should we check for edges in reverse/bad edge links? These will be
    // checked again by hasValidBoundaryConfig

    // Check for valid boundary configuration
    if (! valid(this)) {
        // Unglue
        return false;
    }
    return true;
}

bool TreeDecompSearcher::Bag::hasValidConfig(Triangulation& t) {
    Config c; // Will put configuration of t in here.
    // For each edge of triangulation
    for(auto edge: triangulation.edges()) {
        // TODO
        // prefer if (! edge->isBoundary()) continue;
        if (edge->isBoundary() {
            NEdgeEmbedding ne = edge->embedding(0);
            int teta = ne.tetrahedron();
            int va1 = ne.vertices()[0];
            int va2 = ne.vertices()[1];
            int ea = NEdge::someTable[va1][va2];
            // check orientation?

            ne = edge->embedding(edge->countEmbeddings()-1);
            int tetb = ne.tetrahedron();
            int vb1 = ne.vertices()[0];
            int vb2 = ne.vertices()[1];
            int eb = NEdge::someTable[vb1][vb2];
            // check orientation?

            int fa, fb;
            if ((teta == tetb) && (va1 == vb1) && (va2 == vb2)) {
                fa = NFace::someTable[va1][va2][0];
                fb = NFace::someTable[vb1][vb2][1];
            } else {
                // Find correct face by checking the two faces the edge separates,
                // and selecting the one that is on the boundary.
                fa = NFace::someTable[va1][va2][0];
                if ( ! (t.tetrahedron(teta)->face(fa)->isBoundary()))
                    fa = NFace::someTable[va1][va2][1];

                fb = NFace::someTable[vb1][vb2][0];
                if ( ! (t.tetrahedron(tetb)->face(fb)->isBoundary()))
                    fb = NFace::someTable[vb1][vb2][1];
            }
        }
        TFE tfeA(teta,fa,ea);
        TFE tfeB(tetb,fb,eb);
        c.addTFEPair(tfeA,tfeB);
    }

    for(auto vertex: triangulation.vertices()) {
        if (vertex->isBoundary()) {
            int embeddings = vertex->countEmbeddings();
            int * triToTet = new int[embeddings];
            int * triToFace = new int[embeddings];
            std::set<int> *equiv = new std::set<int>;
            for (int tri = 0; tri < embeddings; ++tri) {
                NVertexEmbedding& ve = vertex.embedding(tri);

                int tet = ve.tetrahedron().index();
                int face = ve.tetrahedron().triangle(ve.vertex());

                equiv->insert(combine(tet, vertex->index()));
                c.equivMap->insert(std::pair<int, std::set<int>*>(combine(tet,
                                vertex->index(), equiv)));


                // Store these for lookup when constructing link boundaries
                triToTet[tri] = tet;
                triToFace[tri] = face;
            }
            const Dim2Triangulation link = vertex->buildLink();
            for(auto b: link.getBoundaryComponents()) {
                for(int i = 0; i < b->countEdges(); ++i) {
                    Dim2Edge *e = b->edge(i);
                    // Since it's a boundary component, it can only have 1
                    // embedding
                    int triangle = e->embedding(0).triangle().index();
                    int tet = triToTet[triangle];
                    int face = triToFace[triangle];
                    // TODO
                    // get int edge from embedding
                    // convert triangle/edge/vertex to tetrahedron/vertex/edge
                    // create TVE
                    // create circular list node
                    // insert
                }
            }
            delete[] triToTet;
            delete[] triToFace;
        }
        // Add to c
    }
    return haveConfig(c);
}

void TreeDecompSearcher::Config::addTFEPair(TFE a, TFE b) {
    pairs.insert(std::pair<TFE, TFE>>(a,b));
    pairs.insert(std::pair<TFE, TFE>>(b,a));
}

}; // namespace
