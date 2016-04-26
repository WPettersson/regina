


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
    getChildConfigs();
    for(auto config : childConfigs) {
        for(auto tet: tetsAddedHere) {
            config.addTet(tet);
        }
        addArcs(config);
    }
    configusFound = true;
}

void TreeDecompSearcher::Bag::findConfigs() {
    getChildTriangulations();
    for(auto tri: childTriangulations) {
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
    if (child[childCount].hasNoConfigs())
        return;
    child[childCount].resetCount();
    Config c = new Config;
    while (true) {
        if (!child[childCount].hasNextConfig()) {
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
                c.undoMerge(child[childCount].contents());
                --childCount;
            }
        }
        // Get next config from child number childCount
        c.mergeWith(child[childCount].getNextConfig());
        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            childConfigs.append(new Config(c));
            c.undoMerge(child[childCount].contents());
        } else {
            ++childCount;
            child[childCount].resetCount();
        }
    }
}


void TreeDecompSearcher::Bag::addArcs(Config c) {
    int arcNo = 0;
    int gluings[numArcs] = {0};
    while (true) {
        if (arcNo == -1)
            return;
        if (arcNo == numArcs) {
            storeThisConfig(new Config(c)); // Not with duplicates
        }
        if ((arcNo == numArcs) || (gluings[arcNo] == 6)) {
            --arcNo;
            c.unGlue(gluings[arcNo]); // What do we need to store to do this
            ++gluings[arcNo];
        }
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        // Find boundary
        int f1=0,f2=0; // Possibly make these arrays?
        while( ! c.onBoundary(t1,f1))
            ++f1;
        while( ! c.onBoundary(t2,f2))
            ++f2;
        // Don't modify config if this fails
        if (c.glue(gluings[arcNo], t1,f1, t2,f2)) {
            ++arcNo;
            gluings[arcNo] = 0;
        } else {
            ++gluings[arcNo];
        }
    }
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
    if (child[childCount].hasNoTriangulations())
        return;
    child[childCount].resetCount();
    Triangulation t = new Triangulation;
    while (true) {
        if (!child[childCount].hasNextConfig()) {
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
                t.undoMerge(child[childCount].triangulationContains());
                --childCount;
            }
        }
        // Get next config from child number childCount
        t.mergeWith(child[childCount].getNextConfig());

        // TODO mergeWith() has to somehow maintain tetrahedron labels/indices
        // Possibly just use a Triangulation*[] array to access tetrahedron i.

        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            childTriangulations.append(new Triangulation(t));
            t.undoMerge(child[childCount].contents());
        } else {
            ++childCount;
            child[childCount].resetCount();
        }
    }
}

void TreeDecompSearcher::Triangulation::mergeWith(Triangulation t) {
    // copy tetrahedra over to this
    // be careful with tetrahedra indexing!
}

void TreeDecompSearcher::Triangulation::undoMerge(std::set<int> tets) {
    // delete any tetrahedra in tets
}

void TreeDecompSearcher::Bag::identifyFaces(Triangulation t) {
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
            t.unGlue(gluings[arcNo]); // What do we need to store to do this
            // Maybe store t1,f1,t2,f2 as an array
            ++gluings[arcNo];
        }
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        // Find boundary
        int f1=0,f2=0;
        // If we make these arrays, we only need to calculate them once
        while( ! c.onBoundary(t1,f1))
            ++f1;
        while( ! c.onBoundary(t2,f2))
            ++f2;
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

bool valid(Triangulation)

bool TreeDecompSearcher::Triangulation::glue(int gluing, int t1, int f1, int
        t2, int f2, bool (valid*)(Triangulation) ) {
    // Glue t1,f1 to t2,f2 according to gluing
    // Should we check for edges in reverse/bad edge links? These will be
    // checked again by hasValidBoundaryConfig
}

bool TreeDecompSearcher::Bag::hasValidBoundaryConfig(Triangulation t) {
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
            for (auto ve: vertex->embeddings()) {
                //NVertexEmbedding ve;
                // add to equiv
            }
            const Dim2Triangulation link = vertex->buildLink();
            for(auto b: link.getBoundaryComponents()) {
                for(int i = 0; i < b->countEdges(); ++i) {
                    Dim2Edge *e = b->edge(i);
                    // go through e->embedding(j=0 ...)
                    // possibly just 1st? if on bdry
                    // get Dim2Triangle and int edge from embedding
                    // convert triangle/edge/vertex to tetrahedron/vertex/edge
                }
            }
        }
        // Add to c
    }
    return haveConfig(c);
}

}; // namespace
