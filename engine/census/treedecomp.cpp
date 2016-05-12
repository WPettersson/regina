
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

// For details on these, see treedecompsearcher.h
int[][] TreeDecompSearcher::FACE_VERTICES = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
int[][] TreeDecompSearcher::FACE_EDGES = {{3,4,5},{1,2,5},{0,2,4},{0,1,3}};
int[][] TreeDecompSearcher::OPP_EDGE = {{-1,6,5,4},{6,-1,3,2},{5,3,-1,1},{4,2,1,-1}};
int[][] TreeDecompSearcher::VERT_SYM_MAP = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
int[][] TreeDecompSearcher::EDGE_SYM_MAP = {{0,1,2},{1,0,2},{0,2,1},{2,0,1},{1,2,0},{2,1,0}};
bool[][] TreeDecompSearcher::EDGE_ORIENT_MAP = {{true,true,true},{true,true,false},{false,true,true},{true,false,false},{false,false,true},{false,false,false}};



TreeDecompSearcher::TreeDecompSearcher(FacePairing pairing) :
    pairing_(pairing) {

    nTets = pairing_->size();
    // get tree decomp
    NTreeDecomposition tree(pairing_, TD_UPPER_GREEDY_FILL_IN);
    tree.makeNice();
    NTreeBag *treeBag = tree.root();
    root = new Bag(treeBag);
    bool* seen = new bool[pairing_->size()] {false};
    root->determineTetrahedra(seen, nTets);
}

void TreeDecompSearcher::Bag::determineTetrahedra(bool* seen, nTets) {
    // Determine which tetrahedra are added to the triangulation at this bag.
    // We do this by a top-down search. The first time we see any tetrahedron
    // in a top-down search is also the last time we would've seen it in a
    // bottom-up search, which ensures correctness.
    toAdd_ = new int[nTets];
    for (int i = 0; i < nTets; ++i)
        toAdd_[i] = NO_TET; // Arbitrary marker for end of array
    int added = 0;
    for(auto element: contents_) {
        if (! seen[element]) {
            seen[element] = true;
            toAdd_[added++] = element;
            // Go through the 4 faces of "element" in pairing_
            // If the the other face has not been seen, then it must appear
            // only below this bag (maybe not obvious that it won't appear in
            // some other "branch", but it's easy to show this).
            for (int f = 0; f < 4 ; ++f) { // 4 = DIM+1
                FacetSpec<3> dest = pairing_->dest(element, f);
                if (seen[dest.simp])
                    continue;
                Arc *a = new Arc(FacetSpec<3>(element,f), dest);
                arcs[numArcs++] = a;
            }
        }
    }
    for(auto child: children) {
        child->determineTetrahedra(seen, nTets);
    }
}

TreeDecompSearcher::Bag::Bag(NTreeBag *treeBag) {
    // At most 4 arcs added per tetrahedron in this bag
    arcs = new Arc*[treeBag->size() * 4];
    contents_ = new int[treeBag->size()];
    for(int i = 0 ; i < treeBag->size(); ++i) {
        contents_[i] = treeBag->element(i);
    }

    NTreeBag *child = treeBag->children();
    while (child != NULL) {
        Bag *c = new Bag(child);
        this->addChild(c);
        child = child->sibling();
    }
}

void TreeDecompSearcher::runSearch() {
    for(auto tri: rootBag.triangulation()) {
        use_(tri, useArgs_);
    }
}

const TriIterator TreeDecompSearcher::Bag::triangulations() {
    if (! trisFound) {
        findTriangulations();
    }
    return triangulations_.const_iterator();
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
    if (! configsFound) {
        findConfigs();
    }
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
    Config c = new Config(numChildren);
    while (true) {
        if (!children[childCount].hasNextConfig()) {
            // Nothing more on this child. Go back one child if possible, else
            // we are done.
            if (childCount == 0) {
                return;
            } else {
                // Undo the last merge of new data
                // Actually, we just remove the last "child"
                c.removeChild(childCount);
                --childCount;
            }
        }
        // Get next config from child number childCount
        c.addChild(children[childCount].getNextConfig(), childCount);
        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            childConfigs.append(new Config(c));
            c.removeChild(childCount);
        } else {
            ++childCount;
            children[childCount].resetCount();
        }
    }
}


void TreeDecompSearcher::Bag::addArcs(Config& c) {
    int arcNo = 0;
    int gluings[numArcs] = {0};
    while (arcNo >= 0) {
        if (arcNo == numArcs) {
            storeConfig(new Config(c)); // Not with duplicates
        }
        if ((arcNo == numArcs) || (gluings[arcNo] == 6)) {
            --arcNo;
            c.unGlue(arcs[arcNo]); // TODO Do we need more info passed in to undo?
            ++gluings[arcNo];
        }
        int t1 = arcs[arcNo].t1;
        int t2 = arcs[arcNo].t2;
        int f1 = arcs[arcNo].f1;
        int f2 = arcs[arcNo].f2;
        if (c.glue(gluings[arcNo], arcs[arcNo])) {
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

bool TreeDecompSearcher::Config::glue(int gluing, Arc& a) {
    // We will store any new pairs here temporarily. At the end of this
    // function, if nothing has gone wrong, we assign them to the Config.
    Pair* newPair[3];
    newPairCount = 0;

    // We temporarily store new link edges here, and if nothing has gone wrong
    // by the end we assign them to the config. Note that for each of the 3
    // pairs of link edges we merge, we need to update the prev & next of both
    // in the pair, which means we need (up to) 4 new LinkEdges per merging.
    LinkEdge* newLinks[12];
    newLinkCount = 0;
    TFE newTFE[12]; //newTFE[i] is the TFE associated with newLinks[i]

    for(int i=0; i < 3; ++i) {
        // Using code blocks/braces to separate the edge and vertex configurations.
        // This first block is for tracking the edge links (and degrees)
        {
            int degreeIncrease = 0;
            int vert = FACE_VERTICES[a.one.facet][i];
            int edge = OPP_EDGE[a.one.facet][vert];
            TVE a = TVE_(a.one.simp, vert, edge);

            vert = FACE_VERTICES[a.two.facet][i];
            edge = OPP_EDGE[a.two.facet][vert];
            TVE b = TVE_(a.two.simp, vert, edge);
            Pair *p = getTVEPair(a);
            Pair *p2 = getTVEPair(b);

            if (p == p2) {
                bool badOrientation = p->o() ^ p2->o() ^
                    EDGE_ORIENT_MAP[a.one.facet][gluing];
                if (badOrientation) {
                    // TODO clear memory
                    for(int j = 0; j < newPairCount; ++j)
                        delete newPair[j];
                    for(int j = 0; j < newLinkCount; ++j)
                        delete newLinks[j];
                    return false;
                }
                int degreeExtra = 0;
                if (a.tet() == b.tet())
                    degreeExtra = 1; // Same tetrahedron

                // Degree of 2 or less is bad. Also, if the degree is 3 and the
                // edge is on 3 distinct tetrahedra we must also fail. To check for
                // 3 distinct tetrahedra, we just add an extra 1 to the degree if
                // we ever spot 2 tetrahedron-edges of the same tetrahedron in the
                // triangulation-edge. Note that we can assume these are adjacent
                // tetrahedron-edges in the triangulation-edge, as otherwise the
                // triangulation-edge must have degree ≥ 4 anyway.
                if ((p->degree() + degreeExtra) < 4) {
                    // TODO clear memory
                    for(int j = 0; j < newPairCount; ++j)
                        delete newPair[j];
                    for(int j = 0; j < newLinkCount; ++j)
                        delete newLinks[j];
                    return false;
                }
                // Else it's ok. This edge is closed, but we don't have to remember
                // that it is closed.
            } else {
                // Two different edges. Create a new Pair object representing the
                // two new open ends, with the right degree value
                TVE newA = p->opp(a);
                TVE newB = p2->opp(b);
                newTVE[2*newPairCount] = newA;
                newTVE[2*newPairCount+1] = newB;
                bool orientation = p->o() ^ p2->o() ^
                    EDGE_SYM_MAP[a.one.facet][gluing];
                if (newA.tet() == newB.tet()) {
                    int deg = p->degree() + p2->degree() + 1;
                    newPair[newPairCount++] = new Pair(newA, newB, orientation, deg);
                } else {
                    int deg = p->degree() + p2->degree();
                    newPair[newPairCount++] = new Pair(newA, newB, orientation, deg);
                }
            }
        }
        // Now for the vertex link tracking
        {
            TFE a = TFE_(a.one.simp, a.one.facet, FACE_EDGES[a.one.facet][i]);
            TFE b = TFE_(a.two.simp, a.two.facet,
                    FACE_EDGES[a.two.facet][EDGE_SYM_MAP[gluing][i]]);
            bool sameOrientation = EDGE_ORIENT_MAP[gluing][i];
            //  True if the two edges we are matching up will have their
            //  orientations aligned. Note that this means that we will have to
            //  "swap" the orientation of one puncture when we walk around it
            //  (draw a picture to see this).


            LinkEdge edgeA = getLinkEdge(a);
            LinkEdge edgeB = getLinkEdge(b);
            // Work out if they are part of same puncture (walk around edgeA,
            // see if we see edgeB.
            bool samePuncture = false;
            TFE nextTFE = getLinkEdge(edgeA.next());
            bool useNext = edgeA.nextO();
            while (nextTFE != a) {
                if (nextTFE  == b) {
                    samePuncture = true;
                    // Check orientations
                    // We either need to be traversing b backwards (useNext is
                    // false) but have a+b with the same orientation, or
                    // traverse b forwards (useNext is true) but have a+b with
                    // different orientations.
                    if (useNext != sameOrientation) {
                        // TODO clear memory
                        for(int j = 0; j < newPairCount; ++j)
                            delete newPair[j];
                        for(int j = 0; j < newLinkCount; ++j)
                            delete newLinks[j];
                    }
                    break;
                }
                LinkEdge le = getLinkEdge(nextTFE);
                nextTFE = le.next(useNext);
                // Yeah, this looks weird. If le.nextO() is true, then we
                // want to keep the current value of useNext, otherwise
                // switch it.
                useNext = (le.nextO(useNext) == useNext);
            }

            if (!samePuncture) {
                // Arc a.one.simp, a.two.simp
                if (areEquiv(a,b)) {
                    // TODO clear memory
                    for(int j = 0; j < newPairCount; ++j)
                        delete newPair[j];
                    for(int j = 0; j < newLinkCount; ++j)
                        delete newLinks[j];
                    return false;
                }
            }

            // Work out what matches with what.
            // First, is either of the links a puncture of size one?
            bool aSolo = (edgeA.next() == a);
            bool bSolo = (edgeB.next() == b);

            // If both have size one, we don't need to do anything.
            if (aSolo && bSolo) {
                // Do nothing
            } else if (aSolo) {
                // If only one has size one, then in the other link, just "skip"
                // the edge being matched.
                LinkEdge nextLEb = getLinkEdge(edgeB.next());
                // Get the next of the next LE after b
                TFE next2b = nextLEb.next(edgeB.nextO());
                bool next2bo = nextLEb.nextO(edgeB.nextO());
                newTFE[newLinkCount] = edgeB.next();
                newLinks[newLinkCount++] = new LinkEdge(edgeB.prev(),
                        edgeB.prevO(), next2b, next2bo);

                LinkEdge prevLEb = getLinkEdge(edgeB.prev());
                // Get prev of prev. Note that we negate nextO
                TFE prev2b = prevLEb.next(! edgeB.prevO());
                bool prev2bo = prevLEb.nextO(! edgeB.prevO());
                newTFE[newLinkCount] = edgeB.prev();
                newLinks[newLinkCount++] = new LinkEdge(prev2b, prev2bo,
                        edgeB.next(), edgeB.nextO());
            } else if (bSolo) {
                LinkEdge nextLEa = getLinkEdge(edgeA.next());
                // Get the next of the next LE after a
                TFE next2a = nextLEa.next(edgeA.nextO());
                bool next2ao = nextLEa.nextO(edgeA.nextO());
                newTFE[newLinkCount] = edgeA.next();
                newLinks[newLinkCount++] = new LinkEdge(edgeA.prev(),
                        edgeA.prevO(), next2a, next2ao);

                LinkEdge prevLEa = getLinkEdge(edgeA.prev());
                // Get prev of prev. Note that we negate nextO
                TFE prev2a = prevLEa.next(! edgeA.prevO());
                bool prev2ao = prevLEa.nextO(! edgeA.prevO());
                newTFE[newLinkCount] = edgeA.prev();
                newLink[newLinkCount++] = new LinkEdge(prev2a, prev2ao,
                        edgeA.next(), edgeA.nextO());
            } else {
                // Otherwise we need to update the neighbours of both link edges
                TFE preva = edgeA.prev();
                bool prevao = edgeA.prevO();
                LinkEdge prevLEa = getLinkEdge(preva);
                // Get prev of prev. Note that we negate nextO
                TFE prev2a = prevLEa.next(! prevao);
                bool prev2ao = prevLEa.nextO(! prevao);

                TFE nexta = edgeA.next();
                bool nextao = edgeA.nextO();
                LinkEdge nextLEa = getLinkEdge(nexta);
                TFE next2a = nextLEa.next(nextao);
                bool next2ao = nextLE.nextO(nextao);

                TFE prevb = edgeB.prev();
                bool prevbo = edgeB.prevO();
                LinkEdge prevLEb = getLinkEdge(prevb);
                // Get prev of prev. Note that we negate nextO
                TFE prev2b = prevLEa.next(! prevbo);
                bool prev2bo = prevLEa.nextO(! prevbo);

                TFE nextb = edgeB.next();
                bool nextbo = edgeB.nextO();
                LinkEdge nextLEb = getLinkEdge(nextb);
                TFE next2b = nextLEb.next(nextbo);
                bool next2bo = next LEb.nextO(nextbo);
                if (sameOrientation) {
                    newTFE[newLinkCount] = edgeA.prev();
                    newLink[newLinkCount++] = new LinkEdge(prev2a, prev2ao,
                            prevb, ! prevbo);

                    newTFE[newLinkCount] = edgeB.prev();
                    newLink[newLinkCount++] = new LinkEdge(prev2b, prev2bo,
                            preva, ! prevao);

                    newTFE[newLinkCount] = edgeA.next();
                    newLink[newLinkCount++] = new LinkEdge(nextb, !nextbo,
                            next2a, next2ao);

                    newTFE[newLinkCount] = edgeB.next();
                    newLink[newLinkCount++] = new LinkEdge(nexta, !nextao,
                            next2b, next2bo);
                } else { // reverse orientation
                    newTFE[newLinkCount] = edgeA.prev();
                    newLink[newLinkCount++] = new LinkEdge(prev2a, prev2ao,
                            nextb, nextbo);

                    newTFE[newLinkCount] = edgeB.next();
                    newLink[newLinkCount++] = new LinkEdge(preva, prevao,
                            next2b, next2bo);

                    newTFE[newLinkCount] = edgeA.next();
                    newLink[newLinkCount++] = new LinkEdge(prevb, prevbo,
                            next2a, next2ao);

                    newTFE[newLinkCount] = edgeB.next();
                    newLink[newLinkCount++] = new LinkEdge(prev2b, prev2bo,
                            nexta, nextao);
                }
            }

            // Link edges done. Lastly, "merge" equivalent vertices.
            // We have two tetrahedra and two faces in a.{one,two}.{simp,facet}
            // and 0 \leq i \leq 2 which allows us to identify distinct
            // vertices of tetrahedra.
            {
                TV a = TV_(arc.one.simp,FACE_VERTICES[arc.one.facet][i]);
                TV b = TV_(arc.two.simp,
                         FACE_VERTICES[arc.two.facet][VERT_SYM_MAP[gluing][i]]);
                int nComponents = mergeTV(a,b);
                if (nComponents == 0) {
                    if (! rootConfig) {
                        // TODO clear memory
                        for(int j = 0; j < newPairCount; ++j)
                            delete newPair[j];
                        for(int j = 0; j < newLinkCount; ++j)
                            delete newLinks[j];
                        return false;
                    }
                }
            }
        }
    }

    for(int i=0; i < 3; ++i) {
        // Joining things now.
        // copy newPair[0 .. newPairCount] into this

    }
}

void TreeDecompSearcher::Bag::getChildTriangulations() {
    if (children.size() == 0) {
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

        // If we are done, store + revert, else move to next child
        if (childCount == (numChildren-1)) {
            if (hasValidConfig(t))
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
        if (! edge->isBoundary())
            continue;
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
        TFE tfeA(teta,fa,ea);
        TFE tfeB(tetb,fb,eb);
        c.addTFEPair(tfeA,tfeB);
    }

    for(auto vertex: triangulation.vertices()) {
        if (vertex->isBoundary()) {
            int embeddings = vertex->countEmbeddings();
            int * triToTet = new int[embeddings];
            int * triToFace = new int[embeddings];
            std::set<TV> *equiv = new std::set<TV>;
            for (int tri = 0; tri < embeddings; ++tri) {
                NVertexEmbedding& ve = vertex.embedding(tri);

                int tet = ve.tetrahedron().index();
                int face = ve.tetrahedron().triangle(ve.vertex());

                equiv->insert(TV_(tet, vertex->index()));
                c.equivMap->insert(std::pair<TV, std::set<TV>*>(TV_(tet,
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

void TreeDecompSearcher::Config::addTVEPair(TVE a, TVE b, bool orientation) {
    Pair *p = new Pair(a,b,orientation);
    pairs.insert(a, p);
    pairs.insert(b, p);
}

TreeDecompSearcher::Pair* TreeDecompSearcher::Config::getPair(TVE a) {
    // If it's not going to be found here, return NULL
    if (! belowHere[a.tet()])
        return NULL;

    // Check if this bag knows about it.
    auto it = pairs.find(a);
    if (it != pairs.end())
        return it.second;

    // Not here, some child must know of it
    for (auto child: children) {
        Pair *p = child->getPair(a);
        if (p != NULL) {
            return p;
        }
    }
    assert(false); // We should never get here.
    return NULL;
}

TreeDecompSearcher::Pair::Pair(TVE a, TVE b, bool orientation, int degree) :
    a_(a), b_(b), orientation_(orientation), degree_(degree) { }

}; // namespace
