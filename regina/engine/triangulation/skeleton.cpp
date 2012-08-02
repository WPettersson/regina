
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include <queue>
#include "maths/nrational.h"
#include "triangulation/ntriangulation.h"

#include <iostream>

namespace regina {

void NTriangulation::calculateSkeleton() const {
    ideal = false;
    valid = true;
    orientable = true;
    standard = true;

#if 0
    checkPermutations();
        // Sets valid to false if gluings are mismatched (which should
        // never happen if the NTetrahedron gluing routines have been
        // used correctly)
#endif

    // Set this now so that any tetrahedron query routines do not try to
    // recursively recompute the skeleton again.
    calculatedSkeleton = true;

    calculateComponents();
        // Sets components, orientable, NComponent.orientable,
        //     NTetrahedron.component
    calculateFaces();
        // Sets faces, NFace.component
    calculateVertices();
        // Sets vertices, NVertex.component, NVertex.linkOrientable.
    calculateEdges();
        // Sets edges, NEdge.component, valid, NEdge.valid
    calculateBoundary();
        // Sets boundaryComponents, NFace.boundaryComponent,
        //     NEdge.boundaryComponent, NVertex.boundaryComponent,
        //     NComponent.boundaryComponents
    calculateVertexLinks();
        // Sets valid, ideal, NVertex.link,
        //     NVertex.linkEulerCharacteristic, NComponent.ideal,
        //     boundaryComponents, NVertex.boundaryComponent
}

void NTriangulation::checkPermutations() const {
    TetrahedronIterator it;

    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++)
        for (int face = 0; face < 4; face++) {
            NTetrahedron * adjacent = (*it) -> adjacentTetrahedron(face);

            if (adjacent) {
                NPerm4 perm = (*it) -> adjacentGluing(face);

                NPerm4 adj_perm = adjacent -> adjacentGluing(perm[face]);

                if (!(perm*adj_perm).isIdentity()) {
                    valid = false;

                    // This printing statement is temporary code 
                    // to be removed once enough people have tested it
                    std::cerr << "ERROR: Permutations of adjacent faces "
                                 "do not match in skeleton.cpp" << std::endl;
                }

                if ((*it) != adjacent -> adjacentTetrahedron(perm[face])) {
                    valid = false;

                    // This printing statement is temporary code 
                    // to be removed once enough people have tested it
                    std::cerr << "ERROR: Adjacency relations do not match"
                                 " in skeleton.cpp" << std::endl;
                }
            }
        }
}

void NTriangulation::calculateComponents() const {
    TetrahedronIterator it;
    NComponent* label;
    NTetrahedron* tet;
    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++)
        (*it)->component = 0;

    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        if (tet->component == 0) {
            label = new NComponent();
            labelComponent(tet, label);
            components.push_back(label);
        }
    }
}

void NTriangulation::labelComponent(NTetrahedron* firstTet,
        NComponent* component) const {
    // Now non-recursive; uses a queue instead.
    // The queue contains tetrahedra from which we need to propagate
    //     component labelling.

    // Use plain C arrays for the queue.  Since each tetrahedron is pushed
    // on at most once, the array size does not need to be very large.

    // Note that we have >= 1 tetrahedron, since firstTet != 0.
    NTetrahedron** queue = new NTetrahedron*[tetrahedra.size()];

    firstTet->component = component;
    component->tetrahedra.push_back(firstTet);
    firstTet->tetOrientation = 1;

    unsigned queueStart = 0, queueEnd = 1;
    queue[0] = firstTet;

    NTetrahedron* tet;
    NTetrahedron* adjTet;
    int face;
    int yourOrientation;
    while (queueStart < queueEnd) {
        tet = queue[queueStart++];

        for (face=0; face<4; face++) {
            adjTet = tet->adjacentTetrahedron(face);
            if (adjTet) {
                yourOrientation = (tet->adjacentGluing(face).
                    sign() == 1 ? -tet->tetOrientation : tet->tetOrientation);
                if (adjTet->component) {
                    if (yourOrientation != adjTet->tetOrientation) {
                        orientable = false;
                        component->orientable = false;
                    }
                } else {
                    adjTet->component = component;
                    component->tetrahedra.push_back(adjTet);
                    adjTet->tetOrientation = yourOrientation;

                    queue[queueEnd++] = adjTet;
                }
            }
        }
    }

    delete[] queue;
}

void NTriangulation::calculateVertices() const {
    TetrahedronIterator it;
    int vertex;
    NTetrahedron* tet;
    NVertex* label;
    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (vertex=0; vertex<4; vertex++)
            tet->vertices[vertex] = 0;
    }

    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (vertex=0; vertex<4; vertex++)
            if (! tet->vertices[vertex]) {
                label = new NVertex(tet->component);
                tet->component->vertices.push_back(label);
                labelVertex(tet, vertex, label);
                vertices.push_back(label);
            }
    }
}

void NTriangulation::labelVertex(NTetrahedron* firstTet, int firstVertex,
        NVertex* label) const {
    // Create a queue using simple arrays.
    // Since each tetrahedron vertex is pushed on at most once, the
    // array size does not need to be very large.

    // Note that we have >= 1 tetrahedron, since firstTet != 0.
    NTetrahedron** queueTet = new NTetrahedron*[tetrahedra.size() * 4];
    int* queueVtx = new int[tetrahedra.size() * 4];

    firstTet->vertices[firstVertex] = label;
    firstTet->vertexMapping[firstVertex] = NPerm4(0, firstVertex);
    firstTet->tmpOrientation[firstVertex] = 1;
    label->embeddings.push_back(NVertexEmbedding(firstTet, firstVertex));

    unsigned queueStart = 0, queueEnd = 1;
    queueTet[0] = firstTet;
    queueVtx[0] = firstVertex;

    NTetrahedron* tet;
    NTetrahedron* adjTet;
    int vertex;
    int adjVertex;
    int adjOrientation;
    int face;
    NPerm4 adjMap;

    while (queueStart < queueEnd) {
        tet = queueTet[queueStart];
        vertex = queueVtx[queueStart];
        queueStart++;

        for (face=0; face<4; face++) {
            if (face == vertex) continue;
            adjTet = tet->adjacentTetrahedron(face);
            if (adjTet) {
                // When we choose an adjacent gluing map, throw in a
                // swap to preserve the "orientation" of the cycle
                // formed by the images of 1, 2 and 3.  Note that this
                // only becomes meaningful if the vertex link is an
                // orientable surface (otherwise there is no consistent
                // way to orient these cycles at all).
                adjMap = tet->adjacentGluing(face) *
                    tet->vertexMapping[vertex] * NPerm4(1, 2);
                adjVertex = adjMap[0];

                // We should actually be inverting NFace::ordering[adjVertex].
                // However, all we care about is the sign of the permutation,
                // so let's save ourselves those extra few CPU cycles.
                if ((NFace::ordering[adjVertex] *
                        tet->adjacentGluing(face) *
                        NFace::ordering[vertex]).sign() > 0)
                    adjOrientation = -(tet->tmpOrientation[vertex]);
                else
                    adjOrientation = tet->tmpOrientation[vertex];

                if (adjTet->vertices[adjVertex]) {
                    if (adjTet->tmpOrientation[adjVertex] != adjOrientation)
                        label->linkOrientable = false;
                } else {
                    adjTet->vertices[adjVertex] = label;
                    adjTet->vertexMapping[adjVertex] = adjMap;
                    adjTet->tmpOrientation[adjVertex] = adjOrientation;
                    label->embeddings.push_back(NVertexEmbedding(adjTet,
                        adjVertex));

                    queueTet[queueEnd] = adjTet;
                    queueVtx[queueEnd] = adjVertex;
                    queueEnd++;
                }
            }
        }
    }

    delete[] queueTet;
    delete[] queueVtx;
}

void NTriangulation::calculateEdges() const {
    TetrahedronIterator it;
    int edge;
    NTetrahedron* tet;
    NEdge* label;
    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (edge=0; edge<6; edge++)
            tet->edges[edge] = 0;
    }

    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (edge=0; edge<6; edge++)
            if (! tet->edges[edge]) {
                label = new NEdge(tet->component);
                tet->component->edges.push_back(label);
                labelEdge(tet, edge, label);
                edges.push_back(label);
            }
    }
}

void NTriangulation::labelEdge(NTetrahedron* firstTet, int firstEdge,
        NEdge* label) const {
    // Since tetrahedron edges are joined together in a loop, the depth-first
    // search is really just a straight line in either direction.
    // We therefore do away with the usual stack/queue and just keep track
    // of the next edge to process in the current direction.

    firstTet->edges[firstEdge] = label;
    firstTet->edgeMapping[firstEdge] = NEdge::ordering[firstEdge];
    label->embeddings.push_back(NEdgeEmbedding(firstTet, firstEdge));

    // The last tetrahedron edge that was successfully processed.
    NTetrahedron* tet;
    int edge;
    NPerm4 tetVertices;

    int exitFace;
    NPerm4 exitPerm;

    // The next tetrahedron edge around from this.
    NTetrahedron* nextTet;
    int nextEdge;
    NPerm4 nextVertices;

    for (int dir = 0; dir < 2; dir++) {
        // Start at the start and walk in one particular direction.
        tet = firstTet;
        edge = firstEdge;
        tetVertices = tet->edgeMapping[edge];

        while (true) {
            // Move through to the next tetrahedron.
            exitFace = tetVertices[dir == 0 ? 2 : 3];
            nextTet = tet->adjacentTetrahedron(exitFace);
            if (! nextTet)
                break;

            exitPerm = tet->adjacentGluing(exitFace);
            nextVertices = exitPerm * tetVertices * NPerm4(2, 3);
            nextEdge = NEdge::edgeNumber[nextVertices[0]][nextVertices[1]];

            if (nextTet->edges[nextEdge]) {
                // We looped right around.
                // Check that we're not labelling the edge in reverse.
                if (nextTet->edgeMapping[nextEdge][0] != nextVertices[0]) {
                    // The edge is being labelled in reverse!
                    label->valid = false;
                    valid = false;
                }
                break;
            }

            // We have a new tetrahedron edge; this needs to be labelled.
            nextTet->edges[nextEdge] = label;
            nextTet->edgeMapping[nextEdge] = nextVertices;

            if (dir == 0)
                label->embeddings.push_back(NEdgeEmbedding(nextTet, nextEdge));
            else
                label->embeddings.push_front(NEdgeEmbedding(nextTet, nextEdge));

            tet = nextTet;
            edge = nextEdge;
            tetVertices = nextVertices;
        }
    }
}

void NTriangulation::calculateFaces() const {
    TetrahedronIterator it;
    int face;
    NTetrahedron* tet;
    NTetrahedron* adjTet;
    NFace* label;
    NPerm4 adjVertices;
    int adjFace;
    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (face=0; face<4; face++)
            tet->faces[face] = 0;
    }

    for (it = tetrahedra.begin(); it != tetrahedra.end(); it++) {
        tet = *it;
        for (face=3; face>=0; face--)
            if (! tet->faces[face]) {
                label = new NFace(tet->component);
                tet->component->faces.push_back(label);
                tet->faces[face] = label;
                tet->faceMapping[face] = NFace::ordering[face];
                label->embeddings[0] = new NFaceEmbedding(tet, face);
                label->nEmbeddings = 1;
                adjTet = tet->adjacentTetrahedron(face);
                if (adjTet) {
                    // Face is not on the boundary.
                    adjFace = tet->adjacentFace(face);
                    adjVertices = (tet->adjacentGluing(face))*
                        tet->faceMapping[face];
                    adjTet->faces[adjFace] = label;
                    adjTet->faceMapping[adjFace] = adjVertices;
                    label->embeddings[1] = new NFaceEmbedding(adjTet, adjFace);
                    label->nEmbeddings = 2;
                }
                faces.push_back(label);
            }
    }
}

void NTriangulation::calculateBoundary() const {
    // Sets boundaryComponents, NFace.boundaryComponent,
    //     NEdge.boundaryComponent, NVertex.boundaryComponent,
    //     NComponent.boundaryComponents
    FaceIterator it;
    NFace* face;
    NBoundaryComponent* label;

    for (it = faces.begin(); it != faces.end(); it++) {
        face = *it;
        if (face->nEmbeddings < 2)
            if (face->boundaryComponent == 0) {
                label = new NBoundaryComponent();
                label->orientable = true;
                labelBoundaryFace(face, label);
                boundaryComponents.push_back(label);
                face->component->boundaryComponents.push_back(label);
            }
    }
}

void NTriangulation::labelBoundaryFace(NFace* firstFace,
        NBoundaryComponent* label) const {
    std::queue<NFace*> faceQueue;
    NFaceEmbedding* emb;

    emb = firstFace->embeddings[0];
    firstFace->boundaryComponent = label;
    label->faces.push_back(firstFace);
    emb->getTetrahedron()->tmpOrientation[emb->getFace()] = 1;
    faceQueue.push(firstFace);

    NTetrahedron* tet;
    NPerm4 tetVertices;
    int tetFace;
    int i,j;
    NVertex* vertex;
    NEdge* edge;

    NFace* face;
    NFace* nextFace;
    int nextFaceNumber;
    NPerm4 nextFacePerm;
    NTetrahedron* nextTet;
    int followFromFace;
    NPerm4 switchPerm;
    int yourOrientation;

    while (! faceQueue.empty()) {
        face = faceQueue.front();
        faceQueue.pop();

        // Run through the edges and vertices on this face.
        emb = face->embeddings[0];
        tet = emb->getTetrahedron();
        tetFace = emb->getFace();
        tetVertices = tet->faceMapping[tetFace];

        // Run through the vertices.
        for (i=0; i<3; i++) {
            vertex = tet->vertices[tetVertices[i]];
            if (vertex->boundaryComponent != label) {
                // A vertex in an invalid triangulation might end up in
                // more than one boundary component.  Push it into all
                // of the relevant boundary components' lists.
                vertex->boundaryComponent = label;
                label->vertices.push_back(vertex);
            }
        }

        // Run through the edges.
        for (i=0; i<3; i++)
            for (j=i+1; j<3; j++) {
                edge = tet->edges[NEdge::edgeNumber[tetVertices[i]]
                    [tetVertices[j]]];
                if (! (edge->boundaryComponent)) {
                    edge->boundaryComponent = label;
                    label->edges.push_back(edge);
                }

                // Label the adjacent boundary face with the same label.
                followFromFace = 6 - tetVertices[i] - tetVertices[j] - tetFace;
                switchPerm = NPerm4(followFromFace, tetFace);
                nextFaceNumber = followFromFace;
                nextFacePerm = NPerm4();
                nextTet = tet;
                while (nextTet->adjacentTetrahedron(nextFaceNumber)) {
                    nextFacePerm = nextTet->adjacentGluing(
                        nextFaceNumber) * nextFacePerm * switchPerm;
                    nextTet = nextTet->adjacentTetrahedron(nextFaceNumber);
                    nextFaceNumber = nextFacePerm[followFromFace];
                }
                nextFace = nextTet->faces[nextFaceNumber];
                // Find the expected orientation of the next face.
                yourOrientation =
                    (nextTet->faceMapping[nextFaceNumber].inverse() *
                    nextFacePerm * switchPerm * tet->faceMapping[tetFace])
                    .sign() == 1 ? -tet->tmpOrientation[tetFace] :
                    tet->tmpOrientation[tetFace];
                if (nextFace->boundaryComponent) {
                    // Check the orientation.
                    if (yourOrientation !=
                            nextTet->tmpOrientation[nextFaceNumber])
                        label->orientable = false;
                } else {
                    // Add this adjacent face to the queue.
                    nextFace->boundaryComponent = label;
                    label->faces.push_back(nextFace);
                    nextTet->tmpOrientation[nextFaceNumber] = yourOrientation;
                    faceQueue.push(nextFace);
                }
            }
    }
}

void NTriangulation::calculateVertexLinks() const {
    // Begin by calculating Euler characteristics.
    // Here we use the formula:  chi = (2 v_int + v_bdry - f) / 2, which
    // is easily proven with a little arithmetic.

    // Note that NVertex::linkEulerCharacteristic is initialised to 0 in
    // the NVertex constructor.

    // Begin by calculating (2 v_int + v_bdry) for each vertex link.
    NEdge* e;
    NVertex* end0;
    NVertex* end1;
    NTetrahedron* tet;
    for (EdgeIterator eit = edges.begin(); eit != edges.end(); eit++) {
        e = *eit;

        // Try to compute e->getVertex(0) and e->getVertex(1), but
        // without calling e->getVertex() which will recursively try to
        // recompute the skeleton.
        const NEdgeEmbedding& emb = e->getEmbeddings().front();
        tet = emb.getTetrahedron();
        end0 = tet->vertices[tet->edgeMapping[emb.getEdge()][0]];
        end1 = tet->vertices[tet->edgeMapping[emb.getEdge()][1]];

        if (e->isBoundary()) {
            // Contribute to v_bdry.
            end0->linkEulerCharacteristic++;
            if (e->valid)
                end1->linkEulerCharacteristic++;
        } else {
            // Contribute to 2 v_int.
            end0->linkEulerCharacteristic += 2;
            if (e->valid)
                end1->linkEulerCharacteristic += 2;
        }
    }

    // Run through each vertex and finalise Euler characteristic, link
    // and more.

    NVertex* vertex;
    for (VertexIterator it = vertices.begin(); it != vertices.end(); it++) {
        vertex = *it;

        // Fix the Euler characteristic (subtract f, divide by two).
        vertex->linkEulerCharacteristic = (vertex->linkEulerCharacteristic
            - static_cast<long>(vertex->getEmbeddings().size())) / 2;

        if (vertex->isBoundary()) {
            // We haven't added ideal vertices to the boundary list yet,
            // so this must be real boundary.
            if (vertex->linkEulerCharacteristic == 1)
                vertex->link = NVertex::DISC;
            else {
                vertex->link = NVertex::NON_STANDARD_BDRY;
                valid = false;
                standard = false;
            }
        } else {
            if (vertex->linkEulerCharacteristic == 2)
                vertex->link = NVertex::SPHERE;
            else {
                if (vertex->linkEulerCharacteristic == 0)
                    vertex->link = (vertex->isLinkOrientable() ?
                        NVertex::TORUS : NVertex::KLEIN_BOTTLE);
                else {
                    vertex->link = NVertex::NON_STANDARD_CUSP;
                    standard = false;
                }

                ideal = true;
                vertex->component->ideal = true;

                NBoundaryComponent* bc = new NBoundaryComponent(vertex);
                bc->orientable = vertex->isLinkOrientable();
                vertex->boundaryComponent = bc;
                boundaryComponents.push_back(bc);
                vertex->component->boundaryComponents.push_back(bc);
            }
        }
    }
}

void NTriangulation::calculateBoundaryProperties() const {
    // Make sure the skeleton has been calculated!
    if (! calculatedSkeleton)
        calculateSkeleton();

    bool localTwoSphereBoundaryComponents = false;
    bool localNegativeIdealBoundaryComponents = false;

    for (BoundaryComponentIterator it = boundaryComponents.begin();
            it != boundaryComponents.end(); it++) {
        if ((*it)->getEulerCharacteristic() == 2)
            localTwoSphereBoundaryComponents = true;
        else if ((*it)->isIdeal() && (*it)->getEulerCharacteristic() < 0)
            localNegativeIdealBoundaryComponents = true;

        // Stop the search if we've found everything we're looking for.
        if (localTwoSphereBoundaryComponents &&
                localNegativeIdealBoundaryComponents)
            break;
    }

    twoSphereBoundaryComponents = localTwoSphereBoundaryComponents;
    negativeIdealBoundaryComponents = localNegativeIdealBoundaryComponents;
}

} // namespace regina
