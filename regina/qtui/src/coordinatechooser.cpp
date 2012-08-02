
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Qt User Interface                                                    *
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

// Regina core includes:
#include "surfaces/nnormalsurfacelist.h"

// UI includes:
#include "coordinatechooser.h"
#include "coordinates.h"
#include "reginaprefset.h"

#include <algorithm>

using regina::NNormalSurfaceList;

void CoordinateChooser::insertSystem(int coordSystem) {
    addItem(tr(Coordinates::name(coordSystem)));
    systems.push_back(coordSystem);
}

void CoordinateChooser::insertAllCreators() {
    insertSystem(NNormalSurfaceList::STANDARD);
    insertSystem(NNormalSurfaceList::AN_STANDARD);
    insertSystem(NNormalSurfaceList::QUAD);
    insertSystem(NNormalSurfaceList::AN_QUAD_OCT);
    if (ReginaPrefSet::global().surfacesSupportOriented) {
        insertSystem(NNormalSurfaceList::ORIENTED);
        insertSystem(NNormalSurfaceList::ORIENTED_QUAD);
    }
}

void CoordinateChooser::insertAllViewers(regina::NNormalSurfaceList* surfaces) {
    if (surfaces->allowsAlmostNormal()) {
        // For lists created with Regina 4.5.1 or earlier, we have
        // already taken out surfaces with multiple octagons.  Make this
        // clear to the user.
        if (surfaces->getFlavour() == NNormalSurfaceList::AN_LEGACY)
            insertSystem(NNormalSurfaceList::AN_LEGACY);
        else {
            insertSystem(NNormalSurfaceList::AN_STANDARD);
            insertSystem(NNormalSurfaceList::AN_QUAD_OCT);
        }
    } else {
        insertSystem(NNormalSurfaceList::STANDARD);
        insertSystem(NNormalSurfaceList::QUAD);

        if (surfaces->allowsOriented()) {
            insertSystem(NNormalSurfaceList::ORIENTED);
            insertSystem(NNormalSurfaceList::ORIENTED_QUAD);
        }
    }

    insertSystem(NNormalSurfaceList::EDGE_WEIGHT);
    insertSystem(NNormalSurfaceList::FACE_ARCS);
}

void CoordinateChooser::setCurrentSystem(int newSystem) {
    std::vector<int>::const_iterator it =
        std::find(systems.begin(), systems.end(), newSystem);

    if (it != systems.end())
        setCurrentIndex(it - systems.begin());
}

