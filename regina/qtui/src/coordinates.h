
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

/*! \file coordinates.h
 *  \brief Assists working in different normal surface coordinate systems.
 */

#ifndef __COORDINATES_H
#define __COORDINATES_H

#include "maths/nlargeinteger.h"

#include <QString>
#include <QObject>

namespace regina {
    class NNormalSurface;
    class NTriangulation;
}

namespace Coordinates {
    /**
     * Return a human-readable name for the given coordinate system.
     */
    const char* name(int coordSystem, bool capitalise = true);

    /**
     * Does the given coordinate system generate almost normal
     * surfaces when used with NNormalSurfaceList::enumerate()?
     *
     * Only coordinate systems that are used for enumerating surfaces
     * (not just viewing surfaces) are relevant here.
     */
    bool generatesAlmostNormal(int coordSystem);

    /**
     * Return the number of coordinate columns in the given coordinate
     * system.
     */
    unsigned long numColumns(int coordSystem, regina::NTriangulation* tri);

    /**
     * Return a column header for the given coordinate of the given
     * coordinate system.
     *
     * The associated triangulation may be passed so that more precise
     * information can be returned, though this routine will behave
     * well without it.
     */
    QString columnName(int coordSystem, unsigned long whichCoord,
        regina::NTriangulation* tri = 0);

    /**
     * Return a column description for the given coordinate of the given
     * coordinate system.
     *
     * The associated triangulation may be passed so that more precise
     * information can be returned, though this routine will behave
     * well without it.
     */
    QString columnDesc(int coordSystem, unsigned long whichCoord,
        const QObject *context, regina::NTriangulation* tri = 0);

    /**
     * Return a particular coordinate of a normal surface in the given
     * coordinate system.
     */
    regina::NLargeInteger getCoordinate(int coordSystem,
        const regina::NNormalSurface& surface, unsigned long whichCoord);
};

#endif
