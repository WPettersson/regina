
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

#include "maths/nmatrixint.h"
#include "snappea/nsnappeatriangulation.h"
#include "snappea/kernel/kernel_prototypes.h"

namespace regina {

const NAbelianGroup* NSnapPeaTriangulation::homologyFilled() const {
    if (! data_)
        return 0;
    if (h1Filled_.known())
        return h1Filled_.value();

    // Fetch the relation matrix from SnapPea.
    regina::snappea::RelationMatrix sRelns;
    regina::snappea::homology_presentation(data_, &sRelns);
    if (! sRelns.relations)
        return (h1Filled_ = 0);

    // Pass the relations to Regina.
    NMatrixInt rRelns(sRelns.num_rows, sRelns.num_columns);
    unsigned i, j;
    for (i = 0; static_cast<int>(i) < sRelns.num_rows; ++i)
        for (j = 0; static_cast<int>(j) < sRelns.num_columns; ++j)
            rRelns.entry(i, j) = sRelns.relations[i][j];

    regina::snappea::free_relations(&sRelns);

    // Let Regina run Smith normal form.
    NAbelianGroup* ans = new NAbelianGroup();
    ans->addGroup(rRelns);
    return (h1Filled_ = ans);
}

} // namespace regina

