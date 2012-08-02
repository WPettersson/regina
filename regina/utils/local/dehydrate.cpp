
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Dehydrate a collection of triangulations                              *
 *                                                                        *
 *  Copyright (c) 2005-2011, Ben Burton                                   *
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

/**
 * The given topology data files are read in one after another.
 * All triangulations are extracted and dehydrated where possible.
 *
 * Final results are written to standard output.  The output lists one
 * triangulation per line; each line is of the form:
 *
 *   <triangulation label> : <dehydration string>
 *
 * If the triangulation cannot be dehydrated, the dehydration string will
 * be written as "---".
 *
 * For further information on dehydration strings, see <i>A Census of
 * Cusped Hyperbolic 3-Manifolds</i>, Callahan, Hildebrand and Weeks,
 * Mathematics of Computation 68/225, 1999.
 */

#include "file/nxmlfile.h"
#include "packet/ncontainer.h"
#include "triangulation/ntriangulation.h"

#include <cstdlib>
#include <iostream>

using namespace regina;

void usage(const char* progName, const std::string& error = std::string()) {
    if (! error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <file.rga> [ <file.rga> ... ]\n";
    std::cerr << std::endl;
    std::cerr << "Results are written to standard output.\n";
    std::cerr << "Errors and warnings are written to standard error.\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    int i = 1;

    // argv[i] is the first filename.
    if (i >= argc)
        usage(argv[0], "At least one data file must be given.");

    // Find the triangulations and dehydrate them.
    bool found = false;
    std::string dehydration;
    for ( ; i < argc; i++) {
        NPacket* tree = readXMLFile(argv[i]);
        if (tree == 0) {
            std::cerr << "ERROR: Could not read data from " << argv[i] << '.'
                << std::endl;
            continue;
        }

        for (NPacket* p = tree; p; p = p->nextTreePacket())
            if (p->getPacketType() == NTriangulation::packetType) {
                found = true;
                dehydration = static_cast<NTriangulation*>(p)->dehydrate();

                std::cout << p->getPacketLabel() << " : ";
                if (dehydration.empty())
                    std::cout << "---";
                else
                    std::cout << dehydration;
                std::cout << std::endl;
            }
    }

    if (! found)
        std::cerr << "WARNING: No triangulations found in the given "
            "input file(s)." << std::endl;

    return 0;
}

