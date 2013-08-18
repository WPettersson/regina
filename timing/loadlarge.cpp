
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Form a census of splitting surface signatures of given order          *
 *                                                                        *
 *  Copyright (c) 2013, Ben Burton                                        *
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

#include <cstdlib>
#include <file/nxmlfile.h>
#include <packet/npacket.h>

void usage(const char* program) {
    std::cerr << "Load a (typically large) file into memory." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: " << program << " [-t] <file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -t : additionally traverse the tree" << std::endl;
    ::exit(1);
}

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3)
        usage(argv[0]);

    bool traverse;
    const char* file;
    if (argc == 2) {
        traverse = false;
        file = argv[1];
    } else if (strcmp(argv[1], "-t"))
        usage(argv[0]);
    else {
        traverse = true;
        file = argv[2];
    }

    std::cout << "Loading..." << std::endl;
    regina::NPacket* tree = regina::readXMLFile(file);
    if (! tree) {
        std::cerr << "ERROR: Could not load file: " << file << std::endl;
        return 1;
    }

    if (traverse) {
        std::cout << "Traversing..." << std::endl;
        for (regina::NPacket* p = tree; p; p = p->nextTreePacket())
            ;
    }

    std::cout << "Cleaning up..." << std::endl;
    delete tree;

    std::cout << "Done." << std::endl;
    return 0;
}