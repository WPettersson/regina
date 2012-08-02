
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Identify the type (binary/XML) and version of a data file             *
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

#include "file/nfileinfo.h"
#include <cstdlib>

void usage(const char* progName) {
    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <file> ...\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc < 2)
        usage(argv[0]);

    regina::NFileInfo* info;
    for (int i = 1; i < argc; i++) {
        if (argc != 2)
            std::cout << "[ " << argv[i] << " ]\n";

        info = regina::NFileInfo::identify(argv[i]);
        if (info)
            info->writeTextLong(std::cout);
        else
            std::cout << "Unknown file format or file could not be opened.\n";
        delete info;

        if (argc != 2)
            std::cout << std::endl;
    }

    return 0;
}

