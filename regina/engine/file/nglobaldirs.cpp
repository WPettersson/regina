
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

#include "file/nglobaldirs.h"

namespace regina {

std::string NGlobalDirs::home_(REGINA_DATADIR);
std::string NGlobalDirs::pythonModule_(REGINA_PYLIBDIR);

std::string NGlobalDirs::home() {
    return home_;
}

std::string NGlobalDirs::pythonModule() {
    return pythonModule_;
}

std::string NGlobalDirs::pythonLibs() {
    return home_ + "/pylib";
}

std::string NGlobalDirs::examples() {
    return home_ + "/examples";
}

std::string NGlobalDirs::engineDocs() {
    return home_ + "/engine-docs";
}

std::string NGlobalDirs::data() {
    return home_ + "/internal/data";
}

void NGlobalDirs::setDirs(const std::string& homeDir,
        const std::string& pythonModuleDir) {
    home_ = homeDir;
    pythonModule_ = pythonModuleDir;
}

} // namespace regina

