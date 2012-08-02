
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

/*! \file foreign/snappea.h
 *  \brief Allows reading and writing SnapPea files.
 */

#ifndef __SNAPPEA_H
#ifndef __DOXYGEN
#define __SNAPPEA_H
#endif

#include <string>
#include "regina-core.h"

namespace regina {

class NTriangulation;

/**
 * \addtogroup foreign Foreign File Formats
 * Reading and writing foreign file formats.
 * @{
 */

/**
 * Reads a triangulation from the given SnapPea file.  A newly allocated
 * triangulation will be returned; it is the user's responsibility to
 * deallocate this when it is finished with.
 *
 * The packet label of the new triangulation will be the manifold name
 * read from the second line of the SnapPea file.  The first line of the
 * SnapPea file must simply be ``<tt>% Triangulation</tt>.
 *
 * If the file could not be read or if the data was not in the correct
 * format, 0 will be returned.
 *
 * \pre The first two lines of the SnapPea file each contain at most
 * 1000 characters.
 *
 * \i18n This routine makes no assumptions about the
 * \ref i18n "character encoding" used in the given file \e name, and
 * simply passes it through unchanged to low-level C/C++ file I/O routines.
 * It assumes however that the \e contents of the file are in UTF-8.
 *
 * @param filename the name of the SnapPea file from which to read.
 * @return a new triangulation containing the data read from the SnapPea
 * file, or 0 on error.
 */
REGINA_API NTriangulation* readSnapPea(const char *filename);

/**
 * Writes the given triangulation to the given file in SnapPea format.
 * All information aside from tetrahedron gluings will be flagged as
 * unknown for SnapPea to recalculate.
 * The manifold name written in the file will be derived from the packet
 * label.
 *
 * \pre The given triangulation is not invalid, and does not contain any
 * boundary faces.
 *
 * \i18n This routine makes no assumptions about the
 * \ref i18n "character encoding" used in the given file \e name, and
 * simply passes it through unchanged to low-level C/C++ file I/O routines.
 * The \e contents of the file will be written using UTF-8.
 *
 * @param filename the name of the SnapPea file to which to write.
 * @param tri the triangulation to write to the SnapPea file.
 * @return \c true if the export was successful, or \c false otherwise.
 */
REGINA_API bool writeSnapPea(const char* filename, NTriangulation& tri);

/**
 * Returns a token derived from the given string.
 * All whitespace characters in the given string will be replaced with
 * an underscore.
 *
 * @param str the string on which to base the token.
 * @return the corresponding token.
 */
REGINA_API std::string stringToToken(const char* str);

/**
 * Returns a token derived from the given string.
 * All whitespace characters in the given string will be replaced with
 * an underscore.
 *
 * @param str the string on which to base the token.
 * @return the corresponding token.
 */
REGINA_API std::string stringToToken(const std::string& str);

/*@}*/

} // namespace regina

#endif

