
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

/*! \file utilities/hashset.h
 *  \brief Includes the Standard Template Library hash set header.
 *  
 *  This header is provided as a convenience to work around the differing
 *  installations of the Standard Template Library extension classes that
 *  are available.
 *
 *  \deprecated Everything related to the \e hash_set and \e hash_map classes
 *  is deprecated, since these classes are not part of the C++ standard.
 *  This material is scheduled to be removed from Regina in version 5.0.
 */

#ifndef __HASHSET_H
#ifndef __DOXYGEN
#define __HASHSET_H
#endif

#warning This header is deprecated, and will be removed in Regina 5.0.

#include <regina-config.h>

#ifdef HAVE_EXT_HASH_SET
    #include <ext/hash_set>
#else
    #include <hash_set>
#endif

#endif

