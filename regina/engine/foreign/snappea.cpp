
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

#include <cctype>
#include <cstring>
#include <fstream>
#include <iomanip>

#include "foreign/snappea.h"
#include "triangulation/ntriangulation.h"

namespace regina {

NTriangulation* readSnapPea(const char* filename) {
    // Open the file.
    std::ifstream in(filename);
    if (!in)
        return 0;

    // Check that this is a SnapPea triangulation.
    char name[1001];
    unsigned len;
    in.getline(name, 1000);
    if (in.fail() || in.eof())
        return 0;
    // Allow junk on the same line following the triangulation marker.
    if (strncmp(name, "% Triangulation", 15) &&
            strncmp(name, "% triangulation", 15))
        return 0;

    // Read in the manifold name.
    in.getline(name, 1000);
    if (in.fail() || in.eof())
        return 0;
    if ((len = strlen(name)) > 0 && name[len - 1] == '\r')
        name[len - 1] = 0;

    // Read in junk.
    std::string tempStr;
    double tempDbl;

    in >> tempStr;         // Solution type
    in >> tempDbl;         // Volume
    in >> tempStr;         // Orientability
    in >> tempStr;         // Chern-Simmon
    if (tempStr[3] == 'k')
        in >> tempDbl;     // Chern-Simmon is known

    unsigned i,j,k;

    // Read in cusp details and ignore them.
    unsigned numOrientCusps, numNonOrientCusps;
    in >> numOrientCusps >> numNonOrientCusps;

    for (i=0; i<numOrientCusps+numNonOrientCusps; i++) {
        in >> tempStr;             // Cusp type
        in >> tempDbl >> tempDbl;  // Filling information
    }

    // Create the new tetrahedra.
    NTriangulation* triang = new NTriangulation();
    triang->setPacketLabel(name);

    unsigned numTet;
    in >> numTet;
    NTetrahedron **tet = new NTetrahedron*[numTet];
    for (i=0; i<numTet; i++)
        tet[i] = triang->newTetrahedron();

    int g[4];
    int p[4][4];

    for (i=0; i<numTet; i++) {
        // Test the state of the input stream.
        if (! in.good()) {
            delete triang;
            delete[] tet;
            return 0;
        }

        // Read in adjacent tetrahedra.
        for (j=0; j<4; j++)
            in >> g[j];

        // Read in gluing permutations.
        for (j=0; j<4; j++) {
            in >> tempStr;
            for (k=0; k<4; k++)
                switch( tempStr[k] ) {
                    case '0': p[j][k] = 0; break;
                    case '1': p[j][k] = 1; break;
                    case '2': p[j][k] = 2; break;
                    case '3': p[j][k] = 3; break;
                    default:
                        delete triang;
                        delete[] tet;
                        return 0;
                }
        }

        // Perform the gluings.
        for (j=0; j<4; j++)
            tet[i]->joinTo(j, tet[g[j]],
                NPerm4(p[j][0], p[j][1], p[j][2], p[j][3]));

        // Read in junk.
        for (j=0; j<4; j++)
            in >> tempStr;
        for (j=0; j<64; j++)
            in >> tempStr;
        for (j=0; j<2; j++)
            in >> tempStr;
    }

    // All done!
    delete[] tet;
    return triang;
}

bool writeSnapPea(const char* filename, NTriangulation& tri) {
    // Open the file.
    std::ofstream out(filename);
    if (!out)
        return 0;

    // Write header information.
    out << "% Triangulation\n";
    if (tri.getPacketLabel().length() == 0)
        out << "Regina_Triangulation\n";
    else
        out << stringToToken(tri.getPacketLabel()) << '\n';

    // Write general details.
    out << "not_attempted 0.0\n";
    out << "unknown_orientability\n";
    out << "CS_unknown\n";

    // Write cusps.
    out << "0 0\n";

    // Write tetrahedra.
    out << tri.getNumberOfTetrahedra() << '\n';

    int i, j;
    for (NTriangulation::TetrahedronIterator it = tri.getTetrahedra().begin();
            it != tri.getTetrahedra().end(); it++) {
        // Although our precondition states that there are no boundary
        // faces, we test for this anyway.  If somebody makes a mistake and
        // calls this routine with a bounded triangulation, we don't want
        // to wind up calling tetrahedronIndex(0) and crashing.
        for (i = 0; i < 4; i++)
            if ((*it)->adjacentTetrahedron(i))
                out << "   " << tri.tetrahedronIndex(
                    (*it)->adjacentTetrahedron(i)) << ' ';
            else
                out << "   -1 ";
        out << '\n';
        for (i = 0; i < 4; i++)
            out << ' ' << (*it)->adjacentGluing(i).toString();
        out << '\n';

        // Incident cusps.
        for (i = 0; i < 4; i++)
            out << "  -1 ";
        out << '\n';

        // Meridians and longitudes.
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 16; j++)
                out << "  0";
            out << '\n';
        }

        // Tetrahedron shape.
        out << "0.0 0.0\n";
    }

    return true;
}

std::string stringToToken(const char* str) {
    std::string ans(str);
    for (std::string::iterator it = ans.begin(); it != ans.end(); it++)
        if (isspace(*it))
            *it = '_';
    return ans;
}

std::string stringToToken(const std::string& str) {
    std::string ans(str);
    for (std::string::iterator it = ans.begin(); it != ans.end(); it++)
        if (isspace(*it))
            *it = '_';
    return ans;
}

} // namespace regina
