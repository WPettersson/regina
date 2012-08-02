
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Count thin edge links for a set of data files                         *
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
 * Reads all files *.rga in the current directory.
 *
 * For each triangulation in each data file, this utility counts the number
 * of thin edge links that appear as vertex normal surfaces in quad space.
 *
 * Results for each data file are written to a CSV file (using space
 * separators); the output filename is based on the original regina data
 * filename.  The output directory must be passed as an additional
 * command-line argument, and this directory must already exist.
 */

#include <file/nxmlfile.h>
#include <packet/ncontainer.h>
#include <surfaces/nnormalsurfacelist.h>
#include <triangulation/ntriangulation.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <dirent.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "popt.h"

using namespace regina;

/**
 * Global variables.
 */
std::string outputDir;

/**
 * Helper struct that allows us to sort files by size.
 */
struct DataFile {
    std::string filename;
    long size;

    DataFile() : size(0) {
    }

    void init(const char* fileToUse) {
        filename = fileToUse;

        // If we can't stat, just pretend the size is zero.
        struct stat info;
        if (stat(fileToUse, &info) == 0)
            size = info.st_size;
        else
            size = 0;
    }

    const bool operator < (const DataFile& rhs) const {
        return (size > rhs.size);
    }
};

/**
 * Parse command-line arguments.
 */
bool parseCmdLineOptions(int argc, const char* argv[]) {
    // Set up the command-line arguments.
    poptOption opts[] = {
        POPT_AUTOHELP
        { 0, 0, 0, 0, 0, 0, 0 }
    };

    poptContext optCon = poptGetContext(0, argc, argv, opts, 0);
    poptSetOtherOptionHelp(optCon, "<output_dir>");

    // Parse the command-line arguments.
    int rc = poptGetNextOpt(optCon);
    if (rc != -1) {
        fprintf(stderr, "%s: %s\n\n",
            poptBadOption(optCon, POPT_BADOPTION_NOALIAS), poptStrerror(rc));
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return false;
    }

    const char** otherOpts = poptGetArgs(optCon);
    if (otherOpts && otherOpts[0]) {
        outputDir = otherOpts[0];
        if (otherOpts[1]) {
            fprintf(stderr, "Only one output directory may be supplied.\n\n");
            poptPrintHelp(optCon, stderr, 0);
            poptFreeContext(optCon);
            return false;
        }
    } else {
        fprintf(stderr, "No output directory was supplied.\n\n");
        poptPrintHelp(optCon, stderr, 0);
        poptFreeContext(optCon);
        return false;
    }

    // All done!
    poptFreeContext(optCon);
    return true;
}

/**
 * Help scandir() identify regina data files.
 */
int isRga(const struct dirent* entry) {
   int len = strlen(entry->d_name);
   return (len > 4 && strcmp(".rga", entry->d_name + len - 4) == 0);
}

/**
 * Main routine for dealing with a single data file.
 */
bool process(const std::string& filename) {
    NPacket* tree = readFileMagic(filename);
    if (! tree)
        return false;

    std::ofstream out((outputDir + "/" + filename + ".dat").c_str());
    if (! out) {
        delete tree;
        return false;
    }

    NTriangulation* t;
    NNormalSurfaceList* s;
    long n, i, links;
    for (NPacket* p = tree; p; p = p->nextTreePacket())
        if (p->getPacketType() == NTriangulation::packetType) {
            t = static_cast<NTriangulation*>(p);
            s = NNormalSurfaceList::enumerate(t, NNormalSurfaceList::QUAD);

            links = 0;
            n = s->getNumberOfSurfaces();
            for (i = 0; i < n; ++i)
                if (s->getSurface(i)->isThinEdgeLink().first)
                    ++links;
            out << t->getNumberOfTetrahedra() << ' ' << links << " \""
                << t->getPacketLabel() << '"' << std::endl;

            delete s;
        }

    delete tree;
    return true;
}

/**
 * Main routine that loops through all data files.
 */
int main(int argc, char* argv[]) {
    // Extract command-line options.
    if (! parseCmdLineOptions(argc, (const char**)argv))
        return 1;

    // Find the list of data files to process.
    struct dirent** entries = 0;

    int nEntries = scandir(".", &entries, &isRga, &alphasort);
    if (nEntries < 0) {
        std::cerr << "ERROR: Could not read directory listing." << std::endl;
        return 1;
    }

    // Sort the entries in descending order by size.
    DataFile* files = new DataFile[nEntries + 1 /* in case nEntries == 0 */];
    int i;
    for (i = 0; i < nEntries; ++i)
        files[i].init(entries[i]->d_name);
    std::sort(files, files + nEntries);

    // Process the files.
    int retVal = 0;
    for (int currEntry = 0; currEntry < nEntries; ++currEntry)
        if (! process(files[currEntry].filename)) {
            std::cerr << "ERROR: Could not process "
                << files[currEntry].filename << "." << std::endl;
            retVal = 1;
        }

    // We should delete the entries array, but we're exiting anyway...
    delete[] files;
    return retVal;
}

