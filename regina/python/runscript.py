#!/bin/false
#
# ----------------------------------------------------------------------
# FOR INTERNAL USE ONLY:  You should not run this program directly.
# Use regina-python instead to run python scripts from the command line.
# ----------------------------------------------------------------------
#
# Regina - A Normal Surface Theory Calculator
# Python Environment Initialisation
#
# Copyright (c) 2003-2011, Ben Burton
# For further details contact Ben Burton (bab@debian.org).
#
# Usage: runscript.py [ import | noimport ]
#                     [ <library> ... ]
#                     [ -- <script> [ <script-arg> ... ]]
#
# Initialises the python environment for working with the Regina
# calculation engine.  Tasks include:
#
#   - Importing the 'regina' module;
#   - If requested, running "from regina import *";
#   - Running a series of provided library scripts;
#   - Running a selected script with the provided command-line arguments,
#       or,
#     printing a startup banner for working in interactive mode.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
# MA 02110-1301, USA.
#

import sys

# --- Import the Regina calculation engine. ---

try:
    import regina
except:
    print 'ERROR: The calculation engine module could not be loaded.'
    sys.exit(1)

# --- Parse the command-line arguments. ---

libNames = []
scriptName = None
scriptArgs = []

if len(sys.argv) < 2:
    print 'ERROR: The import/noimport argument was missing.'
    sys.exit(1)

if sys.argv[1] == 'import':
    from regina import *
elif sys.argv[1] != 'noimport':
    print 'ERROR: The import/noimport argument was incorrect.'
    sys.exit(1)

libsDone = 0
for i in sys.argv[2:]:
    if libsDone:
        if scriptName == None:
            scriptName = i
        else:
            scriptArgs.append(i)
    elif i == '--':
        libsDone = 1
    else:
        libNames.append(i)

# --- Run each library script in turn. ---

for i in libNames:
    try:
        print 'Running ' + i + '...'
        execfile(i)
    except SystemExit:
        pass
    except:
        sys.excepthook(*sys.exc_info())
        print "ERROR: The custom library '" + i + "' could not be executed."
        sys.exit(1)

# --- Run the script or print the banner as required. ---

if scriptName != None:
    sys.argv = [ scriptName ] + scriptArgs
    try:
        execfile(scriptName)
        sys.exit(0)
    except SystemExit:
        pass
    except:
        sys.excepthook(*sys.exc_info())
        print 'ERROR: An error occurred whilst executing the script.'
        sys.exit(1)
else:
    print regina.welcome()

