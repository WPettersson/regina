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
# Copyright (c) 2003-2013, Ben Burton
# For further details contact Ben Burton (bab@debian.org).
#
# Usage: runscript.py [ import | noimport ] [ readline | noreadline ]
#                     [ [-q] <library> ... ]
#                     [ -- <script> [ <script-arg> ... ]]
#
# Initialises the python environment for working with the Regina
# calculation engine.  Tasks include:
#
#   - Importing the 'regina' module;
#   - If requested, running "from regina import *";
#   - If requested, enabling readline tab completion;
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
# As an exception, when this program is distributed through (i) the
# App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or
# (iii) Google Play by Google Inc., then that store may impose any
# digital rights management, device limits and/or redistribution
# restrictions that are required by its terms of service.
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
quiet = False

if len(sys.argv) < 2:
    print 'ERROR: The import/noimport argument was missing.'
    sys.exit(1)
if len(sys.argv) < 3:
    print 'ERROR: The readline/noreadline argument was missing.'
    sys.exit(1)

if sys.argv[1] == 'import':
    from regina import *
elif sys.argv[1] != 'noimport':
    print 'ERROR: The import/noimport argument was incorrect.'
    sys.exit(1)

if sys.argv[2] == 'readline':
    # Enable tab completion through readline, if we can.
    try:
        import rlcompleter, readline
        # readline by default completes an empty string, whereas if 
        # we press tab we want to insert an actual tab character, 
        # so we have our own completion function.

        __internal_python_completer = readline.get_completer()
        def regina_completer(text, state):
          if not text:
            return ('\t', None)[state]
          else:
            return __internal_python_completer(text, state)
        readline.set_completer(regina_completer)

        if 'libedit' in readline.__doc__:
            # Some systems work with libedit, not libreadline, which
            # supports a different set of commands.
            readline.parse_and_bind('bind ^I rl_complete')
        else:
            readline.parse_and_bind('tab: complete')
    except:
        pass
elif sys.argv[2] != 'noreadline':
    print 'ERROR: The readline/noreadline argument was incorrect.'
    sys.exit(1)

libsDone = 0
for i in sys.argv[3:]:
    if i == '-q':
        quiet=True
    elif libsDone:
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
        if not quiet:
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

