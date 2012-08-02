
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
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

// Put this before any Qt/KDE stuff so Python 2.3 "slots" doesn't clash.
#include "pythoninterpreter.h"

#include "regina-config.h"
#include "file/nglobaldirs.h"
#include "packet/npacket.h"


#include "pythonoutputstream.h"

#include <boost/python.hpp>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>

// Python includes:
#include <compile.h>
#include <eval.h>
#include <sysmodule.h>

// Cater for the fact that Py_CompileString has changed behaviour
// between Python 2.2 and Python 2.3.
#if PY_VERSION_HEX >= 0x02030000
    static PyCompilerFlags pyCompFlags = { PyCF_DONT_IMPLY_DEDENT };
    #define PY_COMPILE_SINGLE(cmd) \
        (Py_CompileStringFlags(cmd, "<console>", Py_single_input, &pyCompFlags))
#else
    #define PY_COMPILE_SINGLE(cmd) \
        (Py_CompileString(cmd, "<console>", Py_single_input))
#endif

/**
 * WARNING: We never call Py_Finalize().
 *
 * It can't be called during the program since multiple
 * initialise/finalise sequences cause problems with boost.python
 * modules.
 *
 * It can't be called during global object destruction since this seems to
 * be too late and a crash results (PyThreadState_Get: no current thread).
 *
 * Ideally it is called at program exit, before the destruction of
 * global objects.  This is a hassle though, and since there should be
 * no external resource usage and since memory leaks are reclaimed by the
 * operating system, we simply choose to ignore the problem.
 */

regina::NMutex PythonInterpreter::globalMutex;
bool PythonInterpreter::pythonInitialised = false;

PythonInterpreter::PythonInterpreter(PythonOutputStream* pyStdOut,
        PythonOutputStream* pyStdErr) {
    regina::NMutex::MutexLock lock(globalMutex);

    // Acquire the global interpreter lock.
    if (pythonInitialised)
        PyEval_AcquireLock();
    else {
        PyEval_InitThreads();
        Py_Initialize();
        pythonInitialised = true;
    }

    // Create the new interpreter.
    state = Py_NewInterpreter();

    // Record the main module.
    mainModule = PyImport_AddModule("__main__"); // Borrowed ref.
    mainNamespace = PyModule_GetDict(mainModule); // Borrowed ref.

    // Redirect stdout and stderr if appropriate.
    if (pyStdOut || pyStdErr) {
        boost::python::class_<PythonOutputStream, boost::noncopyable>
                ("PythonOutputStream", boost::python::no_init)
            .def("write", &PythonOutputStream::write)
            .def("flush", &PythonOutputStream::flush);

        boost::python::reference_existing_object::
            apply<PythonOutputStream*>::type conv;

        if (pyStdOut)
            PySys_SetObject("stdout", conv(pyStdOut));
        if (pyStdErr)
            PySys_SetObject("stderr", conv(pyStdErr));
    }

    // Release the global interpreter lock.
    PyEval_SaveThread();
}

PythonInterpreter::~PythonInterpreter() {
    regina::NMutex::MutexLock lock(globalMutex);

    // Acquire the global interpreter lock.
    PyEval_RestoreThread(state);

    // Destroy the interpreter.
    Py_EndInterpreter(state);

    // Release the global interpreter lock.
    PyEval_ReleaseLock();
}

bool PythonInterpreter::executeLine(const std::string& command) {
    /**
     * Part of the code for this routine is based upon pyinterp.cpp by
     * Anthony Martin <amartin@varlog.com>.  The file pyinterp.cpp is
     * licensed under the wxWindows License, which is essentially LGPL
     * with an exception regarding binary object code.  Note that this
     * exception only applies to the original pyinterp.cpp, and is not
     * offered for the PythonInterpreter source code listed here.
     */

    // Join the new command with any previous code.
    std::string fullCommand = currentCode + command;

    // Is it an empty command?
    if (isEmptyCommand(fullCommand)) {
        currentCode.clear();
        return true;
    }

    // Create a char* version of the command to which we might choose to
    // add one or two additional newlines.
    char* cmdBuffer = new char[fullCommand.length() + 3];
    strcpy(cmdBuffer, fullCommand.c_str());

    // Acquire the global interpreter lock.
    PyEval_RestoreThread(state);

    // Attempt to compile the command with no additional newlines.
    PyObject* code = PY_COMPILE_SINGLE(cmdBuffer);
    if (code) {
        // Run the code!
        PyObject* ans = PyEval_EvalCode((PyCodeObject*)code,
            mainNamespace, mainNamespace);
        if (ans)
            Py_DECREF(ans);
        else {
            PyErr_Print();
            PyErr_Clear();
        }

        // Clean up.
        Py_DECREF(code);
        state = PyEval_SaveThread();

        delete[] cmdBuffer;
        currentCode.clear();
        return true;
    }

    // Clear the error so it doesn't stay with the thread state and
    // cause wacky syntax errors later on.  Such behaviour has been
    // observed on macosx/fink.
    PyErr_Clear();

    // Attempt to compile the command with final newline.
    cmdBuffer[fullCommand.length()] = '\n';
    cmdBuffer[fullCommand.length() + 1] = 0;

    code = PY_COMPILE_SINGLE(cmdBuffer);
    if (code) {
        // We're waiting on more code.
        Py_DECREF(code);
        state = PyEval_SaveThread();

        delete[] cmdBuffer;
        currentCode = currentCode + command + '\n';
        return false;
    }

    // Extract the full error details in case we wish to display them later.
    // Note that PyErr_Fetch() has the side-effect of clearing the error also.
    PyObject *errType, *errValue, *errTrace;
    PyErr_Fetch(&errType, &errValue, &errTrace);

    PyObject* errStr1 = PyObject_Str(errValue);

    // Attempt to compile the command with two final newlines.
    cmdBuffer[fullCommand.length() + 1] = '\n';
    cmdBuffer[fullCommand.length() + 2] = 0;

    code = PY_COMPILE_SINGLE(cmdBuffer);
    if (code) {
        // We're waiting on more code.
        Py_DECREF(code);
        Py_XDECREF(errType);
        Py_XDECREF(errValue);
        Py_XDECREF(errTrace);
        Py_XDECREF(errStr1);
        state = PyEval_SaveThread();

        delete[] cmdBuffer;
        currentCode = currentCode + command + '\n';
        return false;
    }

    PyObject* errStr2 = extractErrMsg(); // Also clears the error.

    // Compare the two compile errors.
    if (errStr1 && errStr2) {
        if (PyObject_Compare(errStr1, errStr2)) {
            // Errors are different.  We must be waiting on more code.
            Py_XDECREF(errType);
            Py_XDECREF(errValue);
            Py_XDECREF(errTrace);
            Py_DECREF(errStr1);
            Py_DECREF(errStr2);
            state = PyEval_SaveThread();

            delete[] cmdBuffer;
            currentCode = currentCode + command + '\n';
            return false;
        } else {
            // Errors are the same.  It must be a real compile error.
            PyErr_Restore(errType, errValue, errTrace);
            PyErr_Print();
            PyErr_Clear();

            Py_DECREF(errStr1);
            Py_DECREF(errStr2);
            state = PyEval_SaveThread();

            delete[] cmdBuffer;
            currentCode.clear();
            return true;
        }
    } else {
        // This shouldn't happen either.
        Py_XDECREF(errType);
        Py_XDECREF(errValue);
        Py_XDECREF(errTrace);
        Py_XDECREF(errStr1);
        Py_XDECREF(errStr2);
        state = PyEval_SaveThread();

        pleaseReport("Compile error details are not available.");

        delete[] cmdBuffer;
        currentCode.clear();
        return true;
    }
}

bool PythonInterpreter::importRegina() {
    PyEval_RestoreThread(state);

    // Adjust the python path.
    PyObject* path = PySys_GetObject("path"); // Borrowed reference.
    if (path) {
        PyObject* regModuleDir = PyString_FromString(
            regina::NGlobalDirs::pythonModule().c_str());
        PyList_Append(path, regModuleDir);
        Py_DECREF(regModuleDir);
    }

    // Import the module.
    PyObject* regModule = PyImport_ImportModule("regina"); // New ref.
    if (regModule) {
        PyDict_SetItemString(mainNamespace, "regina", regModule);
        Py_DECREF(regModule);
    } else {
        PyErr_Print();
        PyErr_Clear();
    }

    state = PyEval_SaveThread();
    return (regModule != 0);
}

bool PythonInterpreter::setVar(const char* name, regina::NPacket* value) {
    PyEval_RestoreThread(state);

    boost::python::reference_existing_object::
        apply<regina::NPacket*>::type conv;
    PyObject* pyValue = conv(value);

    if (pyValue) {
        PyObject* nameStr = PyString_FromString(name); // New ref.
        PyDict_SetItem(mainNamespace, nameStr, conv(value));
        Py_DECREF(nameStr);
    }

    state = PyEval_SaveThread();
    return (pyValue != 0);
}

bool PythonInterpreter::compileScript(const char* code) {
    PyEval_RestoreThread(state);

    PyObject* ans = Py_CompileString(const_cast<char*>(code), "<script>",
        Py_file_input);
    if (ans) {
        Py_DECREF(ans);
        state = PyEval_SaveThread();
        return true;
    } else {
        PyErr_Print();
        PyErr_Clear();
        state = PyEval_SaveThread();
        return false;
    }
}

bool PythonInterpreter::runScript(const char* code) {
    PyEval_RestoreThread(state);

    PyObject* ans = PyRun_String(const_cast<char*>(code), Py_file_input,
        mainNamespace, mainNamespace);
    if (ans) {
        Py_DECREF(ans);
        state = PyEval_SaveThread();
        return true;
    } else {
        PyErr_Print();
        PyErr_Clear();
        state = PyEval_SaveThread();
        return false;
    }
}

bool PythonInterpreter::runScript(const char* filename, const char* shortName) {
    PyEval_RestoreThread(state);

    PyObject* script = PyFile_FromString(const_cast<char*>(filename), "r");
    if (script) {
        PyObject* ans = PyRun_File(PyFile_AsFile(script),
            const_cast<char*>(shortName),
            Py_file_input, mainNamespace, mainNamespace);
        Py_DECREF(script);

        if (ans) {
            Py_DECREF(ans);
            state = PyEval_SaveThread();
            return true;
        } else {
            PyErr_Print();
            state = PyEval_SaveThread();
            return false;
        }
    } else {
        state = PyEval_SaveThread();
        return false;
    }
}

bool PythonInterpreter::isEmptyCommand(const std::string& command) {
    for (std::string::const_iterator it = command.begin();
            it != command.end(); it++) {
        if (isspace(*it))
            continue;
        else if (*it == '#')
            return true;
        else
            return false;
    }
    return true;
}

void PythonInterpreter::pleaseReport(const char* msg) {
    std::cerr << "ERROR: " << msg << std::endl;
    std::cerr << "       Please report this anomaly to "
        << PACKAGE_BUGREPORT << std::endl;
    std::cerr << "       since this should never occur.\n";
}

PyObject* PythonInterpreter::extractErrMsg() {
    PyObject *errType, *errValue, *errTrace, *ans;

    PyErr_Fetch(&errType, &errValue, &errTrace);
    ans = PyObject_Str(errValue);

    Py_XDECREF(errType);
    Py_XDECREF(errValue);
    Py_XDECREF(errTrace);
    return ans;
}

