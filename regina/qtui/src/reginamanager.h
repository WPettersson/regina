
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

/*! \file reginamanager.h
 *  \brief The main application class for Regina.
 */

#ifndef __REGINAMANAGER_H_
#define __REGINAMANAGER_H_

#include <QApplication>

class ReginaMain;

/**
 * The main application class for Regina.
 *
 * This class also creates, destroys and manages the individual main windows.
 */

class ReginaManager : public QApplication {
    Q_OBJECT

    public:
        ReginaManager(int &argc, char** argv);
        
        ReginaMain* newWindow(bool firstWindow);

        void aboutToClose(ReginaMain *child);

    protected:
        /**
         * Support QEvent::FileOpen on MacOS.
         */
        bool event(QEvent*);

    private:
        QList<ReginaMain *> children;

};

#endif
