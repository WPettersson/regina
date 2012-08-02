
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

/*! \file ntriangulationcreator.h
 *  \brief Allows the creation of 3-manifold triangulations.
 */

#ifndef __NTRIANGULATIONCREATOR_H
#define __NTRIANGULATIONCREATOR_H

#include "../packetcreator.h"

#include <QStackedWidget>

class QCheckBox;
class QComboBox;
class QLineEdit;

/**
 * An interface for creating 3-manifold triangulations.
 */
class NTriangulationCreator : public PacketCreator {
    private:
        /**
         * Internal components
         */
        QWidget* ui;
        QComboBox* type;
        QStackedWidget* details;

        /**
         * Details for specific triangulation types
         */
        QLineEdit* lstParams;
        QLineEdit* lensParams;
        QLineEdit* loopLen;
        QCheckBox* loopTwisted;
        QLineEdit* augParams;
        QLineEdit* sfsParams;
        QLineEdit* isoSig;
        QLineEdit* dehydrationString;
        QLineEdit* splittingSignature;
        QComboBox* exampleWhich;

    public:
        /**
         * Constructor.
         */
        NTriangulationCreator();

        /**
         * PacketCreator overrides.
         */
        QWidget* getInterface();
        regina::NPacket* createPacket(regina::NPacket* parentPacket,
            QWidget* parentWidget);
};

#endif
