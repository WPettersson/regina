
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

#include "foreign/dehydration.h"
#include "packet/ncontainer.h"
#include "packet/ntext.h"

#include "dehydrationhandler.h"
#include "../packetfilter.h"

#include <QFile>
#include <QMessageBox>

const DehydrationHandler DehydrationHandler::instance;

regina::NPacket* DehydrationHandler::importData(const QString& fileName,
        QWidget* parentWidget) const {
    QString explnSuffix = QObject::tr("<p>The file should be a plain text "
        "file containing one dehydration string per line.  "
        "Dehydration strings are described in detail in "
        "<i>A census of cusped hyperbolic 3-manifolds</i>, "
        "Callahan, Hildebrand and Weeks, published in "
        "<i>Mathematics of Computation</i> <b>68</b>, 1999.</qt>");

    regina::NPacket* ans = regina::readDehydrationList(
        static_cast<const char*>(QFile::encodeName(fileName)));
    if (! ans) {
        QMessageBox::warning(parentWidget, QObject::tr("Error reading file"),
            QObject::tr("<qt>An error occurred whilst attempting to read "
                "from the file %1.").arg(fileName) + explnSuffix);
        return 0;
    }

    regina::NPacket* last = ans->getLastTreeChild();
    if (last == 0) {
        QMessageBox::warning(parentWidget, QObject::tr("No dehydration strings"),
            QObject::tr("<qt>The selected file does not contain any "
                "dehydration strings.") + explnSuffix);
        return 0;
    } else if (last->getPacketType() == regina::NText::packetType) {
        if (last == ans->getFirstTreeChild()) {
            QMessageBox::warning(parentWidget, 
                QObject::tr("Interpretation failed"), QObject::tr("<qt>None "
                "of the dehydration strings found in the selected file could "
                "be interpreted.") + explnSuffix);
            return 0;
        } else {
            QMessageBox::warning(parentWidget, 
                QObject::tr("Interpretation failed"), QObject::tr("<qt>One or "
                "more of the dehydration strings could not be interpreted. "
                "Details of the error(s) can be found in the final text packet "
                "beneath the newly imported tree.") + explnSuffix);
        }
    }

    // All worked out okay.
    ans->setPacketLabel(QObject::tr("Rehydrated Triangulations").
        toAscii().constData());
    return ans;
}

