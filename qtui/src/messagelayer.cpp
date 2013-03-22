
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

#include "messagelayer.h"
#include "reginasupport.h"

#include <QApplication>
#include <QBoxLayout>
#include <QStyle>

MessageLayer::MessageLayer(const char* iconName, const QString& defaultText) {
    QBoxLayout* layout = new QHBoxLayout(this);

    layout->addStretch(1);

    int iconSize = QApplication::style()->pixelMetric(
        QStyle::PM_MessageBoxIconSize);

    QLabel* icon = new QLabel(this);
    icon->setPixmap(ReginaSupport::themeIcon(iconName).pixmap(
        iconSize, iconSize));
    layout->addWidget(icon, 0);

    layout->addSpacing(iconSize / 2 + 2 /* shrug */);

    text = new QLabel(defaultText, this);
    text->setWordWrap(true);
    text->setTextInteractionFlags(Qt::TextBrowserInteraction);
    text->setOpenExternalLinks(true);
    layout->addWidget(text, 4);

    layout->addStretch(1);
}
