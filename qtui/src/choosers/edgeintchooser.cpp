
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
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

// Regina core includes:
#include "triangulation/nedge.h"

// UI includes:
#include "edgeintchooser.h"

#include <algorithm>
#include <QBoxLayout>
#include <QDialogButtonBox>
#include <QLabel>

using regina::NEdge;

EdgeIntChooser::EdgeIntChooser(
        regina::NTriangulation* tri,
        int argMin, int argMax, const QString& argDesc,
        EdgeIntFilterFunc filter, QWidget* parent,
        bool autoUpdate) :
        QComboBox(parent), tri_(tri), filter_(filter),
        argMin_(argMin), argMax_(argMax), argDesc_(argDesc) {
    setMinimumContentsLength(30);
    setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLength);
    if (autoUpdate)
        tri_->listen(this);
    fill();
}

std::pair<NEdge*, int> EdgeIntChooser::selected() {
    if (count() == 0)
        return std::make_pair((NEdge*)(0), (int)(0));
    int curr = currentIndex();
    return (curr < 0 ? std::make_pair((NEdge*)(0), (int)(0)) : options_[curr]);
}

void EdgeIntChooser::select(regina::NEdge* option, int arg) {
    int index = 0;
    std::vector<std::pair<regina::NEdge*, int> >::const_iterator it =
        options_.begin();
    while (it != options_.end()) {
        if (it->first == option && it->second == arg) {
            setCurrentIndex(index);
            return;
        }
        ++it;
        ++index;
    }

    // Not found.
    if (! options_.empty())
        setCurrentIndex(0);
    return;
}

QString EdgeIntChooser::description(regina::NEdge* option, int arg) {
    if (option->getNumberOfEmbeddings() == 1) {
        const regina::NEdgeEmbedding& e0 = option->getEmbedding(0);
        return trUtf8("Edge %1 [%2 %3] — %4 (%5)")
            .arg(tri_->edgeIndex(option))
            .arg(argDesc_)
            .arg(arg)
            .arg(tri_->tetrahedronIndex(e0.getTetrahedron()))
            .arg(e0.getVertices().trunc2().c_str());
    } else {
        const regina::NEdgeEmbedding& e0 = option->getEmbedding(0);
        const regina::NEdgeEmbedding& e1 = option->getEmbedding(1);
        QString base;
        if (option->getNumberOfEmbeddings() == 2)
            base = trUtf8("Edge %1 [%2 %3] — %4 (%5), %6 (%7)");
        else
            base = trUtf8("Edge %1 [%2 %3] — %4 (%5), %6 (%7), ...");
        return base
                .arg(tri_->edgeIndex(option))
                .arg(argDesc_)
                .arg(arg)
                .arg(tri_->tetrahedronIndex(e0.getTetrahedron()))
                .arg(e0.getVertices().trunc2().c_str())
                .arg(tri_->tetrahedronIndex(e1.getTetrahedron()))
                .arg(e1.getVertices().trunc2().c_str());
    }
}

void EdgeIntChooser::fill() {
    regina::NTriangulation::EdgeIterator it;
    int i;
    for (it = tri_->getEdges().begin();
            it != tri_->getEdges().end(); ++it)
        for (i = argMin_; i <= argMax_; ++i)
            if ((! filter_) || (*filter_)(*it, i)) {
                addItem(description(*it, i));
                options_.push_back(std::make_pair(*it, i));
            }
}

EdgeIntDialog::EdgeIntDialog(QWidget* parent,
        regina::NTriangulation* tri,
        int argMin, int argMax, const QString& argDesc,
        EdgeIntFilterFunc filter,
        const QString& title,
        const QString& message,
        const QString& whatsThis) :
        QDialog(parent) {
    setWindowTitle(title);
    setWhatsThis(whatsThis);
    QVBoxLayout* layout = new QVBoxLayout(this);

    QLabel* label = new QLabel(message);
    layout->addWidget(label);

    chooser = new EdgeIntChooser(tri, argMin, argMax, argDesc, filter, this);
    layout->addWidget(chooser);

    QDialogButtonBox* buttonBox = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    layout->addWidget(buttonBox);

    connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
}

std::pair<regina::NEdge*, int> EdgeIntDialog::choose(QWidget* parent,
        regina::NTriangulation* tri,
        int argMin, int argMax, const QString& argDesc,
        EdgeIntFilterFunc filter,
        const QString& title,
        const QString& message,
        const QString& whatsThis) {
    EdgeIntDialog dlg(parent, tri, argMin, argMax, argDesc, filter,
        title, message, whatsThis);
    if (dlg.exec())
        return dlg.chooser->selected();
    else
        return std::make_pair((NEdge*)(0), (int)(0));
}
