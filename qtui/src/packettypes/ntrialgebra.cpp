
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
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

// Regina core includes:
#include "algebra/ngrouppresentation.h"
#include "algebra/nmarkedabeliangroup.h"
#include "maths/numbertheory.h"
#include "triangulation/nhomologicaldata.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "columnlayout.h"
#include "gaprunner.h"
#include "groupwidget.h"
#include "ntrialgebra.h"
#include "reginaprefset.h"
#include "reginasupport.h"

#include <QDir>
#include <QFileInfo>
#include <QHeaderView>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPainter>
#include <QPushButton>
#include <QRegExp>
#include <QScrollArea>
#include <QStyle>
#include <QTreeWidgetItem>
#include <QTextDocument>
#include <QValidator>

using regina::NPacket;
using regina::NTriangulation;

namespace {
    /**
     * How large does r have to be before we start warning the user about
     * Turaev-Viro computation time?
     */
    const unsigned long TV_WARN_LARGE_R = 15;

    /**
     * A regular expression for Turaev-Viro parameters.
     */
    QRegExp reTVParams("^[ \\(]*(\\d+)[ ,]+(\\d+)[ \\)]*$");

    /**
     * A list view item for storing a single Turaev-Viro invariant.
     *
     * These list view items are sorted numerically and drawn with a
     * grid.
     */
    class TuraevViroItem : public QTreeWidgetItem {
        private:
            unsigned long r_;
            unsigned long root_;
            double value_;

        public:
            TuraevViroItem(unsigned long r, unsigned long root, double value) :
                    QTreeWidgetItem(), r_(r), root_(root), value_(value) {
                setText(0, QString::number(r_));
                setText(1, QString::number(root_));
                setText(2, QString::number(value_));

                for (int i = 0; i < 3; ++i)
                    setTextAlignment(i, Qt::AlignRight);
            }

            int compare(const TuraevViroItem* other) const {
                if (r_ < other->r_) return -1;
                if (r_ > other->r_) return 1;
                if (root_ < other->root_) return -1;
                if (root_ > other->root_) return 1;
                return 0;
            }
    };
}

NTriAlgebraUI::NTriAlgebraUI(regina::NTriangulation* packet,
        PacketTabbedUI* useParentUI) :
        PacketTabbedViewerTab(useParentUI,
            ReginaPrefSet::global().tabDim3TriAlgebra) {
    addTab(new NTriHomologyFundUI(packet, this),
        tr("&Homology && Fund. Group"));
    addTab(new NTriTuraevViroUI(packet, this), tr("&Turaev-Viro"));
    addTab(new NTriCellularInfoUI(packet, this), tr("&Cellular Info"));
}

NTriHomologyFundUI::NTriHomologyFundUI(regina::NTriangulation* packet,
        PacketTabbedViewerTab* useParentUI) : PacketViewerTab(useParentUI),
        tri(packet) {
    ui = new QWidget();

    ColumnLayout* master = new ColumnLayout(ui);

    // Homology:

    QGridLayout* homologyGrid = new QGridLayout();//, 7, 4, 0, 5);
    homologyGrid->setRowStretch(0, 1);
    homologyGrid->setRowStretch(6, 1);
    homologyGrid->setColumnStretch(0, 1);
    homologyGrid->setColumnStretch(3, 1);

    QString msg;

    // The text for the following labels differs according to whether or
    // not unicode is enabled.  We therefore set the label texts in
    // refreshLabels(), which is called a little further down.

    labelH1 = new QLabel();
    homologyGrid->addWidget(labelH1, 1, 1);
    H1 = new QLabel(ui);
    homologyGrid->addWidget(H1, 1, 2);
    msg = QObject::tr("The first homology group of this triangulation.");
    labelH1->setWhatsThis(msg);
    H1->setWhatsThis(msg);

    labelH1Rel = new QLabel();
    homologyGrid->addWidget(labelH1Rel, 2, 1);
    H1Rel = new QLabel(ui);
    homologyGrid->addWidget(H1Rel, 2, 2);
    msg = QObject::tr("The relative first homology group of this triangulation "
        "with respect to the boundary.");
    labelH1Rel->setWhatsThis(msg);
    H1Rel->setWhatsThis(msg);

    labelH1Bdry = new QLabel();
    homologyGrid->addWidget(labelH1Bdry, 3, 1);
    H1Bdry = new QLabel(ui);
    homologyGrid->addWidget(H1Bdry, 3, 2);
    msg = QObject::tr("The first homology group of the boundary of this "
        "triangulation.");
    labelH1Bdry->setWhatsThis(msg);
    H1Bdry->setWhatsThis(msg);

    labelH2 = new QLabel();
    homologyGrid->addWidget(labelH2, 4, 1);
    H2 = new QLabel(ui);
    homologyGrid->addWidget(H2, 4, 2);
    msg = QObject::tr("The second homology group of this triangulation.");
    labelH2->setWhatsThis(msg);
    H2->setWhatsThis(msg);

    labelH2Z2 = new QLabel();
    homologyGrid->addWidget(labelH2Z2, 5, 1);
    H2Z2 = new QLabel(ui);
    homologyGrid->addWidget(H2Z2, 5, 2);
    msg = QObject::tr("<qt>The second homology group of this triangulation "
        "with coefficients in Z<sub>2</sub>.</qt>");
    labelH2Z2->setWhatsThis(msg);
    H2Z2->setWhatsThis(msg);

    refreshLabels();

    master->addLayout(homologyGrid, tr("Homology"));

    // Fundamental group:

    QBoxLayout* fundLayout = new QVBoxLayout();

    fgMsg = new QLabel();
    fgMsg->setAlignment(Qt::AlignCenter);
    fundLayout->addWidget(fgMsg);
    fgMsg->hide();

    fgGroup = new GroupWidget(true, true);
    fgGroup->setWhatsThis(tr("A full set of generators and relations "
        "for the fundamental group of this triangulation."));
    connect(fgGroup, SIGNAL(simplified()), this, SLOT(fundGroupSimplified()));
    fundLayout->addWidget(fgGroup, 1);

    master->addLayout(fundLayout, tr("Fundamental Group"));

    connect(&ReginaPrefSet::global(), SIGNAL(preferencesChanged()),
        this, SLOT(updatePreferences()));
}

regina::NPacket* NTriHomologyFundUI::getPacket() {
    return tri;
}

QWidget* NTriHomologyFundUI::getInterface() {
    return ui;
}

void NTriHomologyFundUI::refresh() {
    bool unicode = ReginaPrefSet::global().displayUnicode;

    if (unicode)
        H1->setText(tri->homology().utf8().c_str());
    else
        H1->setText(tri->homology().str().c_str());

    if (tri->isValid()) {
        unsigned long coeffZ2 = tri->homologyH2Z2();

        if (unicode) {
            H1Rel->setText(tri->homologyRel().utf8().c_str());
            H1Bdry->setText(tri->homologyBdry().utf8().c_str());
            H2->setText(tri->homologyH2().utf8().c_str());

            if (coeffZ2 == 0)
                H2Z2->setText("0");
            else if (coeffZ2 == 1)
                H2Z2->setText("\u2124\u2082");
            else
                H2Z2->setText(QString::number(coeffZ2) + " \u2124\u2082");
        } else {
            H1Rel->setText(tri->homologyRel().str().c_str());
            H1Bdry->setText(tri->homologyBdry().str().c_str());
            H2->setText(tri->homologyH2().str().c_str());

            if (coeffZ2 == 0)
                H2Z2->setText("0");
            else if (coeffZ2 == 1)
                H2Z2->setText("Z_2");
            else
                H2Z2->setText(QString::number(coeffZ2) + " Z_2");
        }
    } else {
        QString msg(QObject::tr("Invalid Triangulation"));
        H1Rel->setText(msg);
        H1Bdry->setText(msg);
        H2->setText(msg);
        H2Z2->setText(msg);
    }

    if (tri->countComponents() <= 1) {
        fgMsg->hide();
        fgGroup->refresh(&tri->fundamentalGroup());
        fgGroup->show();
    } else {
        fgGroup->hide();
        fgMsg->setText(tr("<qt>Cannot calculate<p>"
            "(disconnected triangulation)</qt>"));
        fgMsg->show();
    }
}

void NTriHomologyFundUI::fundGroupSimplified() {
    regina::NGroupPresentation* simp = fgGroup->takeSimplifiedGroup();
    if (simp)
        tri->simplifiedFundamentalGroup(simp);
}

void NTriHomologyFundUI::refreshLabels() {
    if (ReginaPrefSet::global().displayUnicode) {
        labelH1->setText(QObject::trUtf8("H\u2081(M):"));
        labelH1Rel->setText(QObject::trUtf8("H\u2081(M, \u2202M):"));
        labelH1Bdry->setText(QObject::trUtf8("H\u2081(\u2202M):"));
        labelH2->setText(QObject::trUtf8("H\u2082(M):"));
        labelH2Z2->setText(QObject::trUtf8("H\u2082(M ; \u2124\u2082):"));
    } else {
        labelH1->setText(QObject::tr("H1(M):"));
        labelH1Rel->setText(QObject::tr("H1(M, bdry M):"));
        labelH1Bdry->setText(QObject::tr("H1(bdry M):"));
        labelH2->setText(QObject::tr("H2(M):"));
        labelH2Z2->setText(QObject::tr("H2(M ; Z_2):"));
    }
}

void NTriHomologyFundUI::updatePreferences() {
    // If we've changed the unicode setting, then we may need some redrawing.
    refreshLabels();
    refresh();
}

NTriTuraevViroUI::NTriTuraevViroUI(regina::NTriangulation* packet,
        PacketTabbedViewerTab* useParentUI) : PacketViewerTab(useParentUI),
        tri(packet) {
    ui = new QWidget();
    QBoxLayout* layout = new QVBoxLayout(ui);

    QBoxLayout* paramsArea = new QHBoxLayout();
    layout->addLayout(paramsArea);
    paramsArea->addStretch(1);

    QString expln = tr("<qt>The (r, root) parameters of a Turaev-Viro "
        "invariant to calculate.  These parameters describe the initial data "
        "for the invariant as described in <i>State sum invariants of "
        "3-manifolds and quantum 6j-symbols</i>, Turaev and Viro, "
        "published in <i>Topology</i> <b>31</b>, no. 4, 1992.<p>"
        "In particular, <i>r</i> and <i>root</i> must both be positive "
        "integers with 0&nbsp;&lt;&nbsp;<i>root</i>&nbsp;&lt;&nbsp;2<i>r</i>, "
        "where <i>root</i> describes a 2<i>r</i>-th root of unity.  "
        "Example parameters are <i>5,3</i>.<p>"
        "Note that only small values of <i>r</i> "
        "should be used, since the time required to calculate the invariant "
        "grows exponentially with <i>r</i>.</qt>");
    paramsLabel = new QLabel(tr("Parameters (r, root):"));
    paramsLabel->setWhatsThis(expln);
    paramsArea->addWidget(paramsLabel);

    params = new QLineEdit(ui);
    params->setValidator(new QRegExpValidator(reTVParams, ui));
    params->setWhatsThis(expln);
    connect(params, SIGNAL(returnPressed()), this, SLOT(calculateInvariant()));
    paramsArea->addWidget(params);

    calculate = new QPushButton(ReginaSupport::themeIcon("system-run"),
        tr("Calculate"));
    // calculate->setFlat(true);
    calculate->setToolTip(tr("Calculate the Turaev-Viro invariant with "
        "these parameters"));
    calculate->setWhatsThis(tr("<qt>Calculate the Turaev-Viro invariant "
        "corresponding to the (r, root) parameters in the nearby text "
        "box.  The result will be added to the list below.<p>"
        "<b>Warning:</b> This calculation can be quite slow for large "
        "values of <i>r</i>, since the processing time grows exponentially "
        "with <i>r</i>.</qt>"));
    connect(calculate, SIGNAL(clicked()), this, SLOT(calculateInvariant()));
    paramsArea->addWidget(calculate);

    paramsArea->addStretch(1);

    QBoxLayout* invArea = new QHBoxLayout();
    layout->addLayout(invArea, 1);
    invArea->addStretch(1);

    invariants = new QTreeWidget(ui);
    invariants->setRootIsDecorated(false);
    invariants->setAlternatingRowColors(true);
    invariants->header()->setStretchLastSection(false);
    invariants->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
    invariants->setSelectionMode(QAbstractItemView::NoSelection);
    invariants->setWhatsThis(tr("A list of all Turaev-Viro invariants "
        "that have been calculated so far for this triangulation.  To "
        "calculate a new invariant, enter the (r, root) parameters into the "
        "text box above and press <i>Calculate</i>."));
    invArea->addWidget(invariants, 1);

    invariants->setColumnCount(3);
    QTreeWidgetItem* header = new QTreeWidgetItem();
    header->setText(0, tr("r"));
    header->setText(1, tr("root"));
    header->setText(2, tr("value"));
    for (int i = 0; i < 3; ++i)
        header->setTextAlignment(i, Qt::AlignCenter);
    invariants->setHeaderItem(header);

    invArea->addStretch(1);

    QLabel* warning = new QLabel(tr("<qt><b>Warning:</b> These are "
        "computed using floating point arithmetic, with no guaranteed level "
        "of accuracy.</qt>"));
    warning->setWordWrap(true);
    warning->setAlignment(Qt::AlignCenter);
    layout->addWidget(warning);
}

regina::NPacket* NTriTuraevViroUI::getPacket() {
    return tri;
}

QWidget* NTriTuraevViroUI::getInterface() {
    return ui;
}

void NTriTuraevViroUI::refresh() {
    paramsLabel->setEnabled(true);
    params->setEnabled(true);
    calculate->setEnabled(true);

    invariants->clear();

    // Since TuraevViroSet is a sorted data type,
    // these items will automatically be inserted in the right order.
    const NTriangulation::TuraevViroSet& invs(tri->allCalculatedTuraevViro());
    for (NTriangulation::TuraevViroSet::const_iterator it = invs.begin();
            it != invs.end(); it++)
        invariants->addTopLevelItem(new TuraevViroItem(
            (*it).first.first, (*it).first.second,
            (*it).second.evaluate((*it).first.second).real()));
}

void NTriTuraevViroUI::calculateInvariant() {
    // Make sure the triangulation is not being edited.
    if (! params->isEnabled())
        return;

    // Run sanity checks.
    if (! (tri->isValid() && tri->isClosed() && ! tri->isEmpty())) {
        ReginaSupport::sorry(ui,
            tr("Turaev-Viro invariants are currently "
            "available only for closed, valid, non-empty triangulations."));
        return;
    }

    if (! reTVParams.exactMatch(params->text())) {
        ReginaSupport::info(ui,
            tr("<qt>The invariant parameters "
            "(<i>r</i>, <i>root</i>) must be two positive integers.</qt>"),
            tr("<qt>These parameters describe the initial data "
            "for the invariant as described in <i>State sum invariants of "
            "3-manifolds and quantum 6j-symbols</i>, Turaev and Viro, "
            "published in <i>Topology</i> <b>31</b>, no. 4, 1992.<p>"
            "In particular, <i>r</i> and <i>root</i> must both be positive "
            "integers with "
            "0&nbsp;&lt;&nbsp;<i>root</i>&nbsp;&lt;&nbsp;2<i>r</i>, "
            "where <i>root</i> describes a 2<i>r</i>-th root of unity.  "
            "Example parameters are <i>5,3</i>.<p>"
            "Note that only small values of <i>r</i> "
            "should be used, since the time required to calculate the "
            "invariant grows exponentially with <i>r</i>.</qt>"));
        return;
    }

    unsigned long r = reTVParams.cap(1).toULong();
    unsigned long root = reTVParams.cap(2).toULong();

    if (r < 3) {
        ReginaSupport::info(ui,
            tr("<qt>The first parameter <i>r</i> must be "
            "at least 3.</qt>"));
        return;
    }

    if (root <= 0 || root >= 2 * r) {
        ReginaSupport::info(ui,
            tr("<qt>The second parameter <i>root</i> "
            "must be strictly between 0 and 2<i>r</i> (it specifies a "
            "2<i>r</i>-th root of unity).</qt>"),
            tr("<qt>Example parameters are <i>5,3</i>.</qt>"));
        return;
    }

    if (regina::gcd(r, root) > 1) {
        ReginaSupport::info(ui,
            tr("<qt>The invariant parameters must have "
            "no common factors.</qt>"),
            tr("<qt>Example parameters are <i>5,3</i>.</qt>"));
        return;
    }

    if (r >= TV_WARN_LARGE_R) {
        QMessageBox msg(QMessageBox::Warning,
            tr("Warning"),
            tr("This calculation could take a very long time."),
            QMessageBox::Yes | QMessageBox::Cancel,
            ui);
        msg.setInformativeText(tr("<qt>The time required "
                "to calculate Turaev-Viro invariants grows exponentially "
                "with <i>r</i>.  It is recommended only to use "
                "r&nbsp;&lt;&nbsp;%1.<p>"
                "Are you sure you wish to "
                "proceed?</qt>").arg(TV_WARN_LARGE_R)),
        msg.setDefaultButton(QMessageBox::Yes);
        if (msg.exec() != QMessageBox::Yes)
            return;
    }

    // Calculate the invariant!
    double value = tri->turaevViro(r, root).evaluate(root).real();
    TuraevViroItem* item = new TuraevViroItem(r, root, value);

    // Insert the invariant in the right place in the table, and delete
    // any previous invariant with the same parameters if it exists.
    TuraevViroItem* curr;
    int cmp;
    for (int i = 0; i < invariants->invisibleRootItem()->childCount(); i++) {
        curr = dynamic_cast<TuraevViroItem*>(
            invariants->invisibleRootItem()->child(i));
        cmp = item->compare(curr);
        if (cmp <= 0) {
            if (cmp == 0)
                delete curr;
            invariants->insertTopLevelItem(i, item);
            return;
        }
    }
    invariants->addTopLevelItem(item);
}

//////////////////////////////////////////////////////////////////////////////
/* RBADD */
/** These routines puts up the interface for the detailed cellular information
        and it is a submenu of the Algebra menu. **/

void NTriCellularInfoUI::refresh() {
    if (tri->isValid()) {
        bool unicode = ReginaPrefSet::global().displayUnicode;

        regina::NHomologicalData minfo(*tri);

        Cells->setText(QObject::tr("%1, %2, %3, %4").
            arg(minfo.countStandardCells(0)).
            arg(minfo.countStandardCells(1)).
            arg(minfo.countStandardCells(2)).
            arg(minfo.countStandardCells(3)));

        DualCells->setText(QObject::tr("%1, %2, %3, %4").
            arg(minfo.countDualCells(0)).
            arg(minfo.countDualCells(1)).
            arg(minfo.countDualCells(2)).
            arg(minfo.countDualCells(3)));

        EulerChar->setText(QString::number(minfo.eulerChar()));

        if (unicode) {
            H0H1H2H3->setText(QObject::trUtf8("H\u2080 = %1,  H\u2081 = %2,  "
                    "H\u2082 = %3,  H\u2083 = %4").
                arg(minfo.homology(0).utf8().c_str()).
                arg(minfo.homology(1).utf8().c_str()).
                arg(minfo.homology(2).utf8().c_str()).
                arg(minfo.homology(3).utf8().c_str()));

            HBdry->setText(
                QObject::tr("H\u2080 = %1,  H\u2081 = %2,  H\u2082 = %3").
                arg(minfo.bdryHomology(0).utf8().c_str()).
                arg(minfo.bdryHomology(1).utf8().c_str()).
                arg(minfo.bdryHomology(2).utf8().c_str()));
        } else {
            H0H1H2H3->setText(
                QObject::tr("H0 = %1,  H1 = %2,  H2 = %3,  H3 = %4").
                arg(minfo.homology(0).str().c_str()).
                arg(minfo.homology(1).str().c_str()).
                arg(minfo.homology(2).str().c_str()).
                arg(minfo.homology(3).str().c_str()));

            HBdry->setText(QObject::tr("H0 = %1,  H1 = %2,  H2 = %3").
                arg(minfo.bdryHomology(0).str().c_str()).
                arg(minfo.bdryHomology(1).str().c_str()).
                arg(minfo.bdryHomology(2).str().c_str()));
        }

        BdryMap->setText(minfo.bdryHomologyMap(1).str().c_str());

        if (! tri->isConnected()) {
            QString msg(QObject::tr("Triangulation is disconnected."));

            TorForOrders->setText(msg);
            TorForSigma->setText(msg);
            TorForLegendre->setText(msg);
            EmbeddingComments->setText(msg);
        } else {
            // 8 principle cases:
            // orientable y/n, boundary y/n, torsion exists y/n
            if (tri->isOrientable()) {
                TorForOrders->setText(
                    minfo.torsionRankVectorString().c_str());
                TorForSigma->setText(
                    minfo.torsionSigmaVectorString().c_str());
                TorForLegendre->setText(
                    minfo.torsionLegendreSymbolVectorString().c_str());
            } else {
                // The torsion linking form routines insist on orientability,
                // so we should avoid calling them.
                QString msg(QObject::tr("Manifold is non-orientable."));

                TorForOrders->setText(msg);
                TorForSigma->setText(msg);
                TorForLegendre->setText(msg);
            }

            // The embeddability comment is good for both orientable and
            // non-orientable triangulations.
            // Encase it in <qt>..</qt> so it can wrap over multiple lines.
            EmbeddingComments->setText(QString("<qt>%1</qt>").arg(
                QString(minfo.embeddabilityComment().c_str()).toHtmlEscaped()));
        }
    } else {
        QString msg(QObject::tr("Invalid Triangulation"));
        Cells->setText(msg);
        DualCells->setText(msg);
        EulerChar->setText(msg);
        H0H1H2H3->setText(msg);
        HBdry->setText(msg);
        BdryMap->setText(msg);
        TorForOrders->setText(msg);
        TorForSigma->setText(msg);
        TorForLegendre->setText(msg);
        EmbeddingComments->setText(msg);
    }
}

NTriCellularInfoUI::NTriCellularInfoUI(regina::NTriangulation* packet,
        PacketTabbedViewerTab* useParentUI) : PacketViewerTab(useParentUI),
        tri(packet) {
    QScrollArea* scroller = new QScrollArea();
    scroller->setWidgetResizable(true);
    scroller->setFrameStyle(QFrame::NoFrame);
    // Transparency must be applied to both the QScrollArea *and* some of its
    // internal components (possibly the widget that holds the viewport?).
    scroller->setStyleSheet("QScrollArea, .QWidget { "
                                "background-color:transparent; "
                            "}");
    ui = scroller;

    QWidget* grid = new QWidget(scroller->viewport());
    scroller->setWidget(grid);

    QGridLayout* homologyGrid = new QGridLayout(grid);//, 11, 4, 0, 5);
    homologyGrid->setRowStretch(0, 1);
    homologyGrid->setRowStretch(11, 1);
    homologyGrid->setColumnStretch(0, 1);
    homologyGrid->setColumnStretch(2, 1); // Give the embeddability comment
                                       // a little room to breathe.
    homologyGrid->setColumnStretch(3, 1);

    QLabel* label;
    QString msg;

    label = new QLabel(QObject::tr("Cells: "), grid);
    homologyGrid->addWidget(label, 1, 1);
    Cells = new QLabel(grid);
    homologyGrid->addWidget(Cells, 1, 2);
    msg = QObject::tr("The number of cells in a proper CW-decomposition of "
               "the compact manifold specified by this triangulation.  "
               "The four numbers displayed here count 0-cells, 1-cells, "
               "2-cells and 3-cells respectively.");
    label->setWhatsThis(msg);
    Cells->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Dual cells: "), grid);
    homologyGrid->addWidget(label, 2, 1);
    DualCells = new QLabel(grid);
    homologyGrid->addWidget(DualCells, 2, 2);
    msg = QObject::tr("The number of cells in the dual CW-decomposition "
                "corresponding to the triangulation of this "
                "compact manifold.  The four numbers displayed here "
                "count 0-cells, 1-cells, 2-cells and 3-cells respectively.");
    label->setWhatsThis(msg);
    DualCells->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Euler characteristic: "), grid);
    homologyGrid->addWidget(label, 3, 1);
    EulerChar = new QLabel(grid);
    homologyGrid->addWidget(EulerChar, 3, 2);
    msg = QObject::tr("The Euler characteristic of this compact manifold.");
    label->setWhatsThis(msg);
    EulerChar->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Homology groups: "), grid);
    homologyGrid->addWidget(label, 4, 1);
    H0H1H2H3 = new QLabel(grid);
    homologyGrid->addWidget(H0H1H2H3, 4, 2);
    msg = QObject::tr("The homology groups of this manifold with coefficients "
               "in the integers.  The groups are listed in order of "
                "increasing dimension.");
    label->setWhatsThis(msg);
    H0H1H2H3->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Boundary homology groups: "), grid);
    homologyGrid->addWidget(label, 5, 1);
    HBdry = new QLabel(grid);
    homologyGrid->addWidget(HBdry, 5, 2);
    msg = QObject::tr("The homology groups of this manifold's boundary with "
               "coefficients in the integers.  The groups are listed "
               "in order of increasing dimension.");
    label->setWhatsThis(msg);
    HBdry->setWhatsThis(msg);

    // The text for the next label differs according to whether or
    // not unicode is enabled.  We therefore set the label text in
    // refreshLabels(), which is called a little further down.

    labelBdryMap = new QLabel(grid);
    homologyGrid->addWidget(labelBdryMap, 6, 1);
    BdryMap = new QLabel(grid);
    homologyGrid->addWidget(BdryMap, 6, 2);
    msg = QObject::tr("<qt>The boundary is a submanifold of the original "
                "manifold.  This item describes some properties of "
                "the induced map on H<sub>1</sub>.</qt>"
                );
    labelBdryMap->setWhatsThis(msg);
    BdryMap->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Torsion form rank vector: "), grid);
    homologyGrid->addWidget(label, 7, 1);
    TorForOrders = new QLabel(grid);
    homologyGrid->addWidget(TorForOrders, 7, 2);
    msg = QObject::tr("<qt>This is the first of the three Kawauchi-Kojima "
               "invariants.  These are invariants of the torsion linking "
               "form on the torsion subgroup of H<sub>1</sub> of an "
               "oriented manifold, and they form a complete set of "
               "invariants.<p>"
               "The torsion form rank vector lists the rank of all the "
               "subgroups of various prime power orders."
               "<p>For further information, see Kawauchi and Kojima's paper "
               "<i>Algebraic classification of linking pairings "
               "on 3-manifolds</i>, Math. Ann. <b>253</b> (1980), "
               "no. 1, 29&ndash;42.</qt>"
                );
    label->setWhatsThis(msg);
    TorForOrders->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Sigma vector: "), grid);
    homologyGrid->addWidget(label, 8, 1);
    TorForSigma = new QLabel(grid);
    homologyGrid->addWidget(TorForSigma, 8, 2);
    msg = QObject::tr("<qt>If H<sub>1</sub> has 2-torsion, this is the "
               "2-torsion sigma-vector.  The sigma-vector is the "
               "second of the Kawauchi-Kojima invariants, and is an "
               "orientation-sensitive invariant."
               "<p>For further information, see Kawauchi and Kojima's paper "
               "<i>Algebraic classification of linking pairings "
               "on 3-manifolds</i>, Math. Ann. <b>253</b> (1980), "
               "no. 1, 29&ndash;42.</qt>"
                );
    label->setWhatsThis(msg);
    TorForSigma->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Legendre symbol vector: "), grid);
    homologyGrid->addWidget(label, 9, 1);
    TorForLegendre = new QLabel(grid);
    homologyGrid->addWidget(TorForLegendre, 9, 2);
    msg = QObject::tr("<qt>If H<sub>1</sub> has odd torsion, this is the "
               "Legendre symbol vector.  The Legendre symbol vector "
               "is the last of the Kawauchi-Kojima invariants, and was "
               "originally constructed by Seifert."
               "<p>For further information, see Kawauchi and Kojima's paper "
               "<i>Algebraic classification of linking pairings "
               "on 3-manifolds</i>, Math. Ann. <b>253</b> (1980), "
               "no. 1, 29&ndash;42.</qt>"
                );
    label->setWhatsThis(msg);
    TorForLegendre->setWhatsThis(msg);

    label = new QLabel(QObject::tr("Comments: "), grid);
    homologyGrid->addWidget(label, 10, 1);
    EmbeddingComments = new QLabel(grid);
    homologyGrid->addWidget(EmbeddingComments, 10, 2);
    msg = QObject::tr("<qt>If the homology allows us to make any deductions "
                "about the embeddability of this manifold in "
                "R<sup>3</sup>, S<sup>3</sup>, S<sup>4</sup> "
                "or a homology sphere, we mention it here.  "
                "Aside from the Kawauchi-Kojima paper, these comments "
                "use C.T.C. Wall's theorem that 3-manifolds embed in "
                "S<sup>5</sup> and some elementary homological "
                "observations.</qt>"
                );
    label->setWhatsThis(msg);
    EmbeddingComments->setWhatsThis(msg);

    refreshLabels();

    connect(&ReginaPrefSet::global(), SIGNAL(preferencesChanged()),
        this, SLOT(updatePreferences()));
}


regina::NPacket* NTriCellularInfoUI::getPacket() {
    return tri;
}

QWidget* NTriCellularInfoUI::getInterface() {
    return ui;
}

void NTriCellularInfoUI::refreshLabels() {
    if (ReginaPrefSet::global().displayUnicode)
        labelBdryMap->setText(QObject::trUtf8(
            "<qt>H\u2081(\u2202M &rarr; M): </qt>"));
    else
        labelBdryMap->setText(QObject::tr("<qt>H1(bdry M &rarr; M): </qt>"));
}

void NTriCellularInfoUI::updatePreferences() {
    // If we've changed the unicode setting, then we may need some redrawing.
    refreshLabels();
    refresh();
}

