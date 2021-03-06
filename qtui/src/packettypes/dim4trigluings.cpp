
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Qt User Interface                                                     *
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
#include "dim4/dim4triangulation.h"
#include "file/nxmlfile.h"
#include "packet/ncontainer.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "dim4eltmovedialog.h"
#include "dim4trigluings.h"
#include "edittableview.h"
#include "reginamain.h"
#include "reginasupport.h"
#include "choosers/dim4boundarycomponentchooser.h"
#include "choosers/facechooser.h"

#include <memory>
#include <QAction>
#include <QHeaderView>
#include <QLabel>
#include <QMessageBox>
#include <QRegExp>
#include <QToolBar>

using regina::NPacket;
using regina::Dim4Triangulation;

namespace {
    /**
     * Represents a destination for a single facet gluing.
     */
    QRegExp reFacetGluing(
        "^\\s*"
        "(\\d+)"
        "(?:\\s*\\(\\s*|\\s+)"
        "([0-4][0-4][0-4][0-4])"
        "\\s*\\)?\\s*$");

    /**
     * Represents a single pentachoron facet.
     */
    QRegExp reFacet("^[0-4][0-4][0-4][0-4]$");
}

Dim4GluingsModel::Dim4GluingsModel(Dim4Triangulation* tri, bool readWrite) :
        tri_(tri), isReadWrite_(readWrite) {
}

void Dim4GluingsModel::rebuild() {
    beginResetModel();
    endResetModel();
}

QModelIndex Dim4GluingsModel::index(int row, int column,
        const QModelIndex& parent) const {
    return createIndex(row, column, quint32(6 * row + column));
}

int Dim4GluingsModel::rowCount(const QModelIndex& parent) const {
    return tri_->size();
}

int Dim4GluingsModel::columnCount(const QModelIndex& parent) const {
    return 6;
}

QVariant Dim4GluingsModel::data(const QModelIndex& index, int role) const {
    regina::Dim4Pentachoron* p = tri_->simplex(index.row());
    if (role == Qt::DisplayRole) {
        // Pentachoron name?
        if (index.column() == 0)
            return (p->description().empty() ? QString::number(index.row()) :
                (QString::number(index.row()) + " (" +
                p->description().c_str() + ')'));

        // Facet gluing?
        int facet = 5 - index.column();
        if (facet >= 0)
            return destString(facet,
                p->adjacentSimplex(facet),
                p->adjacentGluing(facet));
        return QVariant();
    } else if (role == Qt::EditRole) {
        // Pentachoron name?
        if (index.column() == 0)
            return p->description().c_str();

        // Facet gluing?
        int facet = 5 - index.column();
        if (facet >= 0)
            return destString(facet,
                p->adjacentSimplex(facet),
                p->adjacentGluing(facet));
        return QVariant();
    } else
        return QVariant();
}

QVariant Dim4GluingsModel::headerData(int section, Qt::Orientation orientation,
        int role) const {
    if (orientation != Qt::Horizontal)
        return QVariant();
    if (role != Qt::DisplayRole)
        return QVariant();

    switch (section) {
        case 0: return tr("Pentachoron");
        case 1: return tr("Facet 0123");
        case 2: return tr("Facet 0124");
        case 3: return tr("Facet 0134");
        case 4: return tr("Facet 0234");
        case 5: return tr("Facet 1234");
    }
    return QVariant();
}

Qt::ItemFlags Dim4GluingsModel::flags(const QModelIndex& index) const {
    if (isReadWrite_)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
    else
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

bool Dim4GluingsModel::setData(const QModelIndex& index, const QVariant& value,
        int role) {
    regina::Dim4Pentachoron* p = tri_->simplex(index.row());
    if (index.column() == 0) {
        QString newName = value.toString().trimmed();
        if (newName == p->description().c_str())
            return false;

        p->setDescription(newName.toUtf8().constData());
        return true;
    }

    int facet = 5 - index.column();
    if (facet < 0)
        return false;

    int newAdjPent;
    regina::NPerm5 newAdjPerm;

    // Find the proposed new gluing.
    QString text = value.toString().trimmed();

    if (text.isEmpty()) {
        // Boundary facet.
        newAdjPent = -1;
    } else if (! reFacetGluing.exactMatch(text)) {
        // Bad string.
        showError(tr("<qt>The facet gluing should be entered in the "
            "form: <i>pent (facet)</i>.  An example is <i>6 (0342)</i>, "
            "which represents facet 0342 of pentachoron 6.</qt>"));
        return false;
    } else {
        // Real facet.
        newAdjPent = reFacetGluing.cap(1).toInt();
        QString pentFacet = reFacetGluing.cap(2);

        // Check explicitly for a negative pentachoron number
        // since isFacetStringValid() takes an unsigned integer.
        if (newAdjPent < 0 || newAdjPent >= tri_->size()) {
            showError(tr("There is no pentachoron number %1.").
                arg(newAdjPent));
            return false;
        }

        // Do we have a valid gluing?
        QString err = isFacetStringValid(index.row(), facet, newAdjPent,
            pentFacet, &newAdjPerm);
        if (! err.isNull()) {
            showError(err);
            return false;
        }
    }

    // Yes, looks valid.
    // Have we even made a change?
    if (newAdjPent < 0 && ! p->adjacentSimplex(facet))
        return false;
    if (p->adjacentSimplex(facet) &&
            p->adjacentSimplex(facet)->markedIndex() == newAdjPent &&
            newAdjPerm == p->adjacentGluing(facet))
        return false;

    // Yes!  Go ahead and make the change.
    regina::NPacket::ChangeEventSpan span(tri_);

    // First unglue from the old partner if it exists.
    if (p->adjacentSimplex(facet))
        p->unjoin(facet);

    // Are we making the facet boundary?
    if (newAdjPent < 0)
        return true;

    // We are gluing the facet to a new partner.
    int newAdjFacet = newAdjPerm[facet];

    // Does this new partner already have its own partner?
    // If so, better unglue it.
    regina::Dim4Pentachoron* adj = tri_->simplex(newAdjPent);
    if (adj->adjacentSimplex(newAdjFacet))
        adj->unjoin(newAdjFacet);

    // Glue the two facets together.
    p->join(facet, adj, newAdjPerm);
    return true;
}

QString Dim4GluingsModel::isFacetStringValid(unsigned long srcPent,
        int srcFacet, unsigned long destPent, const QString& destFacet,
        regina::NPerm5* gluing) {
    if (destPent >= tri_->size())
        return tr("There is no pentachoron number %1.").arg(destPent);

    if (! reFacet.exactMatch(destFacet))
        return tr("<qt>%1 is not a valid pentachoron facet.  A pentachoron "
            "facet must be described by a sequence of four vertices, each "
            "between 0 and 4 inclusive.  An example is <i>0342</i>.</qt>").
            arg(destFacet);

    if (destFacet[0] == destFacet[1] || destFacet[0] == destFacet[2] ||
            destFacet[0] == destFacet[3] || destFacet[1] == destFacet[2] ||
            destFacet[1] == destFacet[3] || destFacet[2] == destFacet[3])
        return tr("%1 is not a valid pentachoron facet.  The four vertices "
            "forming the facet must be distinct.").arg(destFacet);

    regina::NPerm5 foundGluing = facetStringToPerm(srcFacet, destFacet);
    if (srcPent == destPent && foundGluing[srcFacet] == srcFacet)
        return tr("A facet cannot be glued to itself.");

    // It's valid!
    if (gluing)
        *gluing = foundGluing;

    return QString::null;
}

void Dim4GluingsModel::showError(const QString& message) {
    // We should actually pass the view to KMessageBox, not 0, but we
    // don't have access to any widget from here...
    ReginaSupport::info(0 /* should be the view? */,
        tr("This is not a valid gluing."), message);
}

QString Dim4GluingsModel::destString(int srcFacet,
        regina::Dim4Pentachoron* destPent,
        const regina::NPerm5& gluing) {
    if (! destPent)
        return "";
    else
        return QString::number(destPent->markedIndex()) + " (" +
            (gluing * regina::Dim4Tetrahedron::ordering(srcFacet)).
            trunc4().c_str() + ')';
}

regina::NPerm5 Dim4GluingsModel::facetStringToPerm(int srcFacet,
        const QString& str) {
    int destVertex[5];

    destVertex[4] = 10; // This will be adjusted in a moment.
    for (int i = 0; i < 4; i++) {
        // Use toLatin1() here because we are converting characters,
        // not strings.
        destVertex[i] = str[i].toLatin1() - '0';
        destVertex[4] -= destVertex[i];
    }

    return regina::NPerm5(destVertex[0], destVertex[1], destVertex[2],
        destVertex[3], destVertex[4]) *
        regina::Dim4Tetrahedron::ordering(srcFacet).inverse();
}

Dim4TriGluingsUI::Dim4TriGluingsUI(regina::Dim4Triangulation* packet,
        PacketTabbedUI* useParentUI, bool readWrite) :
        PacketEditorTab(useParentUI), tri(packet) {
    // Set up the table of facet gluings.
    model = new Dim4GluingsModel(packet, readWrite);
    facetTable = new EditTableView();
    facetTable->setSelectionMode(QAbstractItemView::ContiguousSelection);
    facetTable->setModel(model);

    if (readWrite) {
        QAbstractItemView::EditTriggers flags(
            QAbstractItemView::AllEditTriggers);
        flags ^= QAbstractItemView::CurrentChanged;
        facetTable->setEditTriggers(flags);
    } else
        facetTable->setEditTriggers(QAbstractItemView::NoEditTriggers);

    facetTable->setWhatsThis(tr("<qt>A table specifying which pentachoron "
        "facets are identified with which others.<p>"
        "Pentachora are numbered upwards from 0, and the five vertices of each "
        "pentachoron are numbered 0, 1, 2, 3 and 4.  Each row of the table "
        "represents a single pentachoron, and shows the identifications "
        "for each of its five facets.<p>"
        "As an example, if we are looking at the table cell for facet 0123 of "
        "pentachoron 7, a gluing of <i>6 (0241)</i> shows that "
        "that this facet is identified with facet 0241 of pentachoron 6, in "
        "such a way that vertices 0, 1, 2 and 3 of pentachoron 7 "
        "are mapped to vertices 0, 2, 4 and 1 respectively of pentachoron 6.<p>"
        "To change these identifications, simply type your own gluings into "
        "the table.</qt>"));

    facetTable->verticalHeader()->hide();

    //facetTable->setColumnStretchable(0, true);
    //facetTable->setColumnStretchable(1, true);
    //facetTable->setColumnStretchable(2, true);
    //facetTable->setColumnStretchable(3, true);
    //facetTable->setColumnStretchable(4, true);
    //facetTable->setColumnStretchable(5, true);

    ui = facetTable;

    // Set up the triangulation actions.
    QAction* sep;


    actAddPent = new QAction(this);
    actAddPent->setText(tr("&Add Pent"));
    actAddPent->setIcon(ReginaSupport::regIcon("insert"));
    actAddPent->setToolTip(tr("Add a new pentachoron"));
    actAddPent->setEnabled(readWrite);
    actAddPent->setWhatsThis(tr("Add a new pentachoron to this "
        "triangulation."));
    enableWhenWritable.append(actAddPent);
    triActionList.append(actAddPent);
    connect(actAddPent, SIGNAL(triggered()), this, SLOT(addPent()));

    actRemovePent = new QAction(this);
    actRemovePent->setText(tr("&Remove Pent"));
    actRemovePent->setIcon(ReginaSupport::regIcon("delete"));
    actRemovePent->setToolTip(tr("Remove the currently selected pentachora"));
    actRemovePent->setEnabled(false);
    actRemovePent->setWhatsThis(tr("Remove the currently selected "
        "pentachora from this triangulation."));
    connect(actRemovePent, SIGNAL(triggered()), this,
        SLOT(removeSelectedPents()));
    connect(facetTable->selectionModel(),
        SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)),
        this, SLOT(updateRemoveState()));
    triActionList.append(actRemovePent);

    sep = new QAction(this);
    sep->setSeparator(true);
    triActionList.append(sep);

    actSimplify = new QAction(this);
    actSimplify->setText(tr("&Simplify"));
    actSimplify->setIcon(ReginaSupport::regIcon("simplify"));
    actSimplify->setToolTip(tr(
        "Simplify the triangulation as far as possible"));
    actSimplify->setEnabled(readWrite);
    actSimplify->setWhatsThis(tr("Simplify this triangulation to use fewer "
        "pentachora without changing the underlying 4-manifold or its "
        "PL structure.  This triangulation will be modified directly.<p>"
        "Note that there is no guarantee that the smallest possible number of "
        "pentachora will be achieved."));
    connect(actSimplify, SIGNAL(triggered()), this, SLOT(simplify()));
    enableWhenWritable.append(actSimplify);
    triActionList.append(actSimplify);

    QAction* actEltMove = new QAction(this);
    actEltMove->setText(tr("&Elementary Move..."));
    actEltMove->setToolTip(tr(
        "Select an elementary move with which to modify the triangulation"));
    actEltMove->setEnabled(readWrite);
    actEltMove->setWhatsThis(tr("<qt>Perform an elementary move upon this "
        "triangulation.  <i>Elementary moves</i> are modifications local to "
        "a small number of pentachora that do not change the underlying "
        "4-manifold.<p>"
        "A dialog will be presented in which you can select the precise "
        "elementary move to apply.</qt>"));
    enableWhenWritable.append(actEltMove);
    triActionList.append(actEltMove);
    connect(actEltMove, SIGNAL(triggered()), this, SLOT(elementaryMove()));

    sep = new QAction(this);
    sep->setSeparator(true);
    triActionList.append(sep);

    actOrient = new QAction(this);
    actOrient->setText(tr("&Orient"));
    actOrient->setIcon(ReginaSupport::regIcon("orient"));
    actOrient->setToolTip(tr(
        "Relabel vertices of pentachora for consistent orientation"));
    actOrient->setEnabled(readWrite);
    actOrient->setWhatsThis(tr("<qt>Relabel the vertices of each pentachoron "
        "so that all pentachora are oriented consistently, i.e., "
        "so that orientation is preserved across adjacent facets.<p>"
        "If this triangulation includes both orientable and non-orientable "
        "components, only the orientable components will be relabelled.</qt>"));
    triActionList.append(actOrient);
    connect(actOrient, SIGNAL(triggered()), this, SLOT(orient()));

    QAction* actBarycentricSubdivide = new QAction(this);
    actBarycentricSubdivide->setText(tr("&Barycentric Subdivision"));
    actBarycentricSubdivide->setIcon(ReginaSupport::regIcon("barycentric"));
    actBarycentricSubdivide->setToolTip(tr(
        "Perform a barycentric subdivision"));
    actBarycentricSubdivide->setEnabled(readWrite);
    actBarycentricSubdivide->setWhatsThis(tr("Perform a barycentric "
        "subdivision on this triangulation.  The triangulation will be "
        "changed directly.<p>"
        "This operation involves subdividing each pentachoron into "
        "120 smaller pentachora."));
    enableWhenWritable.append(actBarycentricSubdivide);
    triActionList.append(actBarycentricSubdivide);
    connect(actBarycentricSubdivide, SIGNAL(triggered()), this,
        SLOT(barycentricSubdivide()));

    QAction* actIdealToFinite = new QAction(this);
    actIdealToFinite->setText(tr("&Truncate Ideal Vertices"));
    actIdealToFinite->setIcon(ReginaSupport::regIcon("finite"));

    actIdealToFinite->setToolTip(tr(
        "Truncate any ideal vertices"));
    actIdealToFinite->setEnabled(readWrite);
    actIdealToFinite->setWhatsThis(tr("Convert this from an ideal "
        "triangulation to a finite triangulation.  Any vertices whose "
        "links are neither 3-spheres nor 3-balls "
        "will be truncated and converted into boundary tetrahedra.<p>"
        "This triangulation will be modified directly.  If there are no "
        "vertices of this type to truncate, this operation will have no "
        "effect."));
    enableWhenWritable.append(actIdealToFinite);
    triActionList.append(actIdealToFinite);
    connect(actIdealToFinite, SIGNAL(triggered()), this, SLOT(idealToFinite()));

    QAction* actFiniteToIdeal = new QAction(this);
    actFiniteToIdeal->setText(tr("Make &Ideal"));
    actFiniteToIdeal->setIcon(ReginaSupport::regIcon("cone"));
    actFiniteToIdeal->setToolTip(tr(
        "Convert real boundary components into ideal vertices"));
    actFiniteToIdeal->setEnabled(readWrite);
    actFiniteToIdeal->setWhatsThis(tr("Convert this from a finite "
        "triangulation to an ideal triangulation.  Each real boundary "
        "component (formed from one or more boundary tetrahedra) will be "
        "converted into a single ideal vertex.<p>"
        "A side-effect of this operation is that any spherical boundary "
        "components will be filled in with balls.<p>"
        "This triangulation will be modified directly.  If there are no "
        "real boundary components, this operation will have no effect."));
    enableWhenWritable.append(actFiniteToIdeal);
    triActionList.append(actFiniteToIdeal);
    connect(actFiniteToIdeal, SIGNAL(triggered()), this, SLOT(finiteToIdeal()));

    QAction* actDoubleCover = new QAction(this);
    actDoubleCover->setText(tr("&Double Cover"));
    actDoubleCover->setIcon(ReginaSupport::regIcon("doublecover"));
    actDoubleCover->setToolTip(tr(
        "Convert the triangulation to its orientable double cover"));
    actDoubleCover->setEnabled(readWrite);
    actDoubleCover->setWhatsThis(tr("Convert a non-orientable "
        "triangulation into an orientable double cover.  This triangulation "
        "will be modified directly.<p>"
        "If this triangulation is already orientable, it will simply be "
        "duplicated, resulting in a disconnected triangulation."));
    enableWhenWritable.append(actDoubleCover);
    triActionList.append(actDoubleCover);
    connect(actDoubleCover, SIGNAL(triggered()), this, SLOT(doubleCover()));

    sep = new QAction(this);
    sep->setSeparator(true);
    triActionList.append(sep);

    actBoundaryComponents = new QAction(this);
    actBoundaryComponents->setText(tr("Boundar&y Components..."));
    actBoundaryComponents->setIcon(ReginaSupport::regIcon("boundaries"));
    actBoundaryComponents->setToolTip(tr(
        "Build a 3-manifold triangulation from a boundary component"));
    actBoundaryComponents->setWhatsThis(tr("<qt>Build a 3-manifold "
        "triangulation from a boundary component of this triangulation.<p>"
        "If you select a real boundary component, this will construct "
        "a 3-manifold triangulation from its boundary tetrahedra.  "
        "If you select an ideal boundary component, this will construct "
        "a 3-manifold triangulation from the corresponding vertex link.</qt>"));
    triActionList.append(actBoundaryComponents);
    connect(actBoundaryComponents, SIGNAL(triggered()), this,
        SLOT(boundaryComponents()));

    QAction* actVertexLinks = new QAction(this);
    actVertexLinks->setText(tr("&Vertex Links..."));
    actVertexLinks->setIcon(ReginaSupport::regIcon("vtxlinks"));
    actVertexLinks->setToolTip(tr(
        "Build a 3-manifold triangulation from a vertex link"));
    actVertexLinks->setWhatsThis(tr("<qt>Build a 3-manifold triangulation "
        "from the link of a vertex of this triangulation.<p>"
        "If <i>V</i> is a vertex, then the <i>link</i> of <i>V</i> is the "
        "frontier of a small regular neighbourhood of <i>V</i>.  "
        "The tetrahedra that make up this link sit inside "
        "the pentachoron corners that meet together at <i>V</i>.</qt>"));
    triActionList.append(actVertexLinks);
    connect(actVertexLinks, SIGNAL(triggered()), this, SLOT(vertexLinks()));

    QAction* actSplitIntoComponents = new QAction(this);
    actSplitIntoComponents->setText(tr("E&xtract Components"));
    actSplitIntoComponents->setIcon(ReginaSupport::regIcon("components"));
    actSplitIntoComponents->setToolTip(tr(
        "Form a new triangulation for each disconnected component"));
    actSplitIntoComponents->setWhatsThis(tr("<qt>Split a disconnected "
        "triangulation into its individual connected components.  This "
        "triangulation will not be changed &ndash; each "
        "connected component will be added as a new triangulation beneath "
        "it in the packet tree.<p>"
        "If this triangulation is already connected, this operation will "
        "do nothing.</qt>"));
    triActionList.append(actSplitIntoComponents);
    connect(actSplitIntoComponents, SIGNAL(triggered()), this,
        SLOT(splitIntoComponents()));

    // Tidy up.

    refresh();
}

Dim4TriGluingsUI::~Dim4TriGluingsUI() {
    // Make sure the actions, including separators, are all deleted.

    delete model;
}

const QLinkedList<QAction*>& Dim4TriGluingsUI::getPacketTypeActions() {
    return triActionList;
}

void Dim4TriGluingsUI::fillToolBar(QToolBar* bar) {
    bar->addAction(actAddPent);
    bar->addAction(actRemovePent);
    bar->addSeparator();
    bar->addAction(actSimplify);
    bar->addAction(actOrient);
}

regina::NPacket* Dim4TriGluingsUI::getPacket() {
    return tri;
}

QWidget* Dim4TriGluingsUI::getInterface() {
    return ui;
}

void Dim4TriGluingsUI::refresh() {
    model->rebuild();
    updateActionStates();
}

void Dim4TriGluingsUI::endEdit() {
    facetTable->endEdit();
}

void Dim4TriGluingsUI::setReadWrite(bool readWrite) {
    model->setReadWrite(readWrite);

    if (readWrite) {
        QAbstractItemView::EditTriggers flags(
            QAbstractItemView::AllEditTriggers);
        flags ^= QAbstractItemView::CurrentChanged;
        facetTable->setEditTriggers(flags);
    } else
        facetTable->setEditTriggers(QAbstractItemView::NoEditTriggers);

    QLinkedListIterator<QAction*> it(enableWhenWritable);
    while (it.hasNext())
        (it.next())->setEnabled(readWrite);

    updateRemoveState();
    updateActionStates();
}

void Dim4TriGluingsUI::addPent() {
    endEdit();

    tri->newPentachoron();
}

void Dim4TriGluingsUI::removeSelectedPents() {
    endEdit();

    // Gather together all the pentachora to be deleted.
    QModelIndexList sel = facetTable->selectionModel()->selectedIndexes();
    if (sel.empty()) {
        ReginaSupport::warn(ui, tr("No pentachora are selected to remove."));
        return;
    }

    // Selections are contiguous.
    int first, last;
    first = last = sel.front().row();

    int row, i;
    for (i = 1; i < sel.count(); ++i) {
        row = sel[i].row();
        if (row < first)
            first = row;
        if (row > last)
            last = row;
    }

    // Notify the user that pentachora will be removed.
    QMessageBox msgBox(ui);
    msgBox.setWindowTitle(tr("Question"));
    msgBox.setIcon(QMessageBox::Question);
    if (first == last) {
        msgBox.setText(tr("Pentachoron number %1 will be removed.").arg(first));
        msgBox.setInformativeText(tr("Are you sure?"));
    } else {
        msgBox.setText(
            tr("<qt>%1 pentachora (numbers %2&ndash;%3) will be removed.</qt>")
            .arg(last - first + 1).arg(first).arg(last));
        msgBox.setInformativeText(tr("Are you sure?"));
    }
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Yes);
    if (msgBox.exec() != QMessageBox::Yes)
        return;

    // Off we go!
    if (first == 0 && last == tri->size() - 1)
        tri->removeAllSimplices();
    else {
        regina::NPacket::ChangeEventSpan span(tri);
        for (int i = last; i >= first; --i)
            tri->removeSimplexAt(i);
    }
}

void Dim4TriGluingsUI::simplify() {
    endEdit();

    if (! tri->intelligentSimplify())
        ReginaSupport::info(ui,
            tr("I could not simplify the triangulation further."),
            tr("This does not mean that the triangulation is minimal; it "
            "simply means that I could not find a way of reducing it."));
}

void Dim4TriGluingsUI::orient() {
    endEdit();

    if (tri->isOriented()) {
        ReginaSupport::info(ui, tr("This triangulation is already oriented."));
        return;
    }

    bool hasOr = false;
    for (auto c : tri->components())
        if (c->isOrientable()) {
            hasOr = true;
            break;
        }
    if (! hasOr) {
        ReginaSupport::info(ui,
            tr("This triangulation has no orientable components."),
            tr("Non-orientable components cannot be oriented."));
        return;
    }

    tri->orient();
}

void Dim4TriGluingsUI::barycentricSubdivide() {
    endEdit();

    tri->barycentricSubdivision();
}

void Dim4TriGluingsUI::idealToFinite() {
    endEdit();

    if (tri->isValid() && ! tri->isIdeal())
        ReginaSupport::info(ui,
            tr("This triangulation has no ideal vertices."),
            tr("Only ideal vertices can be truncated."));
    else {
        regina::NPacket::ChangeEventSpan span(tri);
        tri->idealToFinite();
        tri->intelligentSimplify();
    }
}

void Dim4TriGluingsUI::finiteToIdeal() {
    endEdit();

    if (! tri->hasBoundaryFacets())
        ReginaSupport::info(ui,
            tr("This triangulation has no real boundary components."),
            tr("Only real boundary components will be converted into "
            "ideal vertices."));
    else {
        regina::NPacket::ChangeEventSpan span(tri);
        tri->finiteToIdeal();
        tri->intelligentSimplify();
    }
}

void Dim4TriGluingsUI::elementaryMove() {
    endEdit();

    (new Dim4EltMoveDialog(ui, tri))->show();
}

void Dim4TriGluingsUI::doubleCover() {
    endEdit();

    tri->makeDoubleCover();
}

void Dim4TriGluingsUI::boundaryComponents() {
    endEdit();

    if (tri->countBoundaryComponents() == 0)
        ReginaSupport::sorry(ui,
            tr("This triangulation does not have any boundary components."));
    else {
        regina::Dim4BoundaryComponent* chosen =
            Dim4BoundaryComponentDialog::choose(ui, tri, 0 /* filter */,
            tr("Boundary Components"),
            tr("Triangulate which boundary component?"),
            tr("<qt>Regina will triangulate whichever "
                "boundary component you choose.<p>"
                "If you select a real boundary component, this will construct "
                "a 3-manifold triangulation from its boundary tetrahedra.  "
                "If you select an ideal boundary component, this will "
                "construct a 3-manifold triangulation from the corresponding "
                "vertex link.</qt>"));
        if (chosen) {
            regina::NTriangulation* ans = new regina::NTriangulation(
                *chosen->triangulation());
            ans->setLabel(tr("Boundary component %1").arg(
                chosen->index()).toUtf8().constData());
            tri->insertChildLast(ans);
            enclosingPane->getMainWindow()->packetView(ans, true, true);
        }
    }
}

void Dim4TriGluingsUI::vertexLinks() {
    endEdit();

    if (tri->countVertices() == 0)
        ReginaSupport::sorry(ui,
            tr("This triangulation does not have any vertices."));
    else {
        regina::Dim4Vertex* chosen =
            FaceDialog<4, 0>::choose(ui, tri, 0 /* filter */,
            tr("Vertex Links"),
            tr("Triangulate the link of which vertex?"),
            tr("<qt>Regina will triangulate the link of whichever "
                "vertex you choose.<p>"
                "If <i>V</i> is a vertex, then the <i>link</i> of "
                "<i>V</i> is the "
                "frontier of a small regular neighbourhood of <i>V</i>.  "
                "The tetrahedra that make up this link sit inside "
                "the pentachoron corners that meet together at "
                "<i>V</i>.</qt>"));
        if (chosen) {
            regina::NTriangulation* ans = new regina::NTriangulation(
                *chosen->buildLink());
            ans->setLabel(tr("Link of vertex %1").arg(
                chosen->index()).toUtf8().constData());
            tri->insertChildLast(ans);
            enclosingPane->getMainWindow()->packetView(ans, true, true);
        }
    }
}

void Dim4TriGluingsUI::splitIntoComponents() {
    endEdit();

    if (tri->countComponents() == 0)
        ReginaSupport::info(ui,
            tr("This triangulation is empty."),
            tr("It has no components."));
    else if (tri->countComponents() == 1)
        ReginaSupport::info(ui,
            tr("This triangulation is connected."),
            tr("It has only one component."));
    else {
        // If there are already children of this triangulation, insert
        // the new triangulations at a deeper level.
        NPacket* base;
        if (tri->firstChild()) {
            base = new regina::NContainer();
            tri->insertChildLast(base);
            base->setLabel(tri->adornedLabel("Components"));
        } else
            base = tri;

        // Make the split.
        size_t nComps = tri->splitIntoComponents(base);

        // Make sure the new components are visible.
        enclosingPane->getMainWindow()->ensureVisibleInTree(base->firstChild());

        // Tell the user what happened.
        ReginaSupport::info(ui,
            tr("%1 components were extracted.").arg(nComps));
    }
}

void Dim4TriGluingsUI::updateRemoveState() {
    // Are we read-write?
    if (model->isReadWrite())
        actRemovePent->setEnabled(
            ! facetTable->selectionModel()->selectedIndexes().empty());
    else
        actRemovePent->setEnabled(false);
}

void Dim4TriGluingsUI::updateActionStates() {
    if (! model->isReadWrite())
        actOrient->setEnabled(false);
    else if (! tri->isOrientable())
        actOrient->setEnabled(false);
    else
        actOrient->setEnabled(! tri->isOriented());

    actBoundaryComponents->setEnabled(! tri->boundaryComponents().empty());
}

