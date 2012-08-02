
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

// Regina core includes:
#include "maths/nmatrixint.h"
#include "surfaces/nnormalsurfacelist.h"
#include "surfaces/nsurfacefilter.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "coordinatechooser.h"
#include "nsurfacecoordinateui.h"
#include "reginasupport.h"
#include "packetchooser.h"
#include "packetfilter.h"
#include "reginamain.h"

#include <QHeaderView>
#include <QLabel>
#include <QLayout>
#include <QStyle>
#include <QTreeView>

using regina::NNormalSurfaceList;
using regina::NPacket;

SurfaceModel::SurfaceModel(regina::NNormalSurfaceList* surfaces,
        bool readWrite) :
        surfaces_(surfaces),
        coordSystem_(surfaces->getFlavour()),
        localName(0),
        isReadWrite(readWrite) {
    nFiltered = surfaces_->getNumberOfSurfaces();
    if (nFiltered == 0)
        realIndex = 0;
    else {
        realIndex = new unsigned[nFiltered];
        for (unsigned i = 0; i < nFiltered; ++i)
            realIndex[i] = i;
    }

    refreshNames();
}

void SurfaceModel::rebuild(int coordSystem) {
    beginResetModel();
    coordSystem_ = coordSystem;
    endResetModel();
}

void SurfaceModel::rebuild(int coordSystem, regina::NSurfaceFilter* filter) {
    beginResetModel();

    coordSystem_ = coordSystem;

    delete[] realIndex;
    nFiltered = 0;
    if (surfaces_->getNumberOfSurfaces() == 0)
        realIndex = 0;
    else {
        realIndex = new unsigned[surfaces_->getNumberOfSurfaces()];
        for (unsigned i = 0; i < surfaces_->getNumberOfSurfaces(); ++i)
            if ((! filter) || filter->accept(*surfaces_->getSurface(i)))
                realIndex[nFiltered++] = i;
    }

    endResetModel();
}

void SurfaceModel::refreshNames() {
    delete[] localName;

    if (surfaces_->getNumberOfSurfaces() == 0)
        localName = 0;
    else {
        localName = new QString[surfaces_->getNumberOfSurfaces()];
        for (unsigned i = 0; i < surfaces_->getNumberOfSurfaces(); ++i)
            localName[i] = surfaces_->getSurface(i)->getName().c_str();
    }
}

void SurfaceModel::commitNames() {
    for (unsigned i = 0; i < surfaces_->getNumberOfSurfaces(); ++i)
        const_cast<regina::NNormalSurface*>(surfaces_->getSurface(i))->
            setName(localName[i].toAscii().constData());
}

void SurfaceModel::setReadWrite(bool readWrite) {
    if (readWrite == isReadWrite)
        return;

    // Presumably a model reset is too severe, but I'm not sure what's
    // actually necessary so let's just be safe.
    beginResetModel();
    readWrite = isReadWrite;
    endResetModel();
}

QModelIndex SurfaceModel::index(int row, int column,
        const QModelIndex& parent) const {
    return createIndex(row, column,
        quint32(columnCount(parent) * row + column));
}

int SurfaceModel::rowCount(const QModelIndex& /* unused parent*/) const {
    return nFiltered;
}

int SurfaceModel::columnCount(const QModelIndex& /* unused parent*/) const {
    return propertyColCount() +
        Coordinates::numColumns(coordSystem_, surfaces_->getTriangulation());
}

QVariant SurfaceModel::data(const QModelIndex& index, int role) const {
    unsigned surfaceIndex = realIndex[index.row()];

    if (role == Qt::DisplayRole) {
        const regina::NNormalSurface* s = surfaces_->getSurface(surfaceIndex);

        if (index.column() == 0)
            return tr("%1.").arg(surfaceIndex);
        else if (index.column() == 1)
            return localName[surfaceIndex];
        else if (index.column() == 2) {
            if (! s->isCompact())
                return QVariant();
            return s->getEulerCharacteristic().stringValue().c_str();
        } else if (surfaces_->isEmbeddedOnly() && index.column() == 3) {
            if (! s->isCompact())
                return QVariant();

            if (s->isOrientable())
                return QString(QChar(0x2713 /* tick */));
                // return tr("Yes");
            else
                return tr("Non-or.");
        } else if (surfaces_->isEmbeddedOnly() && index.column() == 4) {
            if (! s->isCompact())
                return QVariant();

            if (s->isTwoSided())
                return "2";
            else
                return "1";
        } else if ((surfaces_->isEmbeddedOnly() && index.column() == 5) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 3)) {
            if (! s->isCompact()) {
                regina::NMatrixInt* slopes = s->boundarySlopes();
                if (slopes) {
                    QString ans = tr("Spun:");
                    // Display each boundary slope as (nu(L), -nu(M)).
                    for (unsigned i = 0; i < slopes->rows(); ++i)
                        ans += QString(" (%1, %2)")
                            .arg(slopes->entry(i,1).stringValue().c_str())
                            .arg((-slopes->entry(i,0)).stringValue().c_str());
                    return ans;
                } else
                    return tr("Spun");
            } else if (s->hasRealBoundary())
                return tr("Real");
            else
                return QString(QChar(0x2014 /* emdash */));
        } else if ((surfaces_->isEmbeddedOnly() && index.column() == 6) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 4)) {
            const regina::NVertex* v;
            std::pair<const regina::NEdge*, const regina::NEdge*> e;

            if ((v = s->isVertexLink()))
                return tr("Vertex %1").arg(
                    surfaces_->getTriangulation()->vertexIndex(v));
            else if ((e = s->isThinEdgeLink()).first) {
                if (e.second)
                    return tr("Thin edges %1, %2").
                        arg(surfaces_->getTriangulation()->edgeIndex(
                            e.first)).
                        arg(surfaces_->getTriangulation()->edgeIndex(
                            e.second));
                else
                    return tr("Thin edge %1").
                        arg(surfaces_->getTriangulation()->edgeIndex(
                            e.first));
            } else
                return QVariant();
        } else if ((surfaces_->isEmbeddedOnly() && index.column() == 7) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 5)) {
            regina::NLargeInteger tot;
            if (s->isSplitting())
                return tr("Splitting");
            else if ((tot = s->isCentral()) != 0)
                return tr("Central (%1)").arg(tot.longValue());
            else
                return QVariant();
        } else if (surfaces_->allowsAlmostNormal() &&
                ((surfaces_->isEmbeddedOnly() && index.column() == 8) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 6))) {
            regina::NDiscType oct = s->getOctPosition();
            if (oct == regina::NDiscType::NONE)
                return QVariant();
            else {
                regina::NLargeInteger tot = s->getOctCoord(
                    oct.tetIndex, oct.type);
                if (tot == 1) {
                    return tr("K%1: %2 (1 oct)").
                        arg(oct.tetIndex).
                        arg(regina::vertexSplitString[oct.type]);
                } else {
                    return tr("K%1: %2 (%3 octs)").
                        arg(oct.tetIndex).
                        arg(regina::vertexSplitString[oct.type]).
                        arg(tot.stringValue().c_str());
                }
            }
        } else {
            // The default case:
            regina::NLargeInteger ans = Coordinates::getCoordinate(coordSystem_,
                *s, index.column() - propertyColCount());
            if (ans == (long)0)
                return QVariant();
            else if (ans.isInfinite())
                return QString(QChar(0x221e /* &infin; */));
            else
                return ans.stringValue().c_str();
        }
    } else if (role == Qt::EditRole) {
        if (index.column() == 1)
            return localName[surfaceIndex];
        else
            return QVariant();
    } else if (role == Qt::ToolTipRole) {
        int propertyCols = propertyColCount();

        if (index.column() < propertyCols)
            return propertyColDesc(index.column());
        else
            return Coordinates::columnDesc(coordSystem_,
                index.column() - propertyCols, this, surfaces_->getTriangulation());
    } else if (role == Qt::ForegroundRole) {
        const regina::NNormalSurface* s = surfaces_->getSurface(surfaceIndex);

        if (surfaces_->isEmbeddedOnly() && index.column() == 3) {
            if (! s->isCompact())
                return QVariant();

            if (s->isOrientable())
                return Qt::darkGreen;
            else
                return Qt::darkRed;
        } else if (surfaces_->isEmbeddedOnly() && index.column() == 4) {
            if (! s->isCompact())
                return QVariant();

            if (s->isTwoSided())
                return Qt::darkGreen;
            else
                return Qt::darkRed;
        } else if ((surfaces_->isEmbeddedOnly() && index.column() == 5) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 3)) {
            if (! s->isCompact())
                return Qt::darkYellow;
            else if (s->hasRealBoundary())
                return Qt::darkRed;
            else
                return Qt::darkGreen;
        } else if (surfaces_->allowsAlmostNormal() &&
                ((surfaces_->isEmbeddedOnly() && index.column() == 8) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 6))) {
            regina::NDiscType oct = s->getOctPosition();
            if (oct != regina::NDiscType::NONE) {
                if (s->getOctCoord(oct.tetIndex, oct.type) > 1)
                    return Qt::darkRed;
                else
                    return Qt::darkGreen;
            }
        } else
            return QVariant();
    } else if (role == Qt::TextAlignmentRole) {
        if (index.column() == 1)
            return Qt::AlignLeft; // Name
        else if (surfaces_->isEmbeddedOnly() && index.column() == 3)
            return Qt::AlignLeft; // Orientable
        else if ((surfaces_->isEmbeddedOnly() && index.column() == 5) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 3))
            return Qt::AlignLeft; // Boundary
        else if ((surfaces_->isEmbeddedOnly() && index.column() == 6) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 4))
            return Qt::AlignLeft; // Vertex / edge link
        else if ((surfaces_->isEmbeddedOnly() && index.column() == 7) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 5))
            return Qt::AlignLeft; // Splitting / central
        else if (surfaces_->allowsAlmostNormal() &&
                ((surfaces_->isEmbeddedOnly() && index.column() == 8) ||
                ((! surfaces_->isEmbeddedOnly()) && index.column() == 6)))
            return Qt::AlignLeft; // Octahedron
        else
            return Qt::AlignRight; // All other fields are numbers.
    }

    return QVariant();
}

QVariant SurfaceModel::headerData(int section, Qt::Orientation orientation,
        int role) const {
    if (orientation != Qt::Horizontal)
        return QVariant();

    if (role == Qt::DisplayRole) {
        int propCols = propertyColCount();

        if (section < propCols)
            return propertyColName(section);
        else
            return Coordinates::columnName(coordSystem_, section - propCols,
                surfaces_->getTriangulation());
    } else if (role == Qt::ToolTipRole) {
        int propertyCols = propertyColCount();

        if (section < propertyCols)
            return propertyColDesc(section);
        else
            return Coordinates::columnDesc(coordSystem_, section - propertyCols,
                this, surfaces_->getTriangulation());
    } else if (role == Qt::TextAlignmentRole)
        return Qt::AlignCenter;
    else
        return QVariant();
}

Qt::ItemFlags SurfaceModel::flags(const QModelIndex& index) const {
    if (index.column() == 1 && isReadWrite)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
    else
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

bool SurfaceModel::setData(const QModelIndex& index, const QVariant& value,
        int role) {
    if (index.column() == 1 && role == Qt::EditRole) {
        localName[realIndex[index.row()]] = value.toString();
        emit dataChanged(index, index);
        return true;
    } else
        return false;
}

unsigned SurfaceModel::propertyColCount() const {
    return (surfaces_->isEmbeddedOnly() ? 8 : 6) +
        (surfaces_->allowsAlmostNormal() ? 1 : 0);
}

QString SurfaceModel::propertyColName(int whichCol) const {
    if (surfaces_->isEmbeddedOnly()) {
        switch (whichCol) {
            case 0 : return QString(); // Surface number
            case 1 : return tr("Name");
            case 2 : return tr("Euler");
            case 3 : return tr("Orient");
            case 4 : return tr("Sides");
            case 5 : return tr("Bdry");
            case 6 : return tr("Link");
            case 7 : return tr("Type");
        }
        if (whichCol == 8 && surfaces_->allowsAlmostNormal())
            return tr("Octagon");
    } else {
        switch (whichCol) {
            case 0 : return QString(); // Surface number
            case 1 : return tr("Name");
            case 2 : return tr("Euler");
            case 3 : return tr("Bdry");
            case 4 : return tr("Link");
            case 5 : return tr("Type");
        }
        if (whichCol == 6 && surfaces_->allowsAlmostNormal())
            return tr("Octagon");
    }

    return tr("Unknown");
}

QString SurfaceModel::propertyColDesc(int whichCol) const {
    if (surfaces_->isEmbeddedOnly()) {
        switch (whichCol) {
            case 0: return tr("The index of this surface within the "
                "overall list (surfaces are numbered 0,1,2,...)");
            case 1: return tr("Name (this has no special meaning and "
                "can be edited)");
            case 2: return tr("Euler characteristic");
            case 3: return tr("Is this surface orientable?");
            case 4: return tr("1-sided or 2-sided");
            case 5: return tr("What kind of boundary does this "
                "surface have?");
            case 6: return tr("Has this surface been identified as "
                "the link of a particular subcomplex?");
            case 7: return tr("Other interesting properties");
        }
        if (whichCol == 8 && surfaces_->allowsAlmostNormal())
            return tr("The coordinate position containing the octagonal "
                "disc type, and the number of discs of that type");
    } else {
        switch (whichCol) {
            case 0: return tr("The index of this surface within the "
                "overall list (surfaces are numbered 0,1,2,...)");
            case 1: return tr("Name (this has no special meaning and "
                "can be edited)");
            case 2: return tr("Euler characteristic");
            case 3: return tr("What kind of boundary does this "
                "surface have?");
            case 4: return tr("Has this surface been identified as "
                "the link of a particular subcomplex?");
            case 5: return tr("Other interesting properties");
        }
        if (whichCol == 6 && surfaces_->allowsAlmostNormal())
            return tr("The coordinate position containing the octagonal "
                "disc type, and the number of discs of that type");
    }

    return tr("Unknown");
}

NSurfaceCoordinateUI::NSurfaceCoordinateUI(regina::NNormalSurfaceList* packet,
        PacketTabbedUI* useParentUI, bool readWrite) :
        PacketEditorTab(useParentUI), surfaces(packet), appliedFilter(0),
        isReadWrite(readWrite), currentlyResizing(false) {
    // Set up the UI.
    ui = new QWidget();
    QBoxLayout* uiLayout = new QVBoxLayout(ui);
    uiLayout->setContentsMargins(0, 0, 0, 0);
    uiLayout->addSpacing(5);

    QBoxLayout* hdrLayout = new QHBoxLayout();
    uiLayout->addLayout(hdrLayout);

    // Set up the coordinate selector.
    QLabel* label = new QLabel(tr("Display coordinates:"));
    hdrLayout->addWidget(label);
    coords = new CoordinateChooser();
    coords->insertAllViewers(surfaces);
    coords->setCurrentSystem(surfaces->getFlavour());
    connect(coords, SIGNAL(activated(int)), this, SLOT(refreshLocal()));
    hdrLayout->addWidget(coords);
    QString msg = tr("Allows you to view these normal surfaces in a "
        "different coordinate system.");
    label->setWhatsThis(msg);
    coords->setWhatsThis(msg);

    hdrLayout->addStretch(1);

    // Set up the filter selector.
    label = new QLabel(tr("Apply filter:"));
    hdrLayout->addWidget(label);
    filter = new PacketChooser(surfaces->getTreeMatriarch(),
        new SingleTypeFilter<regina::NSurfaceFilter>(), true, 0, ui);
    filter->setAutoUpdate(true);
    connect(filter, SIGNAL(activated(int)), this, SLOT(refreshLocal()));
    hdrLayout->addWidget(filter);
    msg = tr("<qt>Allows you to filter this list so that only normal "
        "surfaces satisfying particular properties are displayed.<p>"
        "To use this feature you need a separate surface filter.  You "
        "can create new surface filters through the <i>Packet Tree</i> "
        "menu.</qt>");
    label->setWhatsThis(msg);
    filter->setWhatsThis(msg);

    // Set up the coordinate table.
    model = new SurfaceModel(surfaces, readWrite);

    table = new QTreeView();
    table->setItemsExpandable(false);
    table->setRootIsDecorated(false);
    table->setAlternatingRowColors(true);
    table->header()->setStretchLastSection(false);
    table->setSelectionMode(QTreeView::SingleSelection);
    table->setWhatsThis(tr("Displays details of the individual normal "
        "surfaces in this list.<p>"
        "Each row represents a single normal (or almost normal) surface.  "
        "As well as various properties of the surface, each row contains "
        "a detailed representation the surface in the currently selected "
        "coordinate system.<p>"
        "For details on what each property means or what each coordinate "
        "represents, hover the mouse over the column header (or refer "
        "to the users' handbook).</qt>"));
    // Add grid lines:
    // For some reason, when we set cell borders, the selection
    // background changes to white (-> white on white, which is unreadable).
    // So reluctantly we also break the native system style and explicitly
    // set a background ourselves.  I would love a better way of dealing
    // with this. :/
    table->setStyleSheet("QTreeView::item:selected { "
                            "background: qlineargradient(x1: 0, y1: 0, "
                                "x2: 0, y2: 1, stop: 0 #6ea1f1, "
                                "stop: 1 #567dbc); "
                         "} "
                         "QTreeView::item { "
                            "border: 1px solid #d9d9d9; "
                            "border-top-color: transparent; "
                            "border-left-color: transparent; "
                         "}");
    table->setModel(model);
    table->header()->resizeSections(QHeaderView::ResizeToContents);
    uiLayout->addWidget(table, 1);

    connect(table->header(), SIGNAL(sectionResized(int, int, int)),
        this, SLOT(columnResized(int, int, int)));

    actCutAlong = new QAction(this);
    actCutAlong->setText(tr("Cu&t Along Surface"));
    actCutAlong->setToolTip(tr("Cut the triangulation along the "
        "selected surface"));
    actCutAlong->setEnabled(false);
    actCutAlong->setWhatsThis(tr("<qt>Cuts open the surround triangulation "
        "along the selected surface.  This triangulation will not "
        "be changed; instead a new cut-open triangulation will be created.<p>"
        "This operation will never change the topology of the underlying "
        "3-manifold beyond just cutting along the surface (as opposed to "
        "the related <i>crushing</i> operation, which might).  However, "
        "because the new surface boundaries are created from real "
        "boundary faces, the resulting number of tetrahedra might be very "
        "large.</qt>"));
    connect(actCutAlong, SIGNAL(triggered()), this, SLOT(cutAlong()));
    surfaceActionList.append(actCutAlong);

    actCrush = new QAction(this);
    actCrush->setText("Crus&h Surface");
    actCrush->setToolTip(tr("Crush the selected surface to a point"));
    actCrush->setEnabled(false);
    actCrush->setWhatsThis(tr("<qt>Crushes the selected surface to a point "
        "within the surrounding triangulation.  This triangulation will not "
        "be changed; instead a new crushed triangulation will be created.<p>"
        "<b>Warning:</b> This routine simply removes all tetrahedra "
        "containing quadrilateral discs and rejoins the others "
        "appropriately.  In some circumstances this might change the "
        "topology of the underlying 3-manifold beyond just slicing along "
        "the surface and shrinking the resulting boundary/boundaries "
        "to points.</qt>"));
    surfaceActionList.append(actCrush);
    connect(actCrush, SIGNAL(triggered()), this, SLOT(crush()));

    connect(table->selectionModel(),
        SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)),
        this, SLOT(updateActionStates()));
    connect(model, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)),
        this, SLOT(notifySurfaceRenamed()));

    // Tidy up.
    refresh();
}

NSurfaceCoordinateUI::~NSurfaceCoordinateUI() {
    // Make sure the actions, including separators, are all deleted.
    surfaceActionList.clear();

    delete model;
}

const QLinkedList<QAction*>& NSurfaceCoordinateUI::getPacketTypeActions() {
    return surfaceActionList;
}

regina::NPacket* NSurfaceCoordinateUI::getPacket() {
    return surfaces;
}

QWidget* NSurfaceCoordinateUI::getInterface() {
    return ui;
}

void NSurfaceCoordinateUI::commit() {
    model->commitNames();
    setDirty(false);
}

void NSurfaceCoordinateUI::refreshLocal() {
    // Update the current filter.
    filter->refreshContents();

    bool filterChanged = false;
    if (filter->selectedPacket() != appliedFilter) {
        filterChanged = true;
        if (appliedFilter)
            appliedFilter->unlisten(this);
        appliedFilter = dynamic_cast<regina::NSurfaceFilter*>(
            filter->selectedPacket());
        if (appliedFilter)
            appliedFilter->listen(this);
    }

    // Rebuild the underlying data model.
    int selectedSystem = coords->getCurrentSystem();
    bool coordsChanged = (model->coordSystem() != selectedSystem);
    if (filterChanged)
        model->rebuild(selectedSystem, appliedFilter);
    else
        model->rebuild(selectedSystem); // Faster if the filter is the same.

    // Tidy up.
    updateActionStates();
    if (coordsChanged) {
        currentlyResizing = true;
        table->header()->resizeSections(QHeaderView::ResizeToContents);
        currentlyResizing = false;
    }
}

void NSurfaceCoordinateUI::refresh() {
    // Refresh the surface names from the underlying packet.
    model->refreshNames();

    // Refresh the table of surfaces.
    refreshLocal();

    setDirty(false);
}

void NSurfaceCoordinateUI::setReadWrite(bool readWrite) {
    isReadWrite = readWrite;

    model->setReadWrite(readWrite);
    updateActionStates();
}

void NSurfaceCoordinateUI::packetToBeDestroyed(NPacket*) {
    // Our currently applied filter is about to be destroyed.
    filter->setCurrentIndex(0); // (i.e., None)
    refreshLocal();
}

void NSurfaceCoordinateUI::cutAlong() {
    if (table->selectionModel()->selectedIndexes().empty()) {
        ReginaSupport::info(ui,
            tr("Please select a normal surface to cut along."));
        return;
    }

    const regina::NNormalSurface* toCutAlong =
        model->surface(table->selectionModel()->selectedIndexes().front());
    if (! toCutAlong->isCompact()) {
        ReginaSupport::info(ui,
            tr("I can only cut along compact surfaces."),
            tr("The surface you have selected is non-compact."));
        return;
    }

    // Go ahead and cut along the surface.
    // Be nice and simplify the triangulation, which could be very large.
    regina::NTriangulation* ans = toCutAlong->cutAlong();
    ans->intelligentSimplify();
    ans->setPacketLabel(surfaces->makeUniqueLabel(tr("Cut-open %1").arg(
        surfaces->getTriangulation()->getPacketLabel().c_str())
          .toAscii().constData()));
    surfaces->insertChildLast(ans);

    enclosingPane->getMainWindow()->packetView(ans, true, true);
}

void NSurfaceCoordinateUI::crush() {
    if (table->selectionModel()->selectedIndexes().empty()) {
        ReginaSupport::info(ui,
            tr("Please select a normal surface to crush."));
        return;
    }

    const regina::NNormalSurface* toCrush =
        model->surface(table->selectionModel()->selectedIndexes().front());
    if (! toCrush->isCompact()) {
        ReginaSupport::info(ui,
            tr("I can only crush compact surfaces."),
            tr("The surface you have selected is non-compact."));
        return;
    }

    // Go ahead and crush it.
    regina::NTriangulation* ans = toCrush->crush();
    ans->setPacketLabel(surfaces->makeUniqueLabel(tr("Crushed %1").arg(
        surfaces->getTriangulation()->getPacketLabel().c_str())
          .toAscii().constData()));
    surfaces->insertChildLast(ans);

    enclosingPane->getMainWindow()->packetView(ans, true, true);
}

void NSurfaceCoordinateUI::updateActionStates() {
    bool canCrushOrCut = isReadWrite &&
        table->selectionModel()->hasSelection() &&
        (! surfaces->allowsAlmostNormal()) && surfaces->isEmbeddedOnly();

    actCutAlong->setEnabled(canCrushOrCut);
    actCrush->setEnabled(canCrushOrCut);
}

void NSurfaceCoordinateUI::columnResized(int section, int, int newSize) {
    int nNonCoordSections = model->propertyColCount();
    if (currentlyResizing || section < nNonCoordSections)
        return;

    // A coordinate column has been resized.
    // Resize all coordinate columns.
    currentlyResizing = true;
    for (int i = nNonCoordSections; i < model->columnCount(QModelIndex()); i++)
        table->setColumnWidth(i, newSize);
    currentlyResizing = false;
}

void NSurfaceCoordinateUI::notifySurfaceRenamed() {
    setDirty(true);
}

