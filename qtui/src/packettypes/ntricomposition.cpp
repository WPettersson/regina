
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
#include "manifold/nmanifold.h"
#include "subcomplex/naugtrisolidtorus.h"
#include "subcomplex/nblockedsfs.h"
#include "subcomplex/nblockedsfsloop.h"
#include "subcomplex/nblockedsfspair.h"
#include "subcomplex/nblockedsfstriple.h"
#include "subcomplex/nl31pillow.h"
#include "subcomplex/nlayeredchain.h"
#include "subcomplex/nlayeredchainpair.h"
#include "subcomplex/nlayeredlensspace.h"
#include "subcomplex/nlayeredloop.h"
#include "subcomplex/nlayeredsolidtorus.h"
#include "subcomplex/nlayeredsurfacebundle.h"
#include "subcomplex/npillowtwosphere.h"
#include "subcomplex/npluggedtorusbundle.h"
#include "subcomplex/nplugtrisolidtorus.h"
#include "subcomplex/nsatblock.h"
#include "subcomplex/nsatregion.h"
#include "subcomplex/nsnappedball.h"
#include "subcomplex/nsnappedtwosphere.h"
#include "subcomplex/nspiralsolidtorus.h"
#include "subcomplex/nstandardtri.h"
#include "subcomplex/ntxicore.h"
#include "triangulation/nisomorphism.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "ntricomposition.h"
#include "../packetchooser.h"
#include "../packeteditiface.h"
#include "../packetfilter.h"
#include "reginasupport.h"

#include <memory>
#include <QApplication>
#include <QClipboard>
#include <QHeaderView>
#include <QLabel>
#include <QLayout>
#include <QMenu>
#include <QPushButton>
#include <QTextDocument>
#include <QTreeWidgetItem>

using regina::NEdge;
using regina::NPacket;
using regina::NPerm4;
using regina::NSatRegion;
using regina::NTriangulation;

NTriCompositionUI::NTriCompositionUI(regina::NTriangulation* packet,
        PacketTabbedUI* useParentUI) : PacketViewerTab(useParentUI),
        tri(packet), comparingTri(0), components(0), lastComponent(0) {
    // Set up the UI.

    ui = new QWidget();
    QBoxLayout* layout = new QVBoxLayout(ui);

    // Set up the isomorphism tester.
    QBoxLayout* wideIsoArea = new QHBoxLayout();
    layout->addLayout(wideIsoArea);

    QBoxLayout* leftIsoArea = new QVBoxLayout();
    wideIsoArea->addLayout(leftIsoArea, 1);

    QString msg = tr("<qt>Compare this with another triangulation to "
        "see whether the triangulations are isomorphic, or whether one is "
        "isomorphic to a subcomplex of the other.<p>"
        "Select the other triangulation in the drop-down box.  The "
        "relationship (if any) between this and the selected triangulation "
        "will be displayed immediately beneath.<p>"
        "If a relationship is found, the specific isomorphism can be "
        "examined through the <i>Details</i> button.");

    QLabel* label = new QLabel(tr("Isomorphism / subcomplex test:"), ui);
    label->setWhatsThis(msg);
    leftIsoArea->addWidget(label);

    QBoxLayout* isoSelectArea = new QHBoxLayout();
    leftIsoArea->addLayout(isoSelectArea);
    label = new QLabel(tr("Compare with T ="), ui);
    label->setWhatsThis(msg);
    isoSelectArea->addWidget(label);
    isoTest = new PacketChooser(tri->root(),
        new SubclassFilter<NTriangulation>(),
        PacketChooser::ROOT_AS_PACKET, true, 0, ui);
    isoTest->setAutoUpdate(true);
    isoTest->setWhatsThis(msg);
    connect(isoTest, SIGNAL(activated(int)), this, SLOT(updateIsoPanel()));
    isoSelectArea->addWidget(isoTest, 1);
    // isoSelectArea->addStretch(1);

    isoResult = new QLabel(tr("Result:"), ui);
    isoResult->setWhatsThis(msg);
    leftIsoArea->addWidget(isoResult);

    isoView = new QPushButton(ReginaSupport::regIcon("packet_view"),
        tr("Details..."), ui);
    // isoView->setFlat(true);
    isoView->setToolTip(tr("View details of isomorphism"));
    isoView->setWhatsThis(tr("View the details of the isomorphism "
        "(if any) between this and the selected triangulation.  The precise "
        "mapping between tetrahedra and tetrahedron vertices will be "
        "displayed in a separate window."));
    connect(isoView, SIGNAL(clicked()), this, SLOT(viewIsomorphism()));
    wideIsoArea->addWidget(isoView);

    // Add a central divider.
    QFrame* divider = new QFrame(ui);
    divider->setFrameStyle(QFrame::HLine | QFrame::Sunken);
    layout->addWidget(divider);

    // Set up the composition viewer.
    msg = tr("<qt>Displays (i) the precise name of the triangulation "
        "and/or underlying 3-manifold if these can be recognised "
        "immediately, (ii) the isomorphism signature from which the "
        "triangulation can be reconstructed, "
        "(iii) the Callahan-Hildebrand-Weeks dehydration "
        "string if the triangulation supports it, and (iv) the details "
        "of any standard combinatorial structures found within the "
        "triangulation.<p>"
        "You can right-click on any line of text to copy it to the "
        "clipboard.<p>"
        "See the users' handbook for further details on the information "
        "listed here.</qt>");

    label = new QLabel(tr("Triangulation composition:"), ui);
    label->setWhatsThis(msg);
    layout->addWidget(label);

    details = new QTreeWidget(ui);
    details->setHeaderHidden(true);
    details->setAlternatingRowColors(true);
    details->setSelectionMode(QAbstractItemView::SingleSelection);
    details->setWhatsThis(msg);
    layout->addWidget(details, 1);

    editIface = new PacketEditTreeWidgetSingleLine(details);
}

NTriCompositionUI::~NTriCompositionUI() {
    delete editIface;
}

regina::NPacket* NTriCompositionUI::getPacket() {
    return tri;
}

QWidget* NTriCompositionUI::getInterface() {
    return ui;
}

void NTriCompositionUI::refresh() {
    updateIsoPanel();

    details->clear();
    components = lastComponent = 0;

    // Try to identify the 3-manifold.
    std::unique_ptr<regina::NStandardTriangulation> standardTri(
        regina::NStandardTriangulation::isStandardTriangulation(tri));
    if (standardTri.get()) {
        addTopLevelSection(
            tr("Triangulation: ") + standardTri->name().c_str());

        std::unique_ptr<regina::NManifold> manifold(standardTri->manifold());
        if (manifold.get())
            addTopLevelSection(
                tr("3-manifold: ") + manifold->name().c_str());
        else
            addTopLevelSection(tr("3-manifold not recognised"));
    } else
        addTopLevelSection(tr("Triangulation not recognised"));

    // Add the isomorphism signature.
    addTopLevelSection(tr("Isomorphism signature: ") + tri->isoSig().c_str());

    // Offer a dehydration string if we have one.
    std::string dehydration = tri->dehydrate();
    if (! dehydration.empty())
        addTopLevelSection(tr("Dehydration: ") + dehydration.c_str());

    // Look for complete closed triangulations.
    findAugTriSolidTori();
    findL31Pillows();
    findLayeredChainPairs();
    findLayeredLensSpaces();
    findLayeredLoops();
    findPlugTriSolidTori();
    findBlockedTriangulations();

    // Look for interesting surfaces.
    findPillowSpheres();
    findSnappedSpheres();

    // Look for bounded subcomplexes.
    findLayeredSolidTori();
    findSnappedBalls();
    findSpiralSolidTori();

    // Expand so that two levels of children are visible.
    bool foundInnerChildren = false;
    QTreeWidgetItem* rootItem = details->invisibleRootItem();
    for (int i=0; i < rootItem->childCount(); i++ ) {

        QTreeWidgetItem* topChild = rootItem->child(i);
        if (topChild->childCount() > 0) {
            topChild->setExpanded(true);
            foundInnerChildren = true;
        }
    }
    details->setRootIsDecorated(foundInnerChildren);
}

void NTriCompositionUI::packetToBeDestroyed(regina::NPacket*) {
    // Our current isomorphism test triangulation is about to be
    // destroyed.
    isoTest->setCurrentIndex(0); // (i.e., None)
    comparingTri = 0; // Don't unlisten, the packet destructor will do that.
    updateIsoPanel();
}

void NTriCompositionUI::updateIsoPanel() {
    // Update the packet chooser in case things have changed.
    isoTest->refreshContents();

    if (isoTest->selectedPacket() != comparingTri) {
        if (comparingTri)
            comparingTri->unlisten(this);
        comparingTri = dynamic_cast<NTriangulation*>(isoTest->selectedPacket());
        if (comparingTri)
            comparingTri->listen(this);
    }

    // Run the isomorphism tests.
    if (comparingTri) {
        if ((isomorphism = tri->isIsomorphicTo(*comparingTri)).get()) {
            isoResult->setText(tr("Result: Isomorphic (this = T)"));
            isoType = IsIsomorphic;
        } else if ((isomorphism = tri->isContainedIn(*comparingTri)).get()) {
            isoResult->setText(tr("Result: Subcomplex (this < T)"));
            isoType = IsSubcomplex;
        } else if ((isomorphism = comparingTri->isContainedIn(*tri)).get()) {
            isoResult->setText(tr("Result: Subcomplex (T < this)"));
            isoType = IsSupercomplex;
        } else {
            isoResult->setText(tr("Result: No relationship"));
            isoType = NoRelationship;
        }
    } else {
        isomorphism.reset();
        isoResult->setText(tr("Result:"));
        isoType = NoRelationship;
    }

    isoView->setEnabled(isomorphism.get());
}

void NTriCompositionUI::viewIsomorphism() {
    if (isoType == NoRelationship || ! comparingTri)
        return;

    QString title, msg;
    QStringList details;

    details += QString("[%1]  &rarr;  [%2]").
        arg(QString(tri->humanLabel().c_str()).toHtmlEscaped()).
        arg(QString(comparingTri->humanLabel().c_str()).toHtmlEscaped());

    if (isoType == IsIsomorphic) {
        title = tr("Details of the isomorphism between "
            "the two triangulations:");
        msg = tr("<qt>The left hand side refers to this "
            "triangulation; the right hand side refers to the selected "
            "triangulation <i>%1</i>.<p>"
            "Each line represents a single tetrahedron and its four "
            "vertices.").
            arg(QString(comparingTri->humanLabel().c_str()).toHtmlEscaped());

        for (unsigned long i = 0; i < tri->size(); i++)
            details += QString("%1 (0123)  &rarr;  %2 (%3)").
                arg(i).
                arg(isomorphism->tetImage(i)).
                arg(isomorphism->facePerm(i).str().c_str())
                ;
    } else {
        title = tr("Details of the isomorphism by which "
            "one triangulation is contained within the other:");
        msg = tr("<qt>The left "
            "hand side refers to this triangulation; the right hand side "
            "refers to the selected "
            "triangulation <i>%1</i>.<p>"
            "Each line represents a single tetrahedron and its four "
            "vertices.").
            arg(QString(comparingTri->humanLabel().c_str()).toHtmlEscaped());

        if (isoType == IsSubcomplex)
            for (unsigned long i = 0; i < tri->size(); i++)
                details += QString("%1 (0123)  &rarr;  %2 (%3)").
                    arg(i).
                    arg(isomorphism->tetImage(i)).
                    arg(isomorphism->facePerm(i).str().c_str())
                    ;
        else
            for (unsigned long i = 0;
                    i < comparingTri->size(); i++)
                details += QString("%2 (%3)  &rarr;  %1 (0123)").
                    arg(i).
                    arg(isomorphism->tetImage(i)).
                    arg(isomorphism->facePerm(i).str().c_str())
                    ;
    }

    if (details.size() == 1)
        details += tr("(no tetrahedra)");

    // Redo this to actually display information as a list?
    ReginaSupport::info(ui,
        title, msg + "<p>" + details.join("<br>") + "<qt>");
}

QTreeWidgetItem* NTriCompositionUI::addTopLevelSection(const QString& text) {
    QTreeWidgetItem* newItem = new QTreeWidgetItem(details);
    newItem->setText(0,text);
    return newItem;
    //if (details->lastItem())
    //    return new QTreeWidgetItem(details, details->lastItem(), text);
    //else
    //    return new QTreeWidgetItem(details, text);
}

QTreeWidgetItem* NTriCompositionUI::addComponentSection(const QString& text) {
    if (! components)
        components = addTopLevelSection(tr("Components"));

    //if (lastComponent)
    //    lastComponent = new QTreeWidgetItem(components, lastComponent, text);
    //else
    //    lastComponent = new QTreeWidgetItem(components, text);
    
    lastComponent = new QTreeWidgetItem(components);
    lastComponent->setText(0,text);

    return lastComponent;
}

void NTriCompositionUI::findAugTriSolidTori() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NAugTriSolidTorus* aug;
    for (unsigned long i = 0; i < nComps; i++) {
        aug = regina::NAugTriSolidTorus::isAugTriSolidTorus(
            tri->component(i));
        if (aug) {
            id = addComponentSection(tr(
                "Augmented triangular solid torus ") + aug->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0,tr("Component %1").arg(i));

            const regina::NTriSolidTorus& core = aug->core();
            details = new QTreeWidgetItem(id, details);
            details->setText(0,tr("Core: tets %1, %2, %3").
                arg(core.tetrahedron(0)->index()).
                arg(core.tetrahedron(1)->index()).
                arg(core.tetrahedron(2)->index()));

            if (aug->hasLayeredChain()) {
                QString chainType;
                if (aug->chainType() ==
                        regina::NAugTriSolidTorus::CHAIN_MAJOR)
                    chainType = tr("major");
                else if (aug->chainType() ==
                        regina::NAugTriSolidTorus::CHAIN_AXIS)
                    chainType = tr("axis");
                else
                    chainType = tr("unknown");

                details = new QTreeWidgetItem(id, details);
                details->setText(0,tr("Attached: layered chain (%1) + "
                    "layered solid torus").
                    arg(chainType));
            } else {
                details = new QTreeWidgetItem(id, details);
                details->setText(0,tr("Attached: 3 layered solid tori"));
            }

            delete aug;
        }
    }
}

void NTriCompositionUI::describeSatRegion(const NSatRegion& region,
        QTreeWidgetItem* parent) {
    QTreeWidgetItem* details;
    QTreeWidgetItem* annuli;

    regina::NSatBlockSpec spec;
    regina::NSatAnnulus ann;
    unsigned long nAnnuli;
    long a, b;
    bool ref, back;
    QString thisAnnulus, adjAnnulus;
    for (b = region.numberOfBlocks() - 1; b >= 0; b--) {
        spec = region.block(b);
        details = new QTreeWidgetItem(parent);
        details->setText(0,tr("Block %1: %2").
            arg(b).arg(spec.block->abbr().c_str()));

        nAnnuli = spec.block->nAnnuli();

        annuli = new QTreeWidgetItem(details);
        annuli->setText(0,tr("Adjacencies:"));

        for (a = nAnnuli - 1; a >= 0; a--) {
            thisAnnulus = tr("Annulus %1/%2").arg(b).arg(a);
            if (! spec.block->hasAdjacentBlock(a))
                (new QTreeWidgetItem(annuli))->setText(0,
                    tr("%1 --> boundary").arg(thisAnnulus));
            else {
                adjAnnulus = tr("Annulus %1/%2").
                    arg(region.blockIndex(spec.block->adjacentBlock(a))).
                    arg(spec.block->adjacentAnnulus(a));
                ref = spec.block->adjacentReflected(a);
                back = spec.block->adjacentBackwards(a);

                if (ref && back)
                    (new QTreeWidgetItem(annuli))->setText(0,
                        tr("%1 --> %2 (reflected, backwards)").
                        arg(thisAnnulus).arg(adjAnnulus));
                else if (ref)
                    (new QTreeWidgetItem(annuli))->setText(0,
                        tr("%1 --> %2 (reflected)").
                        arg(thisAnnulus).arg(adjAnnulus));
                else if (back)
                    (new QTreeWidgetItem(annuli))->setText(0,
                        tr("%1 --> %2 (backwards)").
                        arg(thisAnnulus).arg(adjAnnulus));
                else
                    (new QTreeWidgetItem(annuli))->setText(0,
                        tr("%1 --> %2").
                        arg(thisAnnulus).arg(adjAnnulus));
            }
        }

        if (nAnnuli == 1) {
            annuli = new QTreeWidgetItem(details);
            annuli->setText(0,tr("1 annulus"));
        } else {
            annuli = new QTreeWidgetItem(details);
            annuli->setText(0,tr("%1 annuli").arg(nAnnuli));
        }
        for (a = nAnnuli - 1; a >= 0; a--) {
            thisAnnulus = tr("Annulus %1/%2").arg(b).arg(a);
            ann = spec.block->annulus(a);

            (new QTreeWidgetItem(annuli))->setText(0,
                tr("%1 : Tet %2 (%3%4%5), Tet %6 (%7%8%9)").
                arg(thisAnnulus).
                arg(ann.tet[0]->index()).
                arg(ann.roles[0][0]).
                arg(ann.roles[0][1]).
                arg(ann.roles[0][2]).
                arg(ann.tet[1]->index()).
                arg(ann.roles[1][0]).
                arg(ann.roles[1][1]).
                arg(ann.roles[1][2]));
        }

        if (spec.refVert && spec.refHoriz)
            (new QTreeWidgetItem(details))->setText(0,
                tr("Reflected vertically and horizontally"));
        else if (spec.refVert)
            (new QTreeWidgetItem(details))->setText(0,
                tr("Reflected vertically"));
        else if (spec.refHoriz)
            (new QTreeWidgetItem(details))->setText(0,
                tr("Reflected horizontally"));
        else
            (new QTreeWidgetItem(details))->setText(0,
                tr("No reflections"));

        (new QTreeWidgetItem(details))->setText(0,
            spec.block->str().c_str());
    }
}

void NTriCompositionUI::findBlockedTriangulations() {
    QTreeWidgetItem* id;
    QTreeWidgetItem* details;

    regina::NBlockedSFS* sfs = regina::NBlockedSFS::isBlockedSFS(tri);
    if (sfs) {
        id = addComponentSection(tr("Blocked Seifert Fibred Space"));
        describeSatRegion(sfs->region(), id);
        delete sfs;
    }

    regina::NBlockedSFSLoop* loop =
        regina::NBlockedSFSLoop::isBlockedSFSLoop(tri);
    if (loop) {
        id = addComponentSection(tr("Blocked SFS Loop"));

        details = new QTreeWidgetItem(id);
        details->setText(0,tr("Internal region:"));
        describeSatRegion(loop->region(), details);

        (new QTreeWidgetItem(id))->setText(0, tr("Matching relation: %1").
            arg(matrixString(loop->matchingReln())));

        delete loop;
    }

    regina::NBlockedSFSPair* pair =
        regina::NBlockedSFSPair::isBlockedSFSPair(tri);
    if (pair) {
        id = addComponentSection(tr("Blocked SFS Pair"));

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("Second region:"));
        describeSatRegion(pair->region(1), details);

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("First region:"));
        describeSatRegion(pair->region(0), details);

        (new QTreeWidgetItem(id))->setText(0, tr("Matching relation (first --> second): %1").
            arg(matrixString(pair->matchingReln())));

        delete pair;
    }

    regina::NBlockedSFSTriple* triple =
        regina::NBlockedSFSTriple::isBlockedSFSTriple(tri);
    if (triple) {
        id = addComponentSection(tr("Blocked SFS Triple"));

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("Second end region:"));
        describeSatRegion(triple->end(1), details);

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("First end region:"));
        describeSatRegion(triple->end(0), details);

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("Central region:"));
        describeSatRegion(triple->centre(), details);

        (new QTreeWidgetItem(id))->setText(0,
            tr("Matching relation (centre --> second end): %1").
            arg(matrixString(triple->matchingReln(1))));

        (new QTreeWidgetItem(id))->setText(0,
            tr("Matching relation (centre --> first end): %1").
            arg(matrixString(triple->matchingReln(0))));

        delete triple;
    }

    regina::NLayeredTorusBundle* bundle =
        regina::NLayeredTorusBundle::isLayeredTorusBundle(tri);
    if (bundle) {
        id = addComponentSection(tr("Layered Torus Bundle"));

        (new QTreeWidgetItem(id))->setText(0,
            tr("Layering relation (lower a/b --> upper a/b): %1").
            arg(matrixString(bundle->layeringReln())));

        (new QTreeWidgetItem(id))->setText(0,
            tr("Core relation (upper a/b --> lower a/b): %1").
            arg(matrixString(bundle->core().parallelReln())));

        (new QTreeWidgetItem(id))->setText(0,
            tr("Core T x I triangulation: %1").
            arg(bundle->core().name().c_str()));

        delete bundle;
    }

    regina::NPluggedTorusBundle* pBundle =
        regina::NPluggedTorusBundle::isPluggedTorusBundle(tri);
    if (pBundle) {
        id = addComponentSection(tr("Plugged Torus Bundle"));

        details = new QTreeWidgetItem(id);
        details->setText(0, tr("Saturated region:"));
        describeSatRegion(pBundle->region(), details);

        (new QTreeWidgetItem(id))->setText(0,
            tr("Matching relation (joining region boundaries): %1").
            arg(matrixString(pBundle->matchingReln())));

        (new QTreeWidgetItem(id))->setText(0,
            tr("Thin I-bundle (T x I): %1").
            arg(pBundle->bundle().name().c_str()));

        delete pBundle;
    }
}

void NTriCompositionUI::findL31Pillows() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NL31Pillow* pillow;
    for (unsigned long i = 0; i < nComps; i++) {
        pillow = regina::NL31Pillow::isL31Pillow(tri->component(i));
        if (pillow) {
            id = addComponentSection(tr("L(3,1) pillow ") +
                pillow->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Component %1").arg(i));

            details = new QTreeWidgetItem(id, details);
            details->setText(0, 
                tr("Pillow interior vertex: %1").
                arg(pillow->tetrahedron(0)->vertex(pillow->interiorVertex(0))->
                    index()));

            delete pillow;
        }
    }
}

void NTriCompositionUI::findLayeredChainPairs() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NLayeredChainPair* pair;
    for (unsigned long i = 0; i < nComps; i++) {
        pair = regina::NLayeredChainPair::isLayeredChainPair(
            tri->component(i));
        if (pair) {
            id = addComponentSection(tr("Layered chain pair ") +
                pair->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Component %1").arg(i));

            details = new QTreeWidgetItem(id, details);
            details->setText(0,
                tr("Chain lengths: %1, %2").
                arg(pair->chain(0)->index()).
                arg(pair->chain(1)->index()));

            delete pair;
        }
    }
}

void NTriCompositionUI::findLayeredLensSpaces() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NLayeredLensSpace* lens;
    for (unsigned long i = 0; i < nComps; i++) {
        lens = regina::NLayeredLensSpace::isLayeredLensSpace(
            tri->component(i));
        if (lens) {
            id = addComponentSection(tr("Layered lens space ") +
                lens->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Component %1").arg(i));

            const regina::NLayeredSolidTorus& torus(lens->torus());
            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr(
                "Layered %1-%2-%3 solid torus %4").
                arg(torus.meridinalCuts(0)).
                arg(torus.meridinalCuts(1)).
                arg(torus.meridinalCuts(2)).
                arg(lens->isSnapped() ? tr("snapped shut") :
                    tr("twisted shut")));

            delete lens;
        }
    }
}

void NTriCompositionUI::findLayeredLoops() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NLayeredLoop* loop;
    for (unsigned long i = 0; i < nComps; i++) {
        loop = regina::NLayeredLoop::isLayeredLoop(tri->component(i));
        if (loop) {
            id = addComponentSection(tr("Layered loop ") +
                loop->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Component %1").arg(i));

            if (loop->isTwisted()) {
                details = new QTreeWidgetItem(id, details);
                details->setText(0, tr(
                    "Length %1, twisted").arg(loop->length()));
                details = new QTreeWidgetItem(id, details);
                details->setText(0, tr(
                    "Hinge: edge %1").arg(loop->hinge(0)->index()));
            } else {
                details = new QTreeWidgetItem(id, details);
                details->setText(0, tr(
                    "Length %1, not twisted").arg(loop->length()));
                details = new QTreeWidgetItem(id);
                details->setText(0, tr(
                    "Hinges: edge %1, %2").
                    arg(loop->hinge(0)->index()).
                    arg(loop->hinge(1)->index()));
            }

            delete loop;
        }
    }
}

void NTriCompositionUI::findLayeredSolidTori() {
    unsigned long nTets = tri->size();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NLayeredSolidTorus* torus;
    unsigned long topIndex;
    for (unsigned long i = 0; i < nTets; i++) {
        torus = regina::NLayeredSolidTorus::formsLayeredSolidTorusBase(
            tri->tetrahedron(i));
        if (torus) {
            id = addComponentSection(tr("Layered solid torus ") +
                torus->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Base: tet %1").arg(torus->base()->index()));
            topIndex = torus->topLevel()->index();
            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Top level: tet %1").
                arg(topIndex));

            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr(
                "Weight %1 edge: %2").arg(torus->meridinalCuts(0)).
                arg(edgeString(topIndex, torus->topEdge(0, 0),
                    torus->topEdge(0, 1))));
            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr(
                "Weight %1 edge: %2").arg(torus->meridinalCuts(1)).
                arg(edgeString(topIndex, torus->topEdge(1, 0),
                    torus->topEdge(1, 1))));
            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr(
                "Weight %1 edge: %2").arg(torus->meridinalCuts(2)).
                arg(edgeString(topIndex, torus->topEdge(2, 0),
                    torus->topEdge(2, 1))));

            delete torus;
        }
    }
}

void NTriCompositionUI::findPillowSpheres() {
    unsigned long nTriangles = tri->countTriangles();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    unsigned long i, j;
    regina::NTriangle* f1;
    regina::NTriangle* f2;
    regina::NPillowTwoSphere* pillow;
    for (i = 0; i < nTriangles; i++) {
        f1 = tri->triangle(i);
        for (j = i + 1; j < nTriangles; j++) {
            f2 = tri->triangle(j);
            pillow = regina::NPillowTwoSphere::formsPillowTwoSphere(f1, f2);
            if (pillow) {
                id = addComponentSection(tr("Pillow 2-sphere"));

                details = new QTreeWidgetItem(id);
                details->setText(0, tr("Triangles: %1, %2").
                    arg(i).arg(j));

                details = new QTreeWidgetItem(id, details);
                details->setText(0, tr(
                    "Equator: edges %1, %2, %3").
                     arg(f1->edge(0)->index()).
                     arg(f1->edge(1)->index()).
                     arg(f1->edge(2)->index()));

                delete pillow;
            }
        }
    }
}

void NTriCompositionUI::findPlugTriSolidTori() {
    unsigned long nComps = tri->countComponents();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NPlugTriSolidTorus* plug;
    const regina::NLayeredChain* chain;
    for (unsigned long i = 0; i < nComps; i++) {
        plug = regina::NPlugTriSolidTorus::isPlugTriSolidTorus(
            tri->component(i));
        if (plug) {
            id = addComponentSection(tr("Plugged triangular solid torus ") +
                plug->name().c_str());

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Component %1").arg(i));

            const regina::NTriSolidTorus& core(plug->core());
            details = new QTreeWidgetItem(id, details);
            details->setText(0,
                tr("Core: tets %1, %2, %3").
                arg(core.tetrahedron(0)->index()).
                arg(core.tetrahedron(1)->index()).
                arg(core.tetrahedron(2)->index()));

            QString lengths(tr("Chain lengths: "));
            for (int j = 0; j < 3; j++) {
                chain = plug->chain(j);
                if (chain)
                    lengths += tr("%1 (%2)").arg(chain->index()).
                        arg(plug->chainType(j) ==
                        regina::NPlugTriSolidTorus::CHAIN_MAJOR ?
                        tr("major") : tr("minor"));
                else
                    lengths += "0";
                if (j < 2)
                    lengths += ", ";
            }
            details = new QTreeWidgetItem(id, details);
            details->setText(0, lengths);

            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Equator type: ") +
                (plug->equatorType() ==
                regina::NPlugTriSolidTorus::EQUATOR_MAJOR ?
                tr("major") : tr("minor")));

            delete plug;
        }
    }
}

void NTriCompositionUI::findSnappedBalls() {
    unsigned long nTets = tri->size();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NSnappedBall* ball;
    for (unsigned long i = 0; i < nTets; i++) {
        ball = regina::NSnappedBall::formsSnappedBall(
            tri->tetrahedron(i));
        if (ball) {
            id = addComponentSection(tr("Snapped 3-ball"));

            details = new QTreeWidgetItem(id);
            details->setText(0, tr("Tetrahedron %1").arg(i));

            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Equator: edge %1%2").
                arg(ball->internalFace(0)).arg(ball->internalFace(1)));

            delete ball;
        }
    }
}

void NTriCompositionUI::findSnappedSpheres() {
    unsigned long nTets = tri->size();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    unsigned long i, j;
    regina::NTetrahedron* t1;
    regina::NTetrahedron* t2;
    regina::NSnappedTwoSphere* sphere;
    for (i = 0; i < nTets; i++) {
        t1 = tri->tetrahedron(i);
        for (j = i + 1; j < nTets; j++) {
            t2 = tri->tetrahedron(j);
            sphere = regina::NSnappedTwoSphere::formsSnappedTwoSphere(t1, t2);
            if (sphere) {
                id = addComponentSection(tr("Snapped 2-sphere"));

                details = new QTreeWidgetItem(id);
                details->setText(0, tr("Tetrahedra: %1, %2").
                    arg(i).arg(j));

                const regina::NSnappedBall* ball = sphere->snappedBall(0);
                details = new QTreeWidgetItem(id, details);
                details->setText(0, tr(
                    "Equator: edge %1").arg(
                    ball->tetrahedron()->edge(ball->equatorEdge())->index()));

                delete sphere;
            }
        }
    }
}

void NTriCompositionUI::findSpiralSolidTori() {
    unsigned long nTets = tri->size();

    QTreeWidgetItem* id = 0;
    QTreeWidgetItem* details = 0;

    regina::NSpiralSolidTorus* spiral;
    regina::NTetrahedron* tet;
    int whichPerm;
    unsigned long i, j;
    for (i = 0; i < nTets; i++) {
        tet = tri->tetrahedron(i);
        for (whichPerm = 0; whichPerm < 24 /* size of S4 */; ++whichPerm) {
            if (NPerm4::S4[whichPerm][0] > NPerm4::S4[whichPerm][3])
                continue;

            spiral = regina::NSpiralSolidTorus::formsSpiralSolidTorus(tet,
                NPerm4::S4[whichPerm]);
            if (! spiral)
                continue;
            if (! spiral->isCanonical(tri)) {
                delete spiral;
                continue;
            }

            // We've got one!
            id = addComponentSection(tr("Spiralled solid torus ") +
                spiral->name().c_str());

            unsigned long spiralTets = spiral->size();

            unsigned long* tetIndex = new unsigned long[spiralTets];
            for (j = 0; j < spiralTets; j++)
                tetIndex[j] = spiral->tetrahedron(j)->index();

            QString tetSet(spiralTets == 1 ? tr("Tet: ") : tr("Tets: "));
            for (j = 0; j < spiralTets; j++) {
                if (j > 0)
                    tetSet += ", ";
                tetSet += QString::number(tetIndex[j]);
            }
            details = new QTreeWidgetItem(id);
            details->setText(0, tetSet);

            QString data;
            QTreeWidgetItem* edge;
            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Major edges:"));
            edge = 0;
            for (j = 0; j < spiralTets; j++) {
                data =
                    edgeString(tetIndex[(j + spiralTets - 1) % spiralTets],
                        spiral->vertexRoles(
                        (j + spiralTets - 1) % spiralTets), 2, 3) +
                    " = " +
                    edgeString(tetIndex[j], spiral->vertexRoles(j), 1, 2) +
                    " = " +
                    edgeString(tetIndex[(j + 1) % spiralTets],
                        spiral->vertexRoles((j + 1) % spiralTets), 0, 1);
                if (edge)
                    edge = new QTreeWidgetItem(details, edge);
                else
                    edge = new QTreeWidgetItem(details);

                edge->setText(0, data);
            }

            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Minor edges:"));
            edge = 0;
            for (j = 0; j < spiralTets; j++) {
                data =
                    edgeString(tetIndex[j], spiral->vertexRoles(j), 1, 3) +
                    " = " +
                    edgeString(tetIndex[(j + 1) % spiralTets],
                        spiral->vertexRoles((j + 1) % spiralTets), 0, 2);
                if (edge)
                    edge = new QTreeWidgetItem(details, edge);
                else
                    edge = new QTreeWidgetItem(details);

                edge->setText(0, data);
            }

            details = new QTreeWidgetItem(id, details);
            details->setText(0, tr("Axis edges:"));
            edge = 0;
            for (j = 0; j < spiralTets; j++) {
                data = edgeString(tetIndex[j], spiral->vertexRoles(j),
                    0, 3);
                if (edge)
                    edge = new QTreeWidgetItem(details, edge);
                else
                    edge = new QTreeWidgetItem(details);

                edge->setText(0, data);
            }

            delete spiral;
        }
    }
}

QString NTriCompositionUI::edgeString(unsigned long tetIndex,
        int edge1, int edge2) {
    if (edge1 < 0)
        return tr("None");
    else if (edge2 < 0)
        return QString("%1 (%2%3)").arg(tetIndex).
            arg(NEdge::edgeVertex[edge1][0]).arg(NEdge::edgeVertex[edge1][1]);
    else
        return QString("%1 (%2%3) = %4 (%5%6)").arg(tetIndex).
            arg(NEdge::edgeVertex[edge1][0]).arg(NEdge::edgeVertex[edge1][1]).
            arg(tetIndex).
            arg(NEdge::edgeVertex[edge2][0]).arg(NEdge::edgeVertex[edge2][1]);
}

QString NTriCompositionUI::edgeString(unsigned long tetIndex,
        const regina::NPerm4& roles, int startPreimage, int endPreimage) {
    return QString("%1 (%2%3)").arg(tetIndex).arg(roles[startPreimage]).
        arg(roles[endPreimage]);
}

QString NTriCompositionUI::matrixString(const regina::NMatrix2& matrix) {
    return QString("[ %1 %2 | %3 %4 ]").
        arg(matrix[0][0]).arg(matrix[0][1]).arg(matrix[1][0]).arg(matrix[1][1]);
}

