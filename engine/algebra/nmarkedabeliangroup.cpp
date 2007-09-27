
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2006, Ben Burton                                   *
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

#include "algebra/nmarkedabeliangroup.h"
#include "maths/matrixops.h"
#include "file/nfile.h"
#include <iostream>

namespace regina {

unsigned long rbGetRank(const NMatrixInt& M) // I don't know how to avoid
                                             // using this, unfortunately.
{ // this is needed in NMarkedAbelianGroup::NMarkedAbelianGroup
    NMatrixInt temp(M);
    smithNormalForm(temp);
    unsigned long i=0;
    while ( (i<temp.rows()) && (i<temp.columns()) && (temp.entry(i,i)!=0) )
        i++;
    return i; // i is the rank of M
}

NMarkedAbelianGroup::NMarkedAbelianGroup(const NMatrixInt& M,
        const NMatrixInt& N) :
        OM(M), ON(N), OMR(M.columns(),M.columns()),
        OMC(M.rows(),M.rows()), OMRi(M.columns(),M.columns()),
        OMCi(M.rows(),M.rows()),
        rankOM(rbGetRank(M)),
        ORN(N.rows()-rbGetRank(M),N.columns()),
        ornR(N.columns(),N.columns()),        ornRi(N.columns(),N.columns()),
        ornC(N.rows()-rbGetRank(M),N.rows()-rbGetRank(M)),
        ornCi(N.rows()-rbGetRank(M),N.rows()-rbGetRank(M)),
        SNF_ORN(N.rows()-rbGetRank(M),N.columns()),
        InvFacList(0), InvFacIndex(0),
        snfrank(0), snffreeindex(0), ifNum(0), ifLoc(0) {
    // find SNF(M).
    NMatrixInt tM(M);

    RBMOD_smithNormalForm(tM, OMR, OMRi, OMC, OMCi);
    // now construct OMRi * N, and delete first SNF_OM_firstzero rows,
    // constructing ORN.

    NMatrixRing<NLargeInteger>* prod=OMRi*ON;

    unsigned long i;
    unsigned long j;

    for (i=0;i<ORN.rows();i++)
        for (j=0;j<ORN.columns();j++)
            ORN.entry(i,j) = prod->entry(i+rankOM,j);

    delete(prod);

    NMatrixInt tX(ORN);
    RBMOD_smithNormalForm(tX, ornR, ornRi, ornC, ornCi);
    for (i=0;i<tX.rows();i++)
        for (j=0;j<tX.columns();j++)
            SNF_ORN.entry(i,j)=tX.entry(i,j);
    // now build the list of invariant factors and their row indexes
    // now compute the rank and column indexes ...
    i=0;
    unsigned long totO=0; // number diag entries == 1
    unsigned long totIF=0;// number diag entries > 1
    unsigned long totFR=0;// number diag entries == 0

    while ((i<SNF_ORN.rows()) && (i<SNF_ORN.columns())) {
        if (SNF_ORN.entry(i,i)==1)
            totO++;
        else if (SNF_ORN.entry(i,i)>1) {
            totIF++;
            InvFacIndex.push_back(i);
        } else
            totFR++;
        i++;
    }

    ifNum=totIF;
    ifLoc=totO;

    InvFacList.resize(InvFacIndex.size());
    for (i=0;i<InvFacList.size();i++)
        InvFacList[i]=SNF_ORN.entry(InvFacIndex[i],InvFacIndex[i]);

    snfrank=SNF_ORN.rows()-totO-totIF;
    snffreeindex=totO+totIF;
}


// this goes through the list of invariant factors in order, stops
// at the `index'-th one in increasing order and returns it...
// if the index is out of bounds, it returns zero.
// invariant factors are stored in the matrix SNF_ORN so one has to
// get past the entries labelled with 1, and avoid any possible tail
// entries labelled with 0.
/** Gives the index-th invariant factor, in increasing order.
        returns 0 for an invalid index. */
const NLargeInteger& NMarkedAbelianGroup::getInvariantFactor(
        unsigned long index) const {
    return InvFacList[index];
}

unsigned long NMarkedAbelianGroup::getNumberOfInvariantFactors() const {
    return InvFacList.size();
}

unsigned NMarkedAbelianGroup::getTorsionRank(const NLargeInteger& degree) const {
    unsigned ans = 0;
    for (unsigned long i=0;i<InvFacList.size();i++) {
        if (InvFacList[i] % degree == 0)
            ans++;
    }
    return ans;
}


void NMarkedAbelianGroup::writeTextShort(std::ostream& out) const {
    bool writtenSomething = false;

    if (snfrank > 0) {
        if (snfrank > 1)
            out << snfrank << ' ';
        out << 'Z';
        writtenSomething = true;
    }

    std::vector<NLargeInteger>::const_iterator it = InvFacList.begin();
    NLargeInteger currDegree;
    unsigned currMult = 0;
    while(true) {
        if (it != InvFacList.end()) {
            if ((*it) == currDegree) {
                currMult++;
                it++;
                continue;
            }
        }
        if (currMult > 0) {
            if (writtenSomething)
                out << " + ";
            if (currMult > 1)
                out << currMult << ' ';
            out << "Z_" << currDegree.stringValue();
            writtenSomething = true;
        }
        if (it == InvFacList.end())
            break;
        currDegree = *it;
        currMult = 1;
        it++;
    }

    if (! writtenSomething)
        out << '0';
}

unsigned NMarkedAbelianGroup::getRank() const {
    return snfrank;
}

bool NMarkedAbelianGroup::isTrivial() const {
    return ! ( (snfrank>0) || (InvFacList.size()>0) );
}

bool NMarkedAbelianGroup::operator == (const NMarkedAbelianGroup& other) const {
    return ((InvFacList == other.InvFacList) && (snfrank==other.snfrank));
}


/*
 * The marked abelian group was defined by matrices M and N
 * with M*N==0.  Think of M as m by l and N as l by n.  Then
 * this routine returns the index-th free generator of the
 * ker(M)/img(N) in Z^l.
 */
std::vector<NLargeInteger> NMarkedAbelianGroup::getFreeRep(unsigned long index)
        const {
    std::vector<NLargeInteger> retval(OM.columns(),"0");
    // index corresponds to the (index+snffreeindex)-th column of ornCi
    // we then pad this vector (at the front) with rankOM 0's
    // and apply OMR to it.

    std::vector<NLargeInteger> temp(ornCi.rows()+rankOM,"0");

    for (unsigned long i=0;i<ornCi.rows();i++)
        temp[i+rankOM]=ornCi.entry(i,index+snffreeindex);
    // the above line takes the index+snffreeindex-th column of ornCi and
    // pads it.

    for (unsigned long i=0;i<retval.size();i++)
        for (unsigned long j=0;j<OMR.columns();j++)
            retval[i] += OMR.entry(i,j)*temp[j];
    // the above takes temp and multiplies it by the matrix OMR.

    return retval;
}


/*
 * The marked abelian group was defined by matrices M and N
 * with M*N==0.  Think of M as m by l and N as l by n.  Then
 * this routine returns the index-th torsion generator of the
 * ker(M)/img(N) in Z^l.
 */
std::vector<NLargeInteger> NMarkedAbelianGroup::getTorRep(unsigned long index)
        const {
    std::vector<NLargeInteger> retval(OM.columns(),NLargeInteger::zero);
    // index corresponds to the (InvFacIndex[index])-th column of ornCi
    // we then pad this vector (at the front) with rankOM 0's
    // and apply OMR to it.

    std::vector<NLargeInteger> temp(ornCi.rows()+rankOM,NLargeInteger::zero);

    for (unsigned long i=0;i<ornCi.rows();i++)
        temp[i+rankOM]=ornCi.entry(i,InvFacIndex[index]);
    // the above line takes the index+snffreeindex-th column of ornCi and
    // pads it.

    for (unsigned long i=0;i<retval.size();i++)
        for (unsigned long j=0;j<OMR.columns();j++)
            retval[i] += OMR.entry(i,j)*temp[j];
    // the above takes temp and multiplies it by the matrix OMR.

    return retval;
}


/*
 * The marked abelian group was defined by matrices M and N
 * with M*N==0.  Think of M as m by l and N as l by n.
 * When the group was initialized, it was computed to be isomorphic
 * to some Z^d + Z_{d1} + ... + Z_{dk} where d1 | d2 | ... | dk
 * this routine assumes element is in Z^l, and it returns a vector
 * of length d+k where the first d elements represent which class the
 * vector projects to in Z^d, and the last k elements represent the
 * projections to Z_{d1} + ... + Z_{dk}. Of these last elements, they
 * will be returned mod di respectively. Returns an empty vector if
 * element is not in the kernel of M. element is assumed to have
 * OM.columns()==ON.rows() entries.
 */
std::vector<NLargeInteger> NMarkedAbelianGroup::getSNFisoRep(
        std::vector<NLargeInteger>& element)  const {
    std::vector<NLargeInteger> retval(snfrank+InvFacList.size(),
        NLargeInteger::zero);
    // apply OMRi, crop, then apply ornC, tidy up and return.
    std::vector<NLargeInteger> nullvec(0); // this is returned if element
                                           // not in ker(M)

    std::vector<NLargeInteger> temp(ON.rows(),
        NLargeInteger::zero); // this holds OMRi * element
    // if first rankOM entries are 0, then element is in the kernel

    bool eltinker=true;

    // set up temp.
    for (unsigned long i=0;i<ON.rows();i++)
        for (unsigned long j=0;j<ON.rows();i++)
            temp[i] += OMRi.entry(i,j)*element[j];
    for (unsigned long i=0;i<rankOM;i++)
        if (temp[i] != 0)
            eltinker=false;
    // ON.rows - rankOM == ORN.rows()
    if (eltinker==true) {
        // set up retval. The first snfrank elts are the free generators
        for (unsigned long i=0;i<snfrank;i++)
            for (unsigned long j=rankOM;j<ON.rows();i++)
                retval[i] += ornC.entry(snffreeindex+i,j)*temp[j];
        // the remaining InvFacList.size() elts are torsion generators.
        for (unsigned long i=0;i<ifNum;i++) {
            for (unsigned long j=rankOM;j<ON.rows();i++)
                retval[i+snfrank] += ornC.entry(ifLoc+i,j)*temp[j];
            retval[i+snfrank] = (retval[i+snfrank] % InvFacList[i]);
        }
    } else retval=nullvec;


    return retval;
}


NHomMarkedAbelianGroup::NHomMarkedAbelianGroup(const NMarkedAbelianGroup& dom,
        const NMarkedAbelianGroup& ran,
        const NMatrixInt &mat):
        domain(dom), range(ran), matrix(mat),
        reducedMatrix(0), kernel(0), coKernel(0), image(0),
        reducedKernelLattice(0) {
}


NHomMarkedAbelianGroup::NHomMarkedAbelianGroup(const NHomMarkedAbelianGroup& g):
        ShareableObject(), domain(g.domain), range(g.range), matrix(g.matrix) {
    if (g.reducedMatrix) {
        reducedMatrix = new NMatrixInt(*g.reducedMatrix);
    } else reducedMatrix = 0;

    if (g.kernel) {
        kernel = new NMarkedAbelianGroup(*g.kernel);
    } else kernel = 0;

    if (g.coKernel) {
        coKernel = new NMarkedAbelianGroup(*g.coKernel);
    } else coKernel = 0;

    if (g.image) {
        image = new NMarkedAbelianGroup(*g.image);
    } else image = 0;

    if (g.reducedKernelLattice) {
        reducedKernelLattice = new NMatrixInt(*g.reducedKernelLattice);
    } else reducedKernelLattice = 0;
}

void NHomMarkedAbelianGroup::computeReducedMatrix() {
    if (!reducedMatrix) {
        unsigned long i,j,k;

        NMatrixInt kerMatrix( matrix.rows()-range.getRankOM(),
                matrix.columns()-domain.getRankOM() );
        // kerMatrix = truncate (range.getMRBi() * matrix * domain.getMRBi)
        // to construct this we do it in two steps:
        // step 1) temp1 = truncate columns (matrix * domain.getMRBi )
        // step 2) kerMatrix = truncate rows (range.getMRBi * temp1 )

        NMatrixInt dcckb(domain.getMRB());
        NMatrixInt rcckb(range.getMRBi());

        NMatrixInt temp1( matrix.rows(), matrix.columns()-domain.getRankOM() );
        for (i=0;i<temp1.rows();i++)
            for (j=0;j<temp1.columns();j++)
                for (k=0;k<matrix.columns();k++)
                    temp1.entry(i,j) += matrix.entry(i,k) *
                        dcckb.entry(k,j + domain.getRankOM() );

        for (i=0;i<kerMatrix.rows();i++)
            for (j=0;j<kerMatrix.columns();j++)
                for (k=0;k<rcckb.rows();k++)
                    kerMatrix.entry(i,j) +=
                        rcckb.entry(i+range.getRankOM(), k) * temp1.entry(k,j);

        NMatrixInt redMatrix( kerMatrix.rows()-range.getTorLoc(),
                kerMatrix.columns()-domain.getTorLoc() );

        NMatrixInt dccqb(domain.getNCBi());
        NMatrixInt rccqb(range.getNCB());

        NMatrixInt temp2( kerMatrix.rows(),
            kerMatrix.columns() - domain.getTorLoc() );
        for (i=0;i<temp2.rows();i++)
            for (j=0;j<temp2.columns();j++)
                for (k=0;k<kerMatrix.columns();k++) {
                    temp2.entry(i,j) += kerMatrix.entry(i,k) *
                        dccqb.entry(k,j + domain.getTorLoc() );
                }

        for (i=0;i<redMatrix.rows();i++)
            for (j=0;j<redMatrix.columns();j++)
                for (k=0;k<rccqb.rows();k++)
                    redMatrix.entry(i,j) +=
                        rccqb.entry(i+range.getTorLoc(), k) * temp2.entry(k,j);

        reducedMatrix = new NMatrixInt(redMatrix);
    }
}

void NHomMarkedAbelianGroup::computeReducedKernelLattice() {
    if (!reducedKernelLattice) {
        computeReducedMatrix();

        unsigned long i;

        NMatrixInt redMatrix(*reducedMatrix);

        // the kernel is the dcLpreimage lattice mod the domain lattice.
        // so after computing the dcLpreimage lattice, we need to represent
        // the domain lattice in its coordinates.

        std::vector<NLargeInteger> dcL(range.getRank() +
            range.getNumberOfInvariantFactors() );
        for (i=0; i<dcL.size(); i++)
            if (i<range.getNumberOfInvariantFactors())
                dcL[i]=range.getInvariantFactor(i);
            else
                dcL[i]="0";

        reducedKernelLattice = new NMatrixInt( RBADD_preImageOfLattice(
            redMatrix, dcL ) );
    }
}

void NHomMarkedAbelianGroup::computeKernel() {
    if (!kernel) {
        computeReducedKernelLattice();
        NMatrixInt dcLpreimage( *reducedKernelLattice );

        unsigned long i,j,k;

        NMatrixInt R( dcLpreimage.columns(), dcLpreimage.columns() );
        NMatrixInt Ri( dcLpreimage.columns(), dcLpreimage.columns() );
        NMatrixInt C( dcLpreimage.rows(), dcLpreimage.rows() );
        NMatrixInt Ci( dcLpreimage.rows(), dcLpreimage.rows() );

        RBMOD_smithNormalForm( dcLpreimage, R, Ri, C, Ci );

        // the matrix representing the domain lattice in dcLpreimage
        // coordinates is given by domainLattice * R * (dcLpreimage inverse) * C

        NMatrixInt workMat( dcLpreimage.columns(),
            domain.getNumberOfInvariantFactors() );

        for (i=0;i<workMat.rows();i++)
            for (j=0;j<workMat.columns();j++)
                for (k=0;k<R.columns();k++) {
                    workMat.entry(i,j) += (domain.getInvariantFactor(j) *
                        R.entry(i,k) * C.entry(k,j) ) / dcLpreimage.entry(k,k);
                }

        NMatrixInt dummy( 1, dcLpreimage.columns() );

        NMarkedAbelianGroup retval( dummy, workMat );

        kernel = new NMarkedAbelianGroup(retval);
    }
}



void NHomMarkedAbelianGroup::computeCoKernel() {
    if (!coKernel) {
        computeReducedMatrix();

        NMatrixInt ccrelators( reducedMatrix->rows(),
            reducedMatrix->columns() + range.getNumberOfInvariantFactors() );
        unsigned i,j;
        for (i=0;i<reducedMatrix->rows();i++)
            for (j=0;j<reducedMatrix->columns();j++)
                ccrelators.entry(i,j)=reducedMatrix->entry(i,j);
        for (i=0;i<range.getNumberOfInvariantFactors();i++)
            ccrelators.entry(i,i+reducedMatrix->columns())=
                range.getInvariantFactor(i);

        NMatrixInt ccgenerators( 1, reducedMatrix->rows() );
        NMarkedAbelianGroup retval(ccgenerators, ccrelators);

        coKernel = new NMarkedAbelianGroup(retval);
    }
}




void NHomMarkedAbelianGroup::computeImage() {
    if (!image) {
        computeReducedKernelLattice();
        NMatrixInt dcLpreimage( *reducedKernelLattice );

        unsigned long i,j;

        NMatrixInt imgCCm(1, dcLpreimage.rows() );
        NMatrixInt imgCCn(dcLpreimage.rows(),
            dcLpreimage.columns() + domain.getNumberOfInvariantFactors() );

        for (i=0;i<domain.getNumberOfInvariantFactors();i++)
            imgCCn.entry(i,i) = domain.getInvariantFactor(i);

        for (i=0;i<imgCCn.rows();i++)
            for (j=0;j< dcLpreimage.columns(); j++)
                imgCCn.entry(i,j+domain.getNumberOfInvariantFactors()) =
                    dcLpreimage.entry(i,j);

        NMarkedAbelianGroup retval(imgCCm, imgCCn);

        image = new NMarkedAbelianGroup(retval);
    }
}



bool NHomMarkedAbelianGroup::isEpic() const {
    return getCoKernel().isTrivial();
}

bool NHomMarkedAbelianGroup::isMonic() const {
    return getKernel().isTrivial();
}

bool NHomMarkedAbelianGroup::isIso() const {
    return (getCoKernel().isTrivial() && getKernel().isTrivial());
}

bool NHomMarkedAbelianGroup::isZero() const {
    return getImage().isTrivial();
}



NMarkedAbelianGroup NHomMarkedAbelianGroup::getKernel() const {
    // Cast away const to compute the kernel -- the only reason we're
    // changing data members now is because we delayed calculations
    // until they were really required.
    const_cast<NHomMarkedAbelianGroup*>(this)->computeKernel();
    return *kernel;
}

NMarkedAbelianGroup NHomMarkedAbelianGroup::getImage() const {
    // Cast away const to compute the kernel -- the only reason we're
    // changing data members now is because we delayed calculations
    // until they were really required.
    const_cast<NHomMarkedAbelianGroup*>(this)->computeImage();
    return *image;
}

NMarkedAbelianGroup NHomMarkedAbelianGroup::getCoKernel() const {
    // Cast away const to compute the kernel -- the only reason we're
    // changing data members now is because we delayed calculations
    // until they were really required.
    const_cast<NHomMarkedAbelianGroup*>(this)->computeCoKernel();
    return *coKernel;
}



void NHomMarkedAbelianGroup::writeRedMatrix(std::ostream& out) const {
    // Cast away const to compute the reduced matrix -- the only reason we're
    // changing data members now is because we delayed calculations
    // until they were really required.
    const_cast<NHomMarkedAbelianGroup*>(this)->computeReducedMatrix();

    unsigned long i,j;
    out<<"Reduced Matrix is "<<reducedMatrix->rows()<<" by "
        <<reducedMatrix->columns()<<" corresponding to domain ";
    domain.writeTextShort(out);
    out<<" and range ";
    range.writeTextShort(out);
    out<<"\n";
    for (i=0;i<reducedMatrix->rows();i++) {
        out<<"[";
        for (j=0;j<reducedMatrix->columns();j++) {
            out<<reducedMatrix->entry(i,j);
            if (j+1 < reducedMatrix->columns()) out<<" ";
        }
        out<<"]\n";
    }
}

void NHomMarkedAbelianGroup::writeTextShort(std::ostream& out) const {
    if (isIso())
        out<<"isomorphism";
    else if (isZero())
        out<<"zero map";
    else if (isMonic()) { // monic not epic
        out<<"monic, with cokernel ";
        getCoKernel().writeTextShort(out);
    } else if (isEpic()) { // epic not monic
        out<<"epic, with kernel ";
        getKernel().writeTextShort(out);
    } else { // nontrivial not epic, not monic
        out<<"kernel ";
        getKernel().writeTextShort(out);
        out<<" | cokernel ";
        getCoKernel().writeTextShort(out);
        out<<" | image ";
        getImage().writeTextShort(out);
    }
}

} // namespace regina


