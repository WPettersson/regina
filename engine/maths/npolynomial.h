
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#ifndef __NPOLYNOMIAL_H
#ifndef __DOXYGEN
#define __NPOLYNOMIAL_H
#endif

/*! \file maths/npolynomial.h
 *  \brief Implements single variable polynomials over arbitrary rings.
 */

#include "regina-core.h"
#include <iostream>

namespace regina {

/**
 * \weakgroup maths
 * @{
 */

/**
 * Represents a single-variable polynomial with coefficients of type \a T.
 * All exponents in the polynomial must be non-negative (so you can
 * represent <tt>2+3x</tt> but not <tt>1+1/x</tt>).
 *
 * The type \a T must represent a ring with no zero divisors.
 * In particular, it must:
 *
 * - support basic arithmetic operations;
 * - support assignments of the form <tt>x = int</tt> and
 *   tests of the form <tt>x == int</tt> and <tt>x < int</tt>;
 * - have a default constructor that assigns an explicit value of zero.
 *
 * This means that Regina's numerical types such as NInteger and NRational
 * are supported, but native data types such as int and long are not
 * (since they have no zero-initialising default constructor).
 *
 * The underlying storage method for this class is dense (i.e., all
 * coefficients are explicitly stored, including zero coefficients).
 *
 * \ifacespython In Python, the class NPolynomial refers to the specific
 * template class NPolynomial<NRational>.
 */
template <typename T>
class NPolynomial {
    private:
        size_t degree_;
            /**< The degree of the polynomial.  Here the zero polynomial
                 is considered to have degree zero. */
        T* coeff_;
            /**< An array of length (\a degree_ + 1), storing the
                 coefficients of the polynomial.  Specifically,
                 coeff_[i] stores the coefficient of <tt>x^i</tt>. */

    public:
        /**
         * Creates the zero polynomial.
         */
        NPolynomial();

        /**
         * Creates the polynomial <tt>x^d</tt> for the given degree \a d.
         *
         * @param degree the degree of the new polynomial.
         */
        explicit NPolynomial(size_t degree);

        /**
         * Creates a new copy of the given polynomial.
         *
         * A note for developers: even though this routine is identical to
         * the templated copy constructor, it must be declared and
         * implemented separately.  Otherwise the compiler might create
         * its own (incorrect) copy constructor automatically.
         *
         * @param value the polynomial to clone.
         */
        NPolynomial(const NPolynomial<T>& value);

        /**
         * Creates a new copy of the given polynomial.
         *
         * \pre Objects of type \a T can be assigned values of type \a U.
         *
         * @param value the polynomial to clone.
         */
        template <typename U>
        NPolynomial(const NPolynomial<U>& value);

        /**
         * Creates a new polynomial from the given sequence of coefficients.
         * The coefficients should appear in order from the constant
         * coefficient to the leading coefficient.
         *
         * There is no problem if the leading coefficient (i.e., the
         * last coefficient in the sequence) is zero.
         * An empty sequence will be treated as the zero polynomial.
         *
         * \pre Objects of type \a T can be assigned values from
         * dereferenced iterators of type \a iterator.
         *
         * \ifacespython Instead of a pair of iterators, this routine
         * takes a python list of coefficients.
         *
         * @param begin the beginning of the sequence of coefficients.
         * @param end a past-the-end iterator indicating the end of the
         * sequence of coefficients.
         */
        template <typename iterator>
        NPolynomial(iterator begin, iterator end);

        /**
         * Destroys this polynomial.
         */
        ~NPolynomial();

        /**
         * Sets this to become the zero polynomial.
         */
        void init();

        /**
         * Sets this to become the polynomial <tt>x^d</tt> for the
         * given degree \a d.
         *
         * @param degree the new degree of this polynomial.
         */
        void init(size_t degree);

        /**
         * Sets this to become the polynomial described by the given
         * sequence of coefficients.
         * The coefficients should appear in order from the constant
         * coefficient to the leading coefficient.
         *
         * There is no problem if the leading coefficient (i.e., the
         * last coefficient in the sequence) is zero.
         * An empty sequence will be treated as the zero polynomial.
         *
         * \pre Objects of type \a T can be assigned values from
         * dereferenced iterators of type \a iterator.
         *
         * \ifacespython Instead of a pair of iterators, this routine
         * takes a python list of coefficients.
         *
         * @param begin the beginning of the sequence of coefficients.
         * @param end a past-the-end iterator indicating the end of the
         * sequence of coefficients.
         */
        template <typename iterator>
        void init(iterator begin, iterator end);

        /**
         * Returns the degree of this polynomial.
         * This is the largest exponent with a non-zero coefficient.
         *
         * For the purposes of this class, the zero polynomial is
         * considered to have degree zero.
         *
         * @return the degree of this polynomial.
         */
        size_t degree() const;

        /**
         * Returns whether this is the zero polynomial.
         *
         * @return \c true if and only if this is the zero polynomial.
         */
        bool isZero() const;

        /**
         * Returns whether this polynomial is monic.  A \e monic
         * polynomial is a non-zero polynomial whose leading coefficient
         * is one.
         *
         * @return \c true if and only if this is monic.
         */
        bool isMonic() const;

        /**
         * Returns the leading coefficient of this polynomial.
         * If this is the zero polynomial, then the leading coefficient
         * will be zero.
         *
         * @return the leading coefficient of this polynomial.
         */
        const T& leading() const;

        /**
         * Returns the given coefficient of this polynomial.
         *
         * @param exp the exponent of the term whose coefficient should
         * be returned.  This must be between 0 and degree() inclusive.
         * @return the coefficient of the given term.
         */
        const T& operator [] (size_t exp) const;

        /**
         * Changes the given coefficient of this polynomial.
         *
         * It is fine to set the leading coefficient to zero, though note
         * that degree() will now return a smaller value as a result.
         *
         * It is also fine to set a coefficient whose exponent is larger
         * than the current degree; this time degree() will now return a
         * larger value (unless the given coefficient is zero).  Such an
         * operation is expensive, however, since it will require
         * deallocating and reallocating the full list of coefficients.
         *
         * \ifacespython This set() routine is available, but you can
         * also set coefficients directly using syntax of the form
         * <tt>p[exp] = value</tt>.
         *
         * @param exp the exponent of the term whose coefficient should
         * be changed.
         * @param value the new value of this coefficient.
         */
        void set(size_t exp, const T& value);

        /**
         * Tests whether this and the given polynomial are equal.
         *
         * @param rhs the polynomial to compare with this.
         * @return \c true if and only if this and the given polynomial
         * are equal.
         */
        bool operator == (const NPolynomial<T>& rhs) const;

        /**
         * Tests whether this and the given polynomial are not equal.
         *
         * @param rhs the polynomial to compare with this.
         * @return \c true if and only if this and the given polynomial
         * are not equal.
         */
        bool operator != (const NPolynomial<T>& rhs) const;

        /**
         * Sets this to be a copy of the given polynomial.
         *
         * This and the given polynomial need not have the same degree
         * (but if they do not, then the degree of this polynomial will
         * of course change).
         *
         * A note to developers: although this is identical to the templated
         * assignment operator, it must be declared and implemented separately.
         * See the copy constructor for further details.
         *
         * @param value the polynomial to copy.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator = (const NPolynomial<T>& value);

        /**
         * Sets this to be a copy of the given polynomial.
         *
         * This and the given polynomial need not have the same degree
         * (but if they do not, then the degree of this polynomial will
         * of course change).
         *
         * @param value the polynomial to copy.
         * @return a reference to this polynomial.
         */
        template <typename U>
        NPolynomial& operator = (const NPolynomial<U>& value);

        /**
         * Swaps the contents of this and the given polynomial.
         * This is a fast (constant time) operation.
         *
         * This and the given polynomial do not need to have the same degree.
         *
         * @param other the polynomial whose contents should be swapped
         * with this.
         */
        void swap(NPolynomial<T>& other);

        /**
         * Multiplies this polynomial by the given constant.
         *
         * @param scalar the scalar factor to multiply by.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator *= (const T& scalar);

        /**
         * Divides this polynomial by the given constant.
         *
         * This uses the division operator /= for the coefficient type \a T.
         *
         * \pre The argument \a scalar is non-zero.
         *
         * @param scalar the scalar factor to divide by.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator /= (const T& scalar);

        /**
         * Adds the given polynomial to this.
         *
         * The given polynomial need not have the same degree as this.
         * Note that the degree of this polynomial might change as a
         * result of this operation.
         *
         * @param other the polynomial to add to this.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator += (const NPolynomial<T>& other);

        /**
         * Subtracts the given polynomial from this.
         *
         * The given polynomial need not have the same degree as this.
         * Note that the degree of this polynomial might change as a
         * result of this operation.
         *
         * @param other the polynomial to subtract from this.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator -= (const NPolynomial<T>& other);

        /**
         * Multiplies this by the given polynomial.
         *
         * @param other the polynomial to multiply this by.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator *= (const NPolynomial<T>& other);

        /**
         * Divides this by the given polynomial.
         *
         * More precisely: suppose there exist polynomials \a q and \a r with
         * coefficients of type \a T for which <tt>this = q.other + r</tt>,
         * and where \a r has smaller degree than \a other.  Then we call
         * \a q the \e quotient, and \a r the \e remainder.
         *
         * This routine replaces this polynomial with the quotient \a q,
         * and discards the remainder.  If you need to keep the remainder
         * also, then call divisionAlg() instead.
         *
         * Coefficients are divided using the operator /= on type \a T.
         *
         * If your coefficient type \a T is not a field (e.g., if \a T
         * is NInteger), you must be sure to know in advance that the
         * quotient exists (see the precondition below).  Otherwise the
         * behaviour of this routine is undefined.
         *
         * \pre The given polynomial is not the zero polynomial.
         *
         * \pre The quotient as defined above exists.  If \a T is a field
         * type (e.g., if \a T is NRational) then this is true automatically.
         * If not (e.g., if \a T is NInteger) then this requires some
         * prior knowledge about the arguments.
         *
         * @param other the polynomial to divide this by.
         * @return a reference to this polynomial.
         */
        NPolynomial& operator /= (const NPolynomial<T>& other);

        /**
         * Divides this by the given divisor, and extracts both the
         * quotient and the remainder.
         *
         * More precisely: suppose there exist polynomials \a q and \a r with
         * coefficients of type \a T for which <tt>this = q.divisor + r</tt>,
         * and where \a r has smaller degree than \a divisor.  Then this
         * routine sets the given polynomial \a quotient to \a q, and sets
         * the given polynomial \a remainder to \a r.
         *
         * If you do not need the remainder (e.g., if you know in
         * advance that \a divisor divides into this polynomial exactly),
         * then you can use the division operator /= instead, which will
         * be a little faster.
         *
         * If your coefficient type \a T is not a field (e.g., if \a T
         * is NInteger), you must be sure to know in advance that the
         * quotient exists (see the precondition below).  Otherwise the
         * behaviour of this routine is undefined.
         *
         * Coefficients are divided using the operator /= on type \a T.
         *
         * \pre The given divisor is not the zero polynomial.
         *
         * \pre The quotient as defined above exists.  If \a T is a field
         * type (e.g., if \a T is NRational) then this is true automatically.
         * If not (e.g., if \a T is NInteger) then this requires some
         * prior knowledge about the arguments.
         *
         * \pre Neither \a quotient nor \a remainder is a reference to
         * this polynomial.
         *
         * \ifacespython The arguments \a quotient and \a remainder are
         * missing; instead these are passed back through the return
         * value of the function.  Specifically, this function returns a
         * (\a quotient, \a remainder) pair.
         *
         * @param divisor the polynomial to divide by this.
         * @param quotient a polynomial whose contents will be destroyed and
         * replaced with the quotient \a q, as described above.
         * @param remainder a polynomial whose contents will be destroyed
         * and replaced with the remainder \a r, as described above.
         */
        void divisionAlg(const NPolynomial<T>& divisor,
            NPolynomial<T>& quotient, NPolynomial<T>& remainder) const;

        /**
         * Calculates the greatest common divisor of this and the given
         * polynomial, and finds a linear combination of these
         * polynomials that gives this gcd.
         *
         * The greatest common divisor will be a monic polynomial.
         * The polynomials returned in \a u and \a v will satisfy
         * <tt>u*this + v*other = gcd</tt>.
         *
         * As a special case, gcd(0,0) is considered to be zero.
         *
         * \pre The coefficient type \a T represents a field.  In particular,
         * NRational is supported but NInteger is not.
         *
         * @param other the polynomial whose greatest common divisor with this
         * polynomial we should compute.
         * @param gcd a polynomial whose contents will be destroyed and
         * replaced with the greatest common divisor \a d, as described above.
         * @param u a polynomial whose contents will be destroyed and
         * replaced with \a u, as described above.
         * @param v a polynomial whose contents will be destroyed and
         * replaced with \a v, as described above.
         */
        template <typename U>
        void gcdWithCoeffs(const NPolynomial<U>& other,
            NPolynomial<T>& gcd, NPolynomial<T>& u, NPolynomial<T>& v) const;
};

/**
 * Writes the given polynomial to the given output stream in
 * human-readable form.
 *
 * @param out the output stream to which to write.
 * @param p the polynomial to write.
 * @return a reference to the given output stream.
 */
template <typename T>
std::ostream& operator << (std::ostream& out, const NPolynomial<T>& p);

/*@}*/

template <typename T>
inline NPolynomial<T>::NPolynomial() : degree_(0), coeff_(new T[1]) {
    // The default constructor for T already initialises coeff_[0] to zero.
}

template <typename T>
inline NPolynomial<T>::NPolynomial(size_t degree) :
        degree_(degree), coeff_(new T[degree + 1]) {
    coeff_[degree] = 1;
}

template <typename T>
template <typename iterator>
inline NPolynomial<T>::NPolynomial(iterator begin, iterator end) : coeff_(0) {
    init(begin, end);
}

template <typename T>
inline NPolynomial<T>::NPolynomial(const NPolynomial<T>& value) :
        degree_(value.degree()), coeff_(new T[value.degree() + 1]) {
    for (size_t i = 0; i <= degree_; ++i)
        coeff_[i] = value[i];
}

template <typename T>
template <typename U>
inline NPolynomial<T>::NPolynomial(const NPolynomial<U>& value) :
        degree_(value.degree()), coeff_(new T[value.degree() + 1]) {
    for (size_t i = 0; i <= degree_; ++i)
        coeff_[i] = value[i];
}

template <typename T>
inline NPolynomial<T>::~NPolynomial() {
    delete[] coeff_;
}

template <typename T>
inline void NPolynomial<T>::init() {
    delete[] coeff_;
    degree_ = 0;
    coeff_ = new T[1];
    // coeff_[0] is initialised to 0 automatically.
}

template <typename T>
inline void NPolynomial<T>::init(size_t degree) {
    delete[] coeff_;
    degree_ = degree;
    coeff_ = new T[degree + 1];
    coeff_[degree] = 1;
}

template <typename T>
template <typename iterator>
void NPolynomial<T>::init(iterator begin, iterator end) {
    delete[] coeff_;

    if (begin == end) {
        degree_ = 0;
        coeff_ = new T[1];
        return;
    }

    degree_ = end - begin - 1;
    coeff_ = new T[degree_ + 1];

    size_t i = 0;
    while (begin != end)
        coeff_[i++] = *begin++;

    // The leading coefficient might be zero.
    while (degree_ > 0 && coeff_[degree_] == 0)
        --degree_;
}

template <typename T>
inline size_t NPolynomial<T>::degree() const {
    return degree_;
}

template <typename T>
inline bool NPolynomial<T>::isZero() const {
    return (degree_ == 0 && coeff_[0] == 0);
}

template <typename T>
inline bool NPolynomial<T>::isMonic() const {
    return (coeff_[degree_] == 1);
}

template <typename T>
inline const T& NPolynomial<T>::leading() const {
    return coeff_[degree_];
}

template <typename T>
inline const T& NPolynomial<T>::operator [] (size_t exp) const {
    return coeff_[exp];
}

template <typename T>
void NPolynomial<T>::set(size_t exp, const T& value) {
    if (exp < degree_) {
        coeff_[exp] = value;
    } else if (exp == degree_) {
        if (value == 0) {
            --degree_;
            while (degree_ > 0 && coeff_[degree_] == 0)
                --degree_;
        } else {
            coeff_[exp] = value;
        }
    } else if (! (value == 0)) {
        // The degree will increase.
        T* newCoeff = new T[exp + 1];
        size_t i;
        for (i = 0; i <= degree_; ++i)
            newCoeff[i] = coeff_[i];
        newCoeff[exp] = value;

        delete[] coeff_;
        coeff_ = newCoeff;
        degree_ = exp;
    }
}

template <typename T>
inline bool NPolynomial<T>::operator == (const NPolynomial<T>& rhs) const {
    if (degree_ != rhs.degree_)
        return false;
    for (size_t i = 0; i <= degree_; ++i)
        if (coeff_[i] != rhs.coeff_[i])
            return false;
    return true;
}

template <typename T>
inline bool NPolynomial<T>::operator != (const NPolynomial<T>& rhs) const {
    if (degree_ != rhs.degree_)
        return true;
    for (size_t i = 0; i <= degree_; ++i)
        if (coeff_[i] != rhs.coeff_[i])
            return true;
    return false;
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator = (const NPolynomial<T>& other) {
    // This works even if &other == this, since we don't reallocate if
    // the degrees are equal.
    if (degree_ < other.degree_) {
        delete[] coeff_;
        coeff_ = new T[other.degree_ + 1];
    }
    degree_ = other.degree_;
    for (size_t i = 0; i <= degree_; ++i)
        coeff_[i] = other.coeff_[i];
    return *this;
}

template <typename T>
template <typename U>
NPolynomial<T>& NPolynomial<T>::operator = (const NPolynomial<U>& other) {
    // This works even if &other == this, since we don't reallocate if
    // the degrees are equal.
    if (degree_ < other.degree()) {
        delete[] coeff_;
        coeff_ = new T[other.degree() + 1];
    }
    degree_ = other.degree();
    for (size_t i = 0; i <= degree_; ++i)
        coeff_[i] = other[i];
    return *this;
}

template <typename T>
inline void NPolynomial<T>::swap(NPolynomial<T>& other) {
    std::swap(degree_, other.degree_);
    std::swap(coeff_, other.coeff_);
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator *= (const T& scalar) {
    if (scalar == 0)
        init();
    else {
        for (size_t i = 0; i <= degree_; ++i)
            coeff_[i] *= scalar;
    }
    return *this;
}

template <typename T>
inline NPolynomial<T>& NPolynomial<T>::operator /= (const T& scalar) {
    for (size_t i = 0; i <= degree_; ++i)
        coeff_[i] /= scalar;
    return *this;
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator += (const NPolynomial<T>& other) {
    // This works even if &other == this, since we don't reallocate if
    // the degrees are equal.
    size_t i;
    if (degree_ < other.degree_) {
        T* copy = new T[other.degree_ + 1];
        for (i = 0; i <= degree_; ++i)
            copy[i] = coeff_[i];
        delete[] coeff_;
        coeff_ = copy;
        degree_ = other.degree_;
    }
    for (i = 0; i <= other.degree_; ++i)
        coeff_[i] += other.coeff_[i];

    // We might have zeroed out the leading coefficient.
    while (degree_ > 0 && coeff_[degree_] == 0)
        --degree_;

    return *this;
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator -= (const NPolynomial<T>& other) {
    // This works even if &other == this, since we don't reallocate if
    // the degrees are equal.
    size_t i;
    if (degree_ < other.degree_) {
        T* copy = new T[other.degree_ + 1];
        for (i = 0; i <= degree_; ++i)
            copy[i] = coeff_[i];
        delete[] coeff_;
        coeff_ = copy;
        degree_ = other.degree_;
    }
    for (i = 0; i <= other.degree_; ++i)
        coeff_[i] -= other.coeff_[i];

    // We might have zeroed out the leading coefficient.
    while (degree_ > 0 && coeff_[degree_] == 0)
        --degree_;

    return *this;
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator *= (const NPolynomial<T>& other) {
    if (isZero())
        return *this;
    if (other.isZero()) {
        init();
        return *this;
    }

    // The following code works even if &other == this, since we construct the
    // coefficients of the product in a separate section of memory.
    size_t i, j;
    T* ans = new T[degree_ + other.degree_ + 1];
    for (i = 0; i <= degree_; ++i)
        for (j = 0; j <= other.degree_; ++j)
            ans[i + j] += (coeff_[i] * other.coeff_[j]);

    delete[] coeff_;
    coeff_ = ans;
    degree_ += other.degree_;

    // Both leading coefficients are non-zero, so the degree is correct.
    return *this;
}

template <typename T>
NPolynomial<T>& NPolynomial<T>::operator /= (const NPolynomial<T>& other) {
    // The code below breaks if other and *this are the same object, so
    // treat this case specially.
    if (&other == this) {
        init(0);
        return *this;
    }

    size_t i, j;
    if (other.degree_ == 0) {
        for (i = 0; i <= degree_; ++i)
            coeff_[i] /= other.coeff_[0];
        return *this;
    }

    // The divisor has positive degree.
    if (degree_ == 0)
        return *this;

    // Both this and the divisor have positive degree.
    T* remainder = coeff_;
    coeff_ = new T[degree_ - other.degree_ + 1];
    for (i = degree_; i >= other.degree_; --i) {
        coeff_[i - other.degree_] =
            remainder[i] /= other.coeff_[other.degree_];
        for (j = 0; j <= other.degree_; ++j)
            remainder[j + i - other.degree_] -=
                (coeff_[i - other.degree_] * other.coeff_[j]);
    }
    delete[] remainder;
    degree_ -= other.degree_;

    return *this;
}

template <typename T>
void NPolynomial<T>::divisionAlg(const NPolynomial<T>& divisor,
        NPolynomial<T>& quotient, NPolynomial<T>& remainder) const {
    // The code below breaks if divisor and *this are the same object, so
    // treat this case specially.
    if (&divisor == this) {
        quotient.init(0);
        remainder.init();
        return;
    }

    if (divisor.degree_ > degree_) {
        quotient.init();
        remainder = *this;
        return;
    }

    size_t i, j;
    if (divisor.degree_ == 0) {
        quotient = *this;
        for (i = 0; i <= quotient.degree_; ++i)
            quotient.coeff_[i] /= divisor.coeff_[0];
        remainder.init();
        return;
    }

    // From here we have: 0 < deg(divisor) <= deg(this).
    // In particular, both this and divisor have strictly positive degree.

    quotient.degree_ = degree_ - divisor.degree_;
    delete[] quotient.coeff_;
    quotient.coeff_ = new T[quotient.degree_ + 1];

    remainder = *this;

    for (i = degree_; i >= divisor.degree_; --i) {
        quotient.coeff_[i - divisor.degree_] = remainder.coeff_[i];
        quotient.coeff_[i - divisor.degree_] /= divisor.coeff_[divisor.degree_];
        for (j = 0; j <= divisor.degree_; ++j)
            remainder.coeff_[j + i - divisor.degree_] -=
                (quotient.coeff_[i - divisor.degree_] * divisor.coeff_[j]);
    }

    // Although the degree of the quotient is correct, the remainder
    // might have zero coefficients at any (or all) positions.
    remainder.degree_ = divisor.degree_ - 1;
    while (remainder.degree_ > 0 && remainder.coeff_[remainder.degree_] == 0)
        --remainder.degree_;
}

template <typename T>
template <typename U>
void NPolynomial<T>::gcdWithCoeffs(const NPolynomial<U>& other,
        NPolynomial<T>& gcd, NPolynomial<T>& u, NPolynomial<T>& v) const {
    // Special-case situations where one or both polynomials are zero.
    if (other.isZero()) {
        if (isZero()) {
            // gcd(0, 0) = 0.
            gcd.init();
            u.init();
            v.init();
            return;
        }

        // gcd(this, 0) = this / this.leading()
        gcd = *this;
        gcd /= coeff_[degree_];
        u.init(0);
        u.coeff_[0] /= coeff_[degree_];
        v.init();
        return;
    }
    if (isZero()) {
        // gcd(0, other) = other / other.leading()
        gcd = other;
        gcd /= other[other.degree()];
        u.init();
        v.init(0);
        v.coeff_[0] /= other[other.degree()];
        return;
    }

    // We use Euclid's algorithm to find gcd(this, other).
    //
    // At each stage we maintain the invariants:
    //
    //   u * this + v * other = x
    //   uu * this + vv * other = y
    //   deg(x) >= deg(y)
    //
    // We begin with (x, y, u, v, uu, vv) = (this, other, 1, 0, 0, 1).
    // The iteration step, assuming x = q * y + r, is then:
    //
    //   (x, y, u, v, uu, vv) -> (y, r, uu, vv, u-q*uu, v-q*vv)
    //
    // We finish with (x, y) = (gcd, 0).
    //
    // In the code below we use the given polyomial gcd to store x throughout
    // the algorithm.

    gcd = *this;
    NPolynomial<T> y(other);
    u.init(0);
    v.init();
    NPolynomial<T> uu;
    NPolynomial<T> vv(0);

    if (gcd.degree() < y.degree()) {
        gcd.swap(y);
        u.swap(uu);
        v.swap(vv);
    }

    NPolynomial<T> tmp, q, r;
    while (! y.isZero()) {
        gcd.divisionAlg(y, q, r);

        tmp = q;
        tmp *= uu;
        u -= tmp;

        tmp = q;
        tmp *= vv;
        v -= tmp;

        u.swap(uu);
        v.swap(vv);
        gcd.swap(y);
        y.swap(r);
    }

    // Make the gcd monic.
    // Note that gcd != 0 at this point, so this is indeed possible.
    if (! gcd.isMonic()) {
        T leading(gcd[gcd.degree_]);
        gcd /= leading;
        u /= leading;
        v /= leading;
    }
}

template <typename T>
std::ostream& operator << (std::ostream& out, const NPolynomial<T>& p) {
    if (p.degree() == 0)
        return out << p[0];

    size_t i = p.degree();
    while (true) {
        if (i == 0) {
            // This is the constant term.
            // If non-zero, it is not the first term being output.
            if (p[i] < 0)
                out << " - " << (-p[i]);
            else if (p[i] != 0)
                out << " + " << p[i];
            break;
        }
        if (p[i] != 0) {
            // This is non-zero, and not the constant term.
            if (i == p.degree()) {
                // This is the first term being output.
                if (p[i] == -1)
                    out << "- ";
                else if (p[i] != 1)
                    out << p[i] << ' ';
            } else {
                // This is not the first term being output.
                if (p[i] == -1)
                    out << " - ";
                else if (p[i] < 0)
                    out << " - " << (-p[i]) << ' ';
                else if (p[i] == 1)
                    out << " + ";
                else
                    out << " + " << p[i] << ' ';
            }
            if (i == 1)
                out << 'x';
            else
                out << "x^" << i;
        }
        --i;
    }

    return out;
}

} // namespace regina

#endif
