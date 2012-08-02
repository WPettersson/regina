
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include <algorithm>
#include <functional>
#include "maths/numbertheory.h"
#include "utilities/stlutils.h"

namespace regina {

long reducedMod(long k, long modBase) {
    long ans = k % modBase;
    if (ans < 0) {
        if ((ans + modBase) <= (-ans))
            return (ans + modBase);
    } else if (-(ans - modBase) < ans)
        return (ans - modBase);
    return ans;
}

long gcd(long a, long b) {
    long tmp;
    while (a != b && b != 0) {
        tmp = a;
        a = b;
        b = tmp % b;
    }
    return (a >= 0 ? a : -a);
}

namespace {
    long gcdWithCoeffsInternal(long a, long b, long& u, long& v) {
        // This routine assumes a and b to be non-negative.
        long a_orig = a;
        long b_orig = b;
        u = 1;
        v = 0;
        long uu = 0;
        long vv = 1;
        long tmp1, tmp2, q;
        while (a != b && b != 0) {
            // At each stage:
            // u*(a_orig) + v*(b_orig) = a_curr;
            // uu*(a_orig) + vv*(b_orig) = b_curr.
            tmp1 = u; tmp2 = v;
            u = uu; v = vv;
            q = a / b;
            uu = tmp1 - (q * uu); vv = tmp2 - (q * vv);

            tmp1 = a;
            a = b;
            b = tmp1 % b;
        }

        // a is now our gcd.

        // Put u and v in the correct range.
        if (b_orig == 0)
            return a;

        // We are allowed to add (b_orig/d, -a_orig/d) to (u,v).
        a_orig = -(a_orig / a);
        b_orig = b_orig / a;

        // Now we are allowed to add (b_orig, a_orig), where b_orig >= 0.
        // Add enough copies to put u between 1 and b_orig inclusive.
        long k;
        if (u > 0)
            k = -((u-1) / b_orig);
        else
            k = (b_orig-u) / b_orig;
        if (k) {
            u += k * b_orig;
            v += k * a_orig;
        }

        return a;
    }
}

long gcdWithCoeffs(long a, long b, long& u, long& v) {
    long signA = (a > 0 ? 1 : a == 0 ? 0 : -1);
    long signB = (b > 0 ? 1 : b == 0 ? 0 : -1);
    long ans = gcdWithCoeffsInternal(a >= 0 ? a : -a,
        b >= 0 ? b : -b, u, v);
    u *= signA;
    v *= signB;
    return ans;
}

long lcm(long a, long b) {
    if (a == 0 || b == 0)
        return 0;

    long tmp = regina::gcd(a, b);
    tmp = (a / tmp) * b;
    return (tmp >= 0 ? tmp : -tmp);
}

unsigned long modularInverse(unsigned long n, unsigned long k) {
    if (n == 1)
        return 0;

    long u, v;
    gcdWithCoeffs(n, k % n, u, v);
    // GCD should equal 1, so u*n + k*v = 1.
    // Inverse is v; note that -n < v <= 0.
    // Since n >= 2 now and (n,k) = 1, we know v != 0.
    return v + n;
}

namespace {
    /**
     * Finds the smallest prime factor of the given odd integer.
     * You may specify a known lower bound for this smallest prime factor.
     * If the given integer is prime, 0 is returned.
     *
     * \pre \a n is odd.
     * \pre The smallest prime factor of \a n is known to
     * be at least as large as (and possibly equal to) \a lowerBound.
     * \pre \a lowerBound is odd.
     * 
     * @param n the integer whose smallest prime factor we wish to find.
     * @param lowerBound a known lower bound for this smallest prime factor.
     * @return the smallest prime factor of \a n, or 0 if \a n is prime.
     */
    unsigned long smallestPrimeFactor(unsigned long n,
            unsigned long lowerBound = 1) {
        while (lowerBound * lowerBound <= n) {
            if (n % lowerBound == 0)
                return lowerBound;
            lowerBound += 2;
        }
        return 0;
    }
}

void factorise(unsigned long n, std::list<unsigned long>& factors) {
    // n > 0 is a precondition, but the effects are too unpleasant if
    // it's broken (infinite memory consumption).
    if (n == 0)
        return;

    // First take out all factors of 2.
    while (n % 2 == 0) {
        n = n / 2;
        factors.push_back(2);
    }

    // Run through finding smallest factors.
    unsigned long factor = 3;
    while ((factor = smallestPrimeFactor(n, factor))) {
        factors.push_back(factor);
        n = n / factor;
    }

    // Anything left is prime.
    if (n > 1)
        factors.push_back(n);
}

void primesUpTo(const NLargeInteger& roof, std::list<NLargeInteger>& primes) {
    // First check 2.
    if (roof < 2)
        return;
    primes.push_back(NLargeInteger(2));

    // Run through the rest.
    NLargeInteger current(3);
    while (current <= roof) {
        // Is current prime?
        if (find_if(primes.begin(), primes.end(), regina::stl::compose1(
                bind2nd(std::equal_to<NLargeInteger>(), NLargeInteger::zero),
                bind1st(std::modulus<NLargeInteger>(), current))) ==
                primes.end())
            primes.push_back(current);
        current += 2;
    }
}

} // namespace regina

