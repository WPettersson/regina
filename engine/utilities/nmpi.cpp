
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

#include <cstdlib>
#include "utilities/nmpi.h"

namespace regina {

const NLargeInteger NLargeInteger::zero;
const NLargeInteger NLargeInteger::one(1);
const NLargeInteger NLargeInteger::infinity(true, true);

std::string NLargeInteger::stringValue(int base) const {
    if (infinite)
        return "inf";
    else {
        char* str = mpz_get_str(0, base, data);
        std::string ans(str);
        free(str);
        return ans;
    }
}

NLargeInteger NLargeInteger::gcdWithCoeffs(const NLargeInteger& other,
        NLargeInteger& u, NLargeInteger& v) const {
    NLargeInteger ans;

    // Check for zero coefficients.
    if ((*this) == 0) {
        u = 0L;
        if (other == 0) {
            v = 0L;
            // ans is already zero.
            return ans;
        }
        v = 1;
        ans = other;
        if (ans < 0) {
            v.negate();
            ans.negate();
        }
        return ans;
    }
    if (other == 0) {
        v = 0L;
        u = 1;
        ans = *this;
        if (ans < 0) {
            u.negate();
            ans.negate();
        }
        return ans;
    }
    
    // Neither argument is zero.
    // Run the gcd algorithm.
    mpz_gcdext(ans.data, u.data, v.data, data, other.data);

    // Ensure the gcd is positive.
    if (ans < 0) {
        ans.negate();
        u.negate();
        v.negate();
    }

    // Get u and v in the correct range.
    NLargeInteger addToU(other);
    NLargeInteger addToV(*this);
    addToU.divByExact(ans);
    addToV.divByExact(ans);
    if (addToV < 0)
        addToV.negate();
    else
        addToU.negate();

    // We can add (addToU, addToV) to u and v.
    // We also know that addToV is positive.

    // Add enough copies to make v*sign(other) just non-positive.
    NLargeInteger copies(v);
    if (other > 0) {
        // v must be just non-positive.
        if (v > 0) {
            copies -= 1;
            copies /= addToV;
            copies.negate();
            copies -= 1;
        } else {
            copies /= addToV;
            copies.negate();
        }
    }
    else {
        // v must be just non-negative.
        if (v < 0) {
            copies += 1;
            copies /= addToV;
            copies.negate();
            copies += 1;

        } else {
            copies /= addToV;
            copies.negate();
        }
    }
    addToU *= copies;
    addToV *= copies;
    u += addToU;
    v += addToV;
    return ans;
}

std::ostream& operator << (std::ostream& out, const NLargeInteger& large) {
    if (large.infinite)
        out << "inf";
    else {
        char* str = mpz_get_str(0, 10, large.data);
        out << str;
        free(str);
    }
    return out;
}

} // namespace regina
