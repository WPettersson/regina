
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

#include <cerrno>
#include <cstdlib>
#include "maths/ninteger.h"
#include "utilities/nthread.h"

// We instantiate both variants of the NInteger template at the bottom
// of this file.

namespace regina {

// Initialize const static data:
template <bool supportInfinity>
const NInteger<supportInfinity> NInteger<supportInfinity>::zero;

template <bool supportInfinity>
const NInteger<supportInfinity> NInteger<supportInfinity>::one(1);

template <>
const NInteger<true> NInteger<true>::infinity(false, false);

// The use of errno in this file should be threadsafe, since (as I
// understand it) each thread gets its own errno.  However, there may be
// thread safety issues regarding locales when using strtol(), in
// particular when another thread changes the locale mid-flight.
// I could be wrong about this.

template <bool supportInfinity>
NInteger<supportInfinity>::NInteger(
        const char* value, int base, bool* valid) :
        large_(0) {
    char* endptr;
    errno = 0;
    small_ = strtol(value, &endptr, base);
    if (errno || *endptr) {
        // Something went wrong.  Try again with large integers.
        // Note that in the case of overflow, we may have errno != 0 but
        // *endptr == 0.
        large_ = new mpz_t;
        if (valid)
            *valid = (mpz_init_set_str(large_, value, base) == 0);
        else
            mpz_init_set_str(large_, value, base);
    } else {
        // All good.
        if (valid)
            *valid = true;
    }
}

template <bool supportInfinity>
NInteger<supportInfinity>::NInteger(
        const std::string& value, int base, bool* valid) :
        large_(0) {
    char* endptr;
    errno = 0;
    small_ = strtol(value.c_str(), &endptr, base);
    if (errno || *endptr) {
        // Something went wrong.  Try again with large integers.
        // Note that in the case of overflow, we may have errno != 0 but
        // *endptr == 0.
        large_ = new mpz_t;
        if (valid)
            *valid = (mpz_init_set_str(large_, value.c_str(), base) == 0);
        else
            mpz_init_set_str(large_, value.c_str(), base);
    } else {
        // All good.
        if (valid)
            *valid = true;
    }
}

template <bool supportInfinity>
std::string NInteger<supportInfinity>::stringValue(int base) const {
    if (isInfinite())
        return "inf";
    else if (large_) {
        char* str = mpz_get_str(0, base, large_);
        std::string ans(str);
        free(str);
        return ans;
    } else {
        // Hmm.  std::setbase() only takes 8, 10 or 16 as i understand it.
        // For now, be wasteful and always go through GMP.
        mpz_t tmp;
        mpz_init_set_si(tmp, small_);

        char* str = mpz_get_str(0, base, tmp);
        std::string ans(str);
        free(str);

        mpz_clear(tmp);
        return ans;
    }
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator =(
        const char* value) {
    InfinityPolicy<supportInfinity>::makeFinite();

    char* endptr;
    errno = 0;
    small_ = strtol(value, &endptr, 10 /* base */);
    if (errno || *endptr) {
        // Something went wrong.  Try again with large integers.
        // Note that in the case of overflow, we may have errno != 0 but
        // *endptr == 0.
        bool result;
        if (large_)
            result = (mpz_set_str(large_, value, 10 /* base */) == 0);
        else {
            large_ = new mpz_t;
            result = (mpz_init_set_str(large_, value, 10 /* base */) == 0);
        }
    } else if (large_) {
        // All good, but we must clear out the old large integer.
        clearLarge();
    }
    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator =(
        const std::string& value) {
    InfinityPolicy<supportInfinity>::makeFinite();

    char* endptr;
    errno = 0;
    small_ = strtol(value.c_str(), &endptr, 10 /* base */);
    if (errno || *endptr) {
        // Something went wrong.  Try again with large integers.
        // Note that in the case of overflow, we may have errno != 0 but
        // *endptr == 0.
        bool result;
        if (large_)
            result = (mpz_set_str(large_, value.c_str(), 10 /* base */) == 0);
        else {
            large_ = new mpz_t;
            result = (mpz_init_set_str(large_, value.c_str(),
                10 /* base */) == 0);
        }
    } else if (large_) {
        // All good, but we must clear out the old large integer.
        clearLarge();
    }
    return *this;
}

template <bool supportInfinity>
std::ostream& operator << (std::ostream& out,
        const NInteger<supportInfinity>& i) {
    if (i.isInfinite())
        out << "inf";
    else if (i.large_) {
        char* str = mpz_get_str(0, 10, i.large_);
        out << str;
        free(str);
    } else
        out << i.small_;
    return out;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator +=(long other) {
    if (isInfinite())
        return *this;
    if (! large_) {
        // Use native arithmetic if we can.
        // The following overflow test is taken from Hackers' Delight,
        // sec. 2-21.
        // Here we use (... & LONG_MIN) to extract the sign bit.
        long ans = small_ + other;
        if ((((ans ^ small_) & (ans ^ other))) & LONG_MIN) {
            // Boom.  It's an overflow.
            // Fall back to large integer arithmetic in the next block.
            forceLarge();
        } else {
            // All good: we're done.
            small_ = ans;
            return *this;
        }
    }

    // And now we're down to large integer arithmetic.
    // The following code should work even if other == LONG_MIN (in which case
    // -other == LONG_MIN also), since passing -other to mpz_sub_ui casts it
    // to an unsigned long (and gives it the correct positive value).
    if (other >= 0)
        mpz_add_ui(large_, large_, other);
    else
        mpz_sub_ui(large_, large_, -other);

    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator -=(long other) {
    if (isInfinite())
        return *this;
    if (! large_) {
        // Use native arithmetic if we can.
        // The following overflow test is taken from Hackers' Delight,
        // sec. 2-21.
        // Here we use (... & LONG_MIN) to extract the sign bit.
        long ans = small_ - other;
        if ((((small_ ^ other) & (ans ^ small_))) & LONG_MIN) {
            // Boom.  It's an overflow.
            // Fall back to large integer arithmetic in the next block.
            forceLarge();
        } else {
            // All good: we're done.
            small_ = ans;
            return *this;
        }
    }

    // And now we're down to large integer arithmetic.
    // The following code should work even if other == LONG_MIN (in which case
    // -other == LONG_MIN also), since passing -other to mpz_add_ui casts it
    // to an unsigned long (and gives it the correct positive value).
    if (other >= 0)
        mpz_sub_ui(large_, large_, other);
    else
        mpz_add_ui(large_, large_, -other);

    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator *=(
        const NInteger<supportInfinity>& other) {
    if (isInfinite())
        return *this;
    else if (other.isInfinite()) {
        InfinityPolicy<supportInfinity>::makeInfinite();
        return *this;
    }
    if (large_) {
        if (other.large_)
            mpz_mul(large_, large_, other.large_);
        else
            mpz_mul_si(large_, large_, other.small_);
    } else if (other.large_) {
        large_ = new mpz_t;
        mpz_init(large_);
        mpz_mul_si(large_, other.large_, small_);
    } else {
        long ans = small_ * other.small_;
        // From Hacker's Delight, sec. 2-12:
        // Overflow iff the following conditions hold:
        if ((other.small_ < 0 && small_ == LONG_MIN) ||
                (other.small_ != 0 && ans / other.small_ != small_)) {
            // Overflow.
            large_ = new mpz_t;
            mpz_init_set_si(large_, small_);
            mpz_mul_si(large_, large_, other.small_);
        } else
            small_ = ans;
    }
    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator *=(long other) {
    if (isInfinite())
        return *this;
    if (large_)
        mpz_mul_si(large_, large_, other);
    else {
        long ans = small_ * other;
        // From Hacker's Delight, sec. 2-12:
        // Overflow iff the following conditions hold:
        if ((other < 0 && small_ == LONG_MIN) ||
                (other != 0 && ans / other != small_)) {
            // Overflow.
            large_ = new mpz_t;
            mpz_init_set_si(large_, small_);
            mpz_mul_si(large_, large_, other);
        } else
            small_ = ans;
    }
    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator /=(
        const NInteger<supportInfinity>& other) {
    if (isInfinite())
        return *this;
    if (other.isInfinite())
        return (*this = 0);
    if (other.isZero()) {
        InfinityPolicy<supportInfinity>::makeInfinite();
        return *this;
    }
    if (other.large_) {
        if (large_) {
            mpz_tdiv_q(large_, large_, other.large_);
            return *this;
        }
        // This is a native C/C++ long.
        // One of four things must happen:
        // (i) |other| > |this|, in which case the result = 0;
        // (ii) this = LONG_MIN and OTHER = -1, in which case the result
        // is the large integer -LONG_MIN;
        // (iii) this = LONG_MIN and OTHER is the large integer -LONG_MIN,
        // in which case the result = -1;
        // (iv) other can be converted to a native long, and the result
        // is a native long also.
        //
        // Deal with the problematic LONG_MIN case first.
        if (small_ == LONG_MIN) {
            if (! mpz_cmp_ui(other.large_,
                    LONG_MIN /* casting to unsigned makes this -LONG_MIN */)) {
                small_ = -1;
                return *this;
            }
            if (! mpz_cmp_si(other.large_, -1)) {
                // The result is -LONG_MIN, which requires large integers.
                // Reduce other while we're at it.
                const_cast<NInteger<supportInfinity>&>(other).forceReduce();
                large_ = new mpz_t;
                mpz_init_set_si(large_, LONG_MIN);
                mpz_neg(large_, large_);
                return *this;
            }
            if (mpz_cmp_ui(other.large_,
                    LONG_MIN /* cast to ui makes this -LONG_MIN */) > 0 ||
                    mpz_cmp_si(other.large_, LONG_MIN) < 0) {
                small_ = 0;
                return *this;
            }
            // other is in [ LONG_MIN, -LONG_MIN ) \ {-1}.
            // Reduce it and use native arithmetic.
            const_cast<NInteger<supportInfinity>&>(other).forceReduce();
            small_ /= other.small_;
            return *this;
        }

        // From here we have this in ( LONG_MIN, -LONG_MIN ).
        if (small_ >= 0) {
            if (mpz_cmp_si(other.large_, small_) > 0 ||
                    mpz_cmp_si(other.large_, -small_) < 0) {
                small_ = 0;
                return *this;
            }
        } else {
            // We can negate, since small_ != LONG_MIN.
            if (mpz_cmp_si(other.large_, -small_) > 0 ||
                    mpz_cmp_si(other.large_, small_) < 0) {
                small_ = 0;
                return *this;
            }
        }

        // We can do this all in native longs from here.
        // Opportunistically reduce other, since we know we can.
        const_cast<NInteger<supportInfinity>&>(other).forceReduce();
        small_ /= other.small_;
        return *this;
    } else
        return (*this) /= other.small_;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator /=(long other) {
    if (isInfinite())
        return *this;
    if (other == 0) {
        InfinityPolicy<supportInfinity>::makeInfinite();
        return *this;
    }
    if (large_) {
        if (other >= 0)
            mpz_tdiv_q_ui(large_, large_, other);
        else {
            // The cast to (unsigned long) makes this correct even if
            // other = LONG_MIN.
            mpz_tdiv_q_ui(large_, large_, - other);
            mpz_neg(large_, large_);
        }
    } else if (small_ == LONG_MIN && other == -1) {
        // This is the special case where we must switch from native to
        // large integers.
        large_ = new mpz_t;
        mpz_init_set_si(large_, LONG_MIN);
        mpz_neg(large_, large_);
    } else {
        // We can do this entirely in native arithmetic.
        small_ /= other;
    }
    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::divByExact(
        const NInteger<supportInfinity>& other) {
    if (other.large_) {
        if (large_) {
            mpz_divexact(large_, large_, other.large_);
            return *this;
        }
        // This is a native C/C++ long.
        // Because we are guaranteed other | this, it follows that
        // other must likewise fit within a native long, or else
        // this == LONG_MIN and other == -LONG_MIN.
        // It also follows that the result must fit within a native long,
        // or else this == LONG_MIN and other == -1.
        if (small_ == LONG_MIN) {
            if (! mpz_cmp_ui(other.large_,
                    LONG_MIN /* casting to unsigned makes this -LONG_MIN */)) {
                // The result is -1, since we have LONG_MIN / -LONG_MIN.
                small_ = -1;
                return *this;
            }

            // At this point we know that other fits within a native long.
            // Opportunistically reduce its representation.
            const_cast<NInteger<supportInfinity>&>(other).forceReduce();

            if (other.small_ == -1) {
                // The result is -LONG_MIN, which requires large integers.
                large_ = new mpz_t;
                mpz_init_set_si(large_, LONG_MIN);
                mpz_neg(large_, large_);
            } else {
                // The result will fit within a native long also.
                small_ /= other.small_;
            }
            return *this;
        }

        // Here we know that other always fits within a native long,
        // and so does the result.
        // Opportunisticaly reduce the representation of other, since
        // we know we can.
        const_cast<NInteger<supportInfinity>&>(other).forceReduce();
        small_ /= other.small_;
        return *this;
    } else {
        // other is already a native int.
        // Use the native version of this routine instead.
        return divByExact(other.small_);
    }
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::divByExact(long other) {
    if (large_) {
        if (other >= 0)
            mpz_divexact_ui(large_, large_, other);
        else {
            // The cast to (unsigned long) makes this correct even if
            // other = LONG_MIN.
            mpz_divexact_ui(large_, large_, - other);
            mpz_neg(large_, large_);
        }
    } else if (small_ == LONG_MIN && other == -1) {
        // This is the special case where we must switch from native to
        // large integers.
        large_ = new mpz_t;
        mpz_init_set_si(large_, LONG_MIN);
        mpz_neg(large_, large_);
    } else {
        // We can do this entirely in native arithmetic.
        small_ /= other;
    }
    return *this;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator %=(
        const NInteger<supportInfinity>& other) {
    if (other.large_) {
        if (large_) {
            mpz_tdiv_r(large_, large_, other.large_);
            return *this;
        }

        // We fit into a native long.  Either:
        // (i) |other| > |this|, in which case the result is just this;
        // (ii) |other| == |this|, in which case the result is 0;
        // (iii) |other| < |this|, in which case we can convert
        // everything to native C/C++ integer arithmetic.
        int res = (small_ >= 0 ?
            mpz_cmp_si(other.large_, small_) :
            mpz_cmp_ui(other.large_, - small_) /* ui cast makes this work
                                                 even if small_ = LONG_MIN */);
        if (res > 0)
            return *this;
        if (res == 0) {
            small_ = 0;
            return *this;
        }

        res = (small_ >= 0 ?
            mpz_cmp_si(other.large_, - small_) :
            mpz_cmp_ui(other.large_, small_));

        if (res < 0)
            return *this;
        if (res == 0) {
            small_ = 0;
            return *this;
        }

        // Everything can be made native integer arithmetic.
        // Opportunistically reduce other while we're at it.
        const_cast<NInteger<supportInfinity>&>(other).forceReduce();
        small_ %= other.small_;
        return *this;
    } else
        return (*this) %= other.small_;
}

template <bool supportInfinity>
NInteger<supportInfinity>& NInteger<supportInfinity>::operator %=(long other) {
    // Since |result| < |other|, whatever happens we can fit the result
    // into a native C/C++ long.
    if (large_) {
        // We can safely cast other to an unsigned long, because the rounding
        // rules imply that (this % LONG_MIN) == (this % -LONG_MIN).
        mpz_tdiv_r_ui(large_, large_, other >= 0 ? other : -other);
        forceReduce();
    } else {
        // All native arithmetic from here.
        small_ %= other;
    }
    return *this;
}

template <bool supportInfinity>
void NInteger<supportInfinity>::raiseToPower(unsigned long exp) {
    // TODO: Fix for natives also.
    if (exp == 0)
        (*this) = one;
    else if (! isInfinite())
        mpz_pow_ui(large_, large_, exp);
}

template <bool supportInfinity>
NInteger<supportInfinity> NInteger<supportInfinity>::gcd(
        const NInteger<supportInfinity>& other) const {
    if (large_) {
        if (other.large_) {
            NInteger<supportInfinity> ans;
            ans.large_ = new mpz_t;
            mpz_init(ans.large_);
            mpz_gcd(ans.large_, large_, other.large_);
            mpz_abs(ans.large_, ans.large_);
            return ans;
        } else {
            NInteger<supportInfinity> ans;
            ans.large_ = new mpz_t;
            mpz_init_set_si(ans.large_, other.small_);
            mpz_gcd(ans.large_, ans.large_, large_);
            mpz_abs(ans.large_, ans.large_);
            return ans;
        }
    } else if (other.large_) {
        NInteger<supportInfinity> ans;
        ans.large_ = new mpz_t;
        mpz_init_set_si(ans.large_, small_);
        mpz_gcd(ans.large_, ans.large_, other.large_);
        mpz_abs(ans.large_, ans.large_);
        return ans;
    }

    // Both integers are native.
    long a = small_;
    long b = other.small_;

    if (a == LONG_MIN) {
        if (b == LONG_MIN) {
            // gcd(a,b) = LONG_MIN, which means we can't make it
            // non-negative without switching to large integers.
            NInteger<supportInfinity> ans;
            ans.large_ = new mpz_t;
            mpz_init_set_si(ans.large_, LONG_MIN);
            mpz_neg(ans.large_, ans.large_);
            return ans;
        }
        a >>= 1; // Won't affect the gcd, but allows us to negate.
    } else if (b == LONG_MIN) {
        b >>= 1; // Won't affect the gcd, but allows us to negate.
    }

    if (a < 0) a = -a;
    if (b < 0) b = -b;

    /**
     * Now everything is non-negative.
     * The following code is based on Stein's binary GCD algorithm.
     */
    if (! a)
        return NInteger<supportInfinity>(b);
    if (! b)
        return NInteger<supportInfinity>(a);

    // Compute the largest common power of 2.
    int pow2;
    for (pow2 = 0; ! ((a | b) & 1); ++pow2) {
        a >>= 1;
        b >>= 1;
    }

    // Strip out all remaining powers of 2 from a and b.
    while (! (a & 1))
        a >>= 1;
    while (! (b & 1))
        b >>= 1;

    while (a != b) {
        // INV: a and b are both odd and non-zero.
        if (a < b) {
            b -= a;
            do
                b >>= 1;
            while (! (b & 1));
        } else {
            a -= b;
            do
                a >>= 1;
            while (! (a & 1));
        }
    }
    return NInteger<supportInfinity>(a << pow2);
}

template <bool supportInfinity>
NInteger<supportInfinity> NInteger<supportInfinity>::lcm(
        const NInteger<supportInfinity>& other) const {
    if (isZero() || other.isZero())
        return zero;

    NInteger<supportInfinity> ans(*this);
    ans.divByExact(gcd(other));
    ans *= other;
    return ans;
}

template <bool supportInfinity>
NInteger<supportInfinity> NInteger<supportInfinity>::gcdWithCoeffs(
        const NInteger<supportInfinity>& other,
        NInteger<supportInfinity>& u, NInteger<supportInfinity>& v) const {
    // TODO: Fix for natives.
    NInteger<supportInfinity> ans;

    // Check for zero coefficients.
    if (isZero()) {
        u = 0L;
        if (other.isZero()) {
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
    if (other.isZero()) {
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
    mpz_gcdext(ans.large_, u.large_, v.large_, large_, other.large_);

    // Ensure the gcd is positive.
    if (ans < 0) {
        ans.negate();
        u.negate();
        v.negate();
    }

    // Get u and v in the correct range.
    NInteger<supportInfinity> addToU(other);
    NInteger<supportInfinity> addToV(*this);
    addToU.divByExact(ans);
    addToV.divByExact(ans);
    if (addToV < 0)
        addToV.negate();
    else
        addToU.negate();

    // We can add (addToU, addToV) to u and v.
    // We also know that addToV is positive.

    // Add enough copies to make v*sign(other) just non-positive.
    NInteger<supportInfinity> copies(v);
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

template <bool supportInfinity>
NInteger<supportInfinity> NInteger<supportInfinity>::divisionAlg(
        const NInteger<supportInfinity>& divisor,
        NInteger<supportInfinity>& remainder) const {
    // TODO: Fix for natives.
    if (divisor == zero) {
        remainder = *this;
        return zero;
    }

    // Preconditions state that nothing is infinite, and we've dealt with d=0.
    // Pass it to GMP.
    NInteger<supportInfinity> quotient;
    mpz_fdiv_qr(quotient.large_, remainder.large_, large_, divisor.large_);

    // The remainder can still be negative (though this will only happen
    // if the divisor is also negative).  In this case we still have
    // more to do.
    if (remainder < zero) {
        remainder -= divisor;
        quotient += 1;
    }

    return quotient;
}

template <bool supportInfinity>
int NInteger<supportInfinity>::legendre(
        const NInteger<supportInfinity>& p) const {
    // TODO: Fix this for natives.
    return mpz_legendre(large_, p.large_);
}

// Instantiate the templates!
template class NInteger<true>;
template class NInteger<false>;

template std::ostream& operator << (std::ostream&, const NInteger<true>&);
template std::ostream& operator << (std::ostream&, const NInteger<false>&);

} // namespace regina
