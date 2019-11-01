/**
 * A good portion of this class is taken from the boost library.
 * I wanted to learn some template metaprogramming and also I didn't want
 * to download boost, so I remade it
 */
#pragma once

#include <limits>
#include <type_traits>
#include <iostream>
#include <cassert>

namespace Simplex {
  class division_by_zero{};
  class arithmetic_with_infinity{};

  // loosely inspired by boost::rational
  // some code was copy/pasted and reformatted
  // again I don't really own any of this stuff
  /**
   * Rational supports the usual arithmetic operators and comparisons, as well
   * as a handful of math functions. They're represented internally as two
   * IntTypes, wtih the addition of two special objects positive infinity
   * and negative infinity which compare as you'd expect with any other
   * rational.
   */
  template <typename IntType>
  class rational {
    typedef IntType int_type;
    int_type num, den;
    short inf = 0;
    int_type gcd(int_type, int_type) const;
  public:
    rational(): num(0), den(1) {}

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value, int> = 0
    >
    rational(const NewType &n): num(n), den(1) {}

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value, int> = 0
    >
    rational(const NewType &num, const NewType &den): num(num), den(den)
    { normalize(); }

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value, int> = 0
    >
    rational(const NewType &num, const NewType &den, short inf)
      : num(num), den(den), inf{inf}
    { normalize(); }

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value
                           and std::is_integral<NewType>::value, int> = 0
    >
    rational(const rational<NewType> &r): num(r.num), den(r.den), inf{r.inf} {}

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value, int> = 0
    >
    rational &operator=(const NewType &n)
    { return assign(n, 1); }

    template <typename NewType,
              std::enable_if_t<std::is_integral<NewType>::value, int> = 0
    >
    rational &assign(const NewType &num, const NewType &den, short inf=0)
    { this->num = num; this->den = den; this->inf = inf; return *this; }

    rational &operator+=(const rational &r);
    rational &operator-=(const rational &r);
    rational &operator*=(const rational &r);
    rational &operator/=(const rational &r);

    rational operator+(const rational &r) const;
    rational operator-(const rational &r) const;
    rational operator*(const rational &r) const;
    rational operator/(const rational &r) const;

    rational operator+() const;
    rational operator-() const;

    const rational &operator++();
    const rational &operator--();

    rational operator++(int);
    rational operator--(int);

    constexpr bool operator!() const;

    inline bool operator< (const rational &other) const;
    inline bool operator> (const rational &other) const;
    inline bool operator<=(const rational &other) const;
    inline bool operator>=(const rational &other) const;

    inline bool operator==(const rational &other) const;
    inline bool operator!=(const rational &other) const;

    // normalized form is when gcd(num, den) = 1, den > 0.
    // rationals are always normalized
    void normalize();

    const int_type &numerator() const { return num; }
    const int_type &denominator() const { return den; }

    // warning: attempting to do math with either infinity will instantly
    // throw an exception
    // TODO: can probably cache these values somewhere

    // will compare greater than any other rational
    static rational pos_inf() { return rational(1, 0, 1); }
    // will compare less than any other rational
    static rational neg_inf() { return rational(-1, 0, -1); }
  };

  template <typename I>
  I rational<I>::gcd(I u, I v) const
  {
    if (u < 0) u = -u;
    if (v < 0) v = -v;
    if (v) while ((u %= v) && (v %= u));
    return u + v;
  }

  template <typename I>
  rational<I> &rational<I>::operator+=(const rational<I> &r)
  {
    if (inf) throw arithmetic_with_infinity{};
    if (r.inf) throw arithmetic_with_infinity{};

    // This calculation avoids overflow, and minimises the number of expensive
    // calculations. Thanks to Nickolay Mladenov for this algorithm.
    I r_num = r.num;
    I r_den = r.den;

    I g = gcd(den, r_den);
    den /= g;
    num = num * (r_den / g) + r_num * den; g = gcd(num, g); num /= g;
    den *= r_den/g;

    return *this;
  }

  template <typename I>
  rational<I> &rational<I>::operator-=(const rational<I> &r)
  {
    if (inf) throw arithmetic_with_infinity{};
    if (r.inf) throw arithmetic_with_infinity{};

    I r_num = r.num;
    I r_den = r.den;

    I g = gcd(den, r_den);
    den /= g;
    num = num * (r_den / g) - r_num * den;
    g = gcd(num, g);
    num /= g;
    den *= r_den/g;

    return *this;
  }

  template <typename I>
  rational<I> &rational<I>::operator*=(const rational<I> &r)
  {
    if (inf) throw arithmetic_with_infinity{};
    if (r.inf) throw arithmetic_with_infinity{};

    I r_num = r.num;
    I r_den = r.den;

    // Avoid overflow and preserve normalization
    I gcd1 = gcd(num, r_den);
    I gcd2 = gcd(r_num, den);
    num = (num/gcd1) * (r_num/gcd2);
    den = (den/gcd2) * (r_den/gcd1);
    return *this;
  }

  template <typename I>
  rational<I> &rational<I>::operator/=(const rational<I> &r)
  {
    if (inf) throw arithmetic_with_infinity{};
    if (r.inf) throw arithmetic_with_infinity{};

    // Protect against self-modification
    I r_num = r.num;
    I r_den = r.den;

    // Avoid repeated construction
    I zero(0);

    // Trap division by zero
    if (r_num == zero)
      throw division_by_zero{};
    if (num == zero)
      return *this;

    // Avoid overflow and preserve normalization
    I gcd1 = gcd(num, r_num);
    I gcd2 = gcd(r_den, den);
    num = (num/gcd1) * (r_den/gcd2);
    den = (den/gcd2) * (r_num/gcd1);

    if (den < zero) {
        num = -num;
        den = -den;
    }
    return *this;
  }

  template <typename I>
  rational<I> rational<I>::operator+(const rational<I> &r) const
  { rational<I> s(*this); return s += r; }
  
  template <typename I>
  rational<I> rational<I>::operator-(const rational<I> &r) const
  { rational<I> s(*this); return s -= r; }

  template <typename I>
  rational<I> rational<I>::operator*(const rational<I> &r) const
  { rational<I> s(*this); return s *= r; }

  template <typename I>
  rational<I> rational<I>::operator/(const rational<I> &r) const
  { rational<I> s(*this); return s /= r; }

  template <typename I>
  rational<I> rational<I>::operator+() const { return *this; }

  template <typename I>
  rational<I> rational<I>::operator-() const
  { return rational<I>(-num, den, -inf); }

  template <typename I>
  const rational<I> &rational<I>::operator++()
  { 
    if (inf) throw arithmetic_with_infinity{};
    num += den; return *this;
  }

  template <typename I>
  const rational<I> &rational<I>::operator--()
  {
    if (inf) throw arithmetic_with_infinity{};
    num -= den; return *this;
  }

  template <typename I>
  rational<I> rational<I>::operator++(int)
  { rational<I> t(*this); ++(*this); return t; }

  template <typename I>
  rational<I> rational<I>::operator--(int)
  { rational<I> t(*this); --(*this); return t; }

  template <typename I>
  constexpr bool rational<I>::operator!() const { return !num; }

  template <typename int_type>
  inline bool rational<int_type>::operator<(const rational<int_type> &r) const
  {
    if (inf < r.inf) return true;
    if (inf > r.inf) return false;
    if (inf) return false;

    // Avoid repeated construction
    int_type const zero(0);

    // This should really be a class-wide invariant.  The reason for these
    // checks is that for 2's complement systems, INT_MIN has no corresponding
    // positive, so negating it during normalization keeps it INT_MIN, which
    // is bad for later calculations that assume a positive denominator.
    assert(this->den > zero and "Divide by zero ?? ");
    assert(r.den > zero and "Divide by zero ?? ");

    // Determine relative order by expanding each value to its simple continued
    // fraction representation using the Euclidian GCD algorithm.
    struct { int_type  n, d, q, r; }
     ts = { this->num, this->den, static_cast<int_type>(this->num / this->den),
     static_cast<int_type>(this->num % this->den) },
     rs = { r.num, r.den, static_cast<int_type>(r.num / r.den),
     static_cast<int_type>(r.num % r.den) };
    unsigned reverse = 0u;

    // Normalize negative moduli by repeatedly adding the (positive) denominator
    // and decrementing the quotient.  Later cycles should have all positive
    // values, so this only has to be done for the first cycle.  (The rules of
    // C++ require a nonnegative quotient & remainder for a nonnegative dividend
    // & positive divisor.)
    while ( ts.r < zero )  { ts.r += ts.d; --ts.q; }
    while ( rs.r < zero )  { rs.r += rs.d; --rs.q; }

    // Loop through and compare each variable's continued-fraction components
    for (;;) {
      // The quotients of the current cycle are the continued-fraction
      // components.  Comparing two c.f. is comparing their sequences,
      // stopping at the first difference.
      if (ts.q != rs.q) {
          // Since reciprocation changes the relative order of two variables,
          // and c.f. use reciprocals, the less/greater-than test reverses
          // after each index.  (Start w/ non-reversed @ whole-number place.)
          return reverse ? ts.q > rs.q : ts.q < rs.q;
      }

      // Prepare the next cycle
      reverse ^= 1u;

      if ( (ts.r == zero) || (rs.r == zero) ) {
        // At least one variable's c.f. expansion has ended
        break;
      }

      ts.n = ts.d;         ts.d = ts.r;
      ts.q = ts.n / ts.d;  ts.r = ts.n % ts.d;
      rs.n = rs.d;         rs.d = rs.r;
      rs.q = rs.n / rs.d;  rs.r = rs.n % rs.d;
    }

    // Compare infinity-valued components for otherwise equal sequences
    if (ts.r == rs.r) {
      // Both remainders are zero, so the next (and subsequent) c.f.
      // components for both sequences are infinity.  Therefore, the sequences
      // and their corresponding values are equal.
      return false;
    } else {
      // Exactly one of the remainders is zero, so all following c.f.
      // components of that variable are infinity, while the other variable
      // has a finite next c.f. component.  So that other variable has the
      // lesser value (modulo the reversal flag!).
      return (ts.r != zero) != static_cast<bool>(reverse);
    }
  }

  template <typename I>
  inline bool rational<I>::operator> (const rational<I> &other) const
  { return other < *this; }

  template <typename I>
  inline bool rational<I>::operator<=(const rational<I> &other) const
  { return !(other < *this); }

  template <typename I>
  inline bool rational<I>::operator>=(const rational<I> &other) const
  { return !(*this < other); }

  template <typename I>
  inline bool rational<I>::operator==(const rational<I> &other) const
  {
    return (this->inf and this->inf == other.inf)
       or  (!this->inf and !other.inf and
             this->num == other.num   and this->den == other.den);
  }

  template <typename I>
  inline bool rational<I>::operator!=(const rational<I> &other) const
  { return !(*this == other); }

  // normalized form is when gcd(num, den) = 1, den > 0.
  // rationals are always normalized
  template <typename I>
  void rational<I>::normalize() {
    if (inf) return;

    // Avoid repeated construction
    I zero(0);

    if (den == zero) throw division_by_zero{};

    // Handle the case of zero separately, to avoid division by zero
    if (num == zero) { den = I(1); return; }

    I g = gcd(num, den);

    num /= g;
    den /= g;

    // Ensure that the denominator is positive
    if (den < zero) {
        num = -num;
        den = -den;
    }

    // ...But acknowledge that the previous step doesn't always work.
    // (Nominally, this should be done before the mutating steps, but this
    // member function is only called during the constructor, so we never have
    // to worry about zombie objects.)
    if (den < zero) throw -100; // ???
  }

  template <typename I>
  inline const rational<I> conj(const rational<I> &x) { return x; }

  template <typename I>
  inline const rational<I> real(const rational<I> &x) { return x; }

  template <typename I>
  inline const rational<I> imag(const rational<I> &x) { return 0; }

  template <typename I>
  inline const rational<I> abs(const rational<I> &x) { return x > 0? x : -x; }

  template <typename I>
  inline const rational<I> abs2(const rational<I> &x) { return x * x; }

  template <typename I>
  inline const rational<I> sqrt(const rational<I> &x) { 
    std::cout << x << std::endl;
    I num_root(-1), den_root(-1);
    for (I i(0); i * i <= x.numerator(); ++i) {
      if ((x.numerator() % i == 0) and (x.numerator() / i == i))
        num_root = i;
    }
    if (num_root == -1)
      throw -213498; // tried to take square root of non square
    for (I i(0); i * i <= x.denominator(); ++i) {
      if ((x.denominator() % i == 0) and (x.denominator() / i == i))
        den_root = i;
    }
    if (num_root == -1)
      throw -213498; // tried to take square root of non square
    return rational<I>(num_root, den_root);
  }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator<(const T &a, const rational<I> &b)
  { return b > a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator<=(const T &a, const rational<I> &b)
  { return b >= a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator==(const T &a, const rational<I> &b)
  { return b == a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator>=(const T &a, const rational<I> &b)
  { return b <= a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator>(const T &a, const rational<I> &b)
  { return b < a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline bool operator!=(const T &a, const rational<I> &b)
  { return b != a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline rational<I> operator+(const T &a, const rational<I> &r)
  { rational<I> s(r); return s += a; }
  
  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline rational<I> operator-(const T &a, const rational<I> &r)
  { rational<I> s(r); return s -= a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline rational<I> operator*(const T &a, const rational<I> &r)
  { rational<I> s(r); return s *= a; }

  template <typename T, typename I,
            std::enable_if_t<std::is_integral<T>::value, int> = 0
  >
  inline rational<I> operator/(const T &a, const rational<I> &r)
  { rational<I> s(r); return s /= a; }

  template <typename I>
  std::ostream &operator<<(std::ostream &os, const Simplex::rational<I> &r)
  { os << r.numerator() << " / " << r.denominator(); return os; }

}

#include "../Eigen/Core"

namespace Eigen {
  template <typename IntType>
  struct NumTraits<Simplex::rational<IntType>> 
    : Eigen::GenericNumTraits<Simplex::rational<IntType>>
  {
    typedef Simplex::rational<IntType> Real;
    typedef Simplex::rational<IntType> NonInteger;
    typedef Simplex::rational<IntType> Literal;
    typedef Simplex::rational<IntType> Nested;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }
    static inline Simplex::rational<IntType> highest() {
      return std::numeric_limits<IntType>::max();
    }
    static inline Simplex::rational<IntType> lowest() {
      return std::numeric_limits<IntType>::max();
    }
    static inline int digits10() {
      return std::numeric_limits<IntType>::digits10;
    }

    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 0,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 3,
      MulCost = 3
    };
  };
}


