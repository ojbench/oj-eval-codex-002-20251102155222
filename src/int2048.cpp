// Implementation for sjtu::int2048
#include "include/int2048.h"

namespace sjtu {

// utility
void int2048::normalize() {
  while (!d.empty() && d.back() == 0) d.pop_back();
  if (d.empty()) neg = false;
}

int int2048::cmp_abs(const int2048 &x, const int2048 &y) {
  if (x.d.size() != y.d.size()) return x.d.size() < y.d.size() ? -1 : 1;
  for (size_t i = x.d.size(); i-- > 0;) {
    if (x.d[i] != y.d[i]) return x.d[i] < y.d[i] ? -1 : 1;
  }
  return 0;
}

void int2048::add_abs_to(std::vector<unsigned int> &a, const std::vector<unsigned int> &b) {
  unsigned long long carry = 0;
  size_t n = a.size(), m = b.size();
  size_t L = n > m ? n : m;
  if (a.size() < L) a.resize(L, 0);
  for (size_t i = 0; i < L; ++i) {
    unsigned long long ai = (i < n ? a[i] : 0);
    unsigned long long bi = (i < m ? b[i] : 0);
    unsigned long long sum = ai + bi + carry;
    a[i] = static_cast<unsigned int>(sum % BASE);
    carry = sum / BASE;
  }
  if (carry) a.push_back(static_cast<unsigned int>(carry));
}

void int2048::sub_abs_to(std::vector<unsigned int> &a, const std::vector<unsigned int> &b) {
  // assume |a| >= |b|
  long long carry = 0;
  size_t n = a.size(), m = b.size();
  for (size_t i = 0; i < n; ++i) {
    long long ai = a[i];
    long long bi = (i < m ? b[i] : 0);
    long long diff = ai - bi + carry;
    if (diff < 0) {
      diff += BASE;
      carry = -1;
    } else {
      carry = 0;
    }
    a[i] = static_cast<unsigned int>(diff);
  }
  while (!a.empty() && a.back() == 0) a.pop_back();
}

// FFT multiplication with base 1e3 digits expansion
static void fft(std::vector<std::complex<double>> &a, bool invert) {
  size_t n = a.size();
  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  const double PI = 3.141592653589793238462643383279502884;
  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = 2 * PI / len * (invert ? -1 : 1);
    std::complex<double> wlen(cos(ang), sin(ang));
    for (size_t i = 0; i < n; i += len) {
      std::complex<double> w(1);
      for (size_t j = 0; j < len / 2; ++j) {
        auto u = a[i + j];
        auto v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) {
    for (size_t i = 0; i < n; ++i) a[i] /= static_cast<double>(n);
  }
}

std::vector<unsigned int> int2048::mul_naive(const std::vector<unsigned int> &a, const std::vector<unsigned int> &b) {
  if (a.empty() || b.empty()) return {};
  std::vector<unsigned long long> tmp(a.size() + b.size(), 0);
  for (size_t i = 0; i < a.size(); ++i) {
    unsigned long long carry = 0;
    for (size_t j = 0; j < b.size() || carry; ++j) {
      unsigned __int128 cur = tmp[i + j] + (unsigned __int128)a[i] * (j < b.size() ? b[j] : 0) + carry;
      tmp[i + j] = (unsigned long long)(cur % BASE);
      carry = (unsigned long long)(cur / BASE);
    }
  }
  std::vector<unsigned int> res(tmp.size());
  for (size_t i = 0; i < tmp.size(); ++i) res[i] = (unsigned int)tmp[i];
  while (!res.empty() && res.back() == 0) res.pop_back();
  return res;
}

std::vector<unsigned int> int2048::mul_fft(const std::vector<unsigned int> &a, const std::vector<unsigned int> &b) {
  // expand to base 1e3 digits
  if (a.empty() || b.empty()) return {};
  const int BASE_FFT = 1000;
  std::vector<int> A;
  std::vector<int> B;
  A.reserve(a.size() * 3);
  B.reserve(b.size() * 3);
  auto push_base1000 = [](std::vector<int> &vec, unsigned int x) {
    vec.push_back(x % 1000);
    x /= 1000;
    vec.push_back(x % 1000);
    x /= 1000;
    vec.push_back(x % 1000);
  };
  for (unsigned int v : a) push_base1000(A, v);
  for (unsigned int v : b) push_base1000(B, v);
  size_t n = 1;
  size_t need = A.size() + B.size() - 1;
  while (n < need) n <<= 1;
  std::vector<std::complex<double>> fa(n), fb(n);
  for (size_t i = 0; i < A.size(); ++i) fa[i] = std::complex<double>(A[i], 0);
  for (size_t i = A.size(); i < n; ++i) fa[i] = 0;
  for (size_t i = 0; i < B.size(); ++i) fb[i] = std::complex<double>(B[i], 0);
  for (size_t i = B.size(); i < n; ++i) fb[i] = 0;
  fft(fa, false);
  fft(fb, false);
  for (size_t i = 0; i < n; ++i) fa[i] *= fb[i];
  fft(fa, true);
  std::vector<long long> c(need);
  for (size_t i = 0; i < need; ++i) {
    double vr = fa[i].real();
    long long rounded = (long long)(vr + (vr >= 0 ? 0.5 : -0.5));
    c[i] = rounded;
  }
  // carry in base 1e3
  long long carry2 = 0;
  for (size_t i = 0; i < c.size(); ++i) {
    long long cur = c[i] + carry2;
    if (cur >= 0) {
      carry2 = cur / BASE_FFT;
      c[i] = cur % BASE_FFT;
    } else {
      long long k = (-cur + BASE_FFT - 1) / BASE_FFT;
      cur += k * BASE_FFT;
      carry2 = -k;
      c[i] = cur;
    }
  }
  while (carry2) { c.push_back(carry2 % BASE_FFT); carry2 /= BASE_FFT; }
  while (!c.empty() && c.back() == 0) c.pop_back();
  // regroup to base 1e9 (3 digits of base1e3)
  if (c.empty()) return {};
  if (c.size() % 3 != 0) c.resize((c.size() + 2) / 3 * 3, 0);
  std::vector<unsigned int> res(c.size() / 3);
  for (size_t i = 0; i < res.size(); ++i) {
    unsigned long long x = (unsigned long long)c[3 * i] + (unsigned long long)c[3 * i + 1] * 1000ull + (unsigned long long)c[3 * i + 2] * 1000000ull;
    res[i] = (unsigned int)x; // should be < 1e9
  }
  while (!res.empty() && res.back() == 0) res.pop_back();
  return res;
}

void int2048::divmod_abs(const int2048 &A, const int2048 &B, int2048 &Q, int2048 &R) {
  // assumes B != 0 and both non-negative
  R.d = A.d;
  R.neg = false;
  Q.d.clear();
  Q.neg = false;
  if (cmp_abs(A, B) < 0) {
    R.normalize();
    return;
  }
  size_t n = A.d.size();
  size_t m = B.d.size();
  Q.d.assign(n - m + 1, 0);

  // Normalize divisor and dividend for estimation (Knuth-like)
  unsigned long long f = BASE / ((unsigned long long)(B.d.back()) + 1);
  auto mul_small_vec = [&](std::vector<unsigned int> &x, unsigned long long k) {
    unsigned long long carry = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      unsigned __int128 cur = (unsigned __int128)x[i] * k + carry;
      x[i] = (unsigned int)(cur % BASE);
      carry = (unsigned long long)(cur / BASE);
    }
    if (carry) x.push_back((unsigned int)carry);
  };
  std::vector<unsigned int> u = R.d; // dividend
  std::vector<unsigned int> v = B.d; // divisor
  if (f != 1) {
    mul_small_vec(u, f);
    mul_small_vec(v, f);
  }
  if (u.size() == n) u.push_back(0); // ensure one extra limb for simplicity
  for (size_t k = n - m + 1; k-- > 0;) {
    // estimate qhat using top 2-3 digits
    unsigned long long u2 = (k + m < u.size() ? u[k + m] : 0);
    unsigned long long u1 = u[k + m - 1];
    unsigned long long u0 = (m >= 2 ? u[k + m - 2] : 0);
    unsigned long long v1 = v[m - 1];
    unsigned long long v0 = (m >= 2 ? v[m - 2] : 0);
    unsigned long long dividend = u2 * (unsigned long long)BASE + u1;
    unsigned long long qhat = dividend / v1;
    if (qhat >= BASE) qhat = BASE - 1;
    // refine qhat
    while (true) {
      unsigned long long p1 = v1 * qhat;
      unsigned long long p2 = v0 * qhat;
      unsigned long long left = (dividend - p1) * (unsigned long long)BASE + u0;
      if (p2 <= left) break;
      --qhat;
      if (qhat == 0) break;
    }
    // subtract qhat * v from u at position k
    long long borrow = 0;
    unsigned long long carry3 = 0;
    for (size_t j = 0; j < m; ++j) {
      unsigned __int128 prod = (unsigned __int128)v[j] * qhat + carry3;
      unsigned long long pj = (unsigned long long)(prod % BASE);
      carry3 = (unsigned long long)(prod / BASE);
      long long cur = (long long)u[k + j] - (long long)pj + borrow;
      if (cur < 0) {
        cur += BASE;
        borrow = -1;
      } else {
        borrow = 0;
      }
      u[k + j] = (unsigned int)cur;
    }
    long long cur = (long long)u[k + m] - (long long)carry3 + borrow;
    if (cur < 0) {
      // too much subtracted; add back once
      --qhat;
      unsigned long long c2 = 0;
      for (size_t j = 0; j < m; ++j) {
        unsigned long long sum = (unsigned long long)u[k + j] + v[j] + c2;
        if (sum >= BASE) {
          u[k + j] = (unsigned int)(sum - BASE);
          c2 = 1;
        } else {
          u[k + j] = (unsigned int)sum;
          c2 = 0;
        }
      }
      u[k + m] = (unsigned int)((unsigned long long)u[k + m] + c2);
    } else {
      u[k + m] = (unsigned int)cur;
    }
    Q.d[k] = (unsigned int)qhat;
  }
  // undo normalization for remainder
  R.d.assign(u.begin(), u.begin() + m);
  // divide R by f
  if (f != 1) {
    unsigned long long rem = 0;
    for (size_t i = R.d.size(); i-- > 0;) {
      unsigned __int128 cur = (unsigned __int128)rem * BASE + R.d[i];
      R.d[i] = (unsigned int)(cur / f);
      rem = (unsigned long long)(cur % f);
    }
  }
  Q.normalize();
  R.normalize();
}

// constructors
int2048::int2048() { neg = false; }

int2048::int2048(long long v) {
  if (v < 0) {
    neg = true;
    unsigned long long x = -(unsigned long long)v;
    while (x) {
      d.push_back((unsigned int)(x % BASE));
      x /= BASE;
    }
  } else {
    neg = false;
    unsigned long long x = (unsigned long long)v;
    while (x) {
      d.push_back((unsigned int)(x % BASE));
      x /= BASE;
    }
  }
  normalize();
}

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &o) = default;

// Integer1: read/print
void int2048::read(const std::string &s) {
  d.clear();
  neg = false;
  size_t i = 0;
  while (i < s.size() && isspace((unsigned char)s[i])) ++i;
  if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
    neg = (s[i] == '-');
    ++i;
  }
  while (i < s.size() && s[i] == '0') ++i; // skip leading zeros
  std::vector<int> digits;
  for (; i < s.size(); ++i) {
    if (s[i] >= '0' && s[i] <= '9') digits.push_back(s[i] - '0');
  }
  if (digits.empty()) {
    neg = false;
    return;
  }
  // parse in chunks of BASE_DIG
  for (int pos = (int)digits.size(); pos > 0; pos -= BASE_DIG) {
    int l = std::max(0, pos - BASE_DIG);
    unsigned int val = 0;
    for (int j = l; j < pos; ++j) val = val * 10 + digits[j];
    d.push_back(val);
  }
  normalize();
}

void int2048::print() {
  if (d.empty()) {
    std::cout << 0;
    return;
  }
  if (neg) std::cout << '-';
  std::cout << d.back();
  char buf[16];
  for (size_t i = d.size() - 1; i-- > 0;) {
    std::snprintf(buf, sizeof(buf), "%09u", d[i]);
    std::cout << buf;
  }
}

// Integer1: add/minus
int2048 &int2048::add(const int2048 &b) {
  if (b.d.empty()) return *this;
  if (d.empty()) { *this = b; return *this; }
  if (neg == b.neg) {
    add_abs_to(d, b.d);
  } else {
    int cmp = cmp_abs(*this, b);
    if (cmp == 0) {
      d.clear();
      neg = false;
    } else if (cmp > 0) {
      sub_abs_to(d, b.d);
      // sign remains this->neg
    } else {
      std::vector<unsigned int> tmp = b.d;
      sub_abs_to(tmp, d);
      d.swap(tmp);
      neg = b.neg;
    }
  }
  normalize();
  return *this;
}

int2048 add(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::minus(const int2048 &b) {
  if (b.d.empty()) return *this;
  if (d.empty()) { *this = b; this->neg = !b.neg; return *this; }
  if (neg != b.neg) {
    add_abs_to(d, b.d);
  } else {
    int cmp = cmp_abs(*this, b);
    if (cmp == 0) {
      d.clear();
      neg = false;
    } else if (cmp > 0) {
      sub_abs_to(d, b.d);
      // sign remains this->neg
    } else {
      std::vector<unsigned int> tmp = b.d;
      sub_abs_to(tmp, d);
      d.swap(tmp);
      neg = !b.neg; // opposite sign
    }
  }
  normalize();
  return *this;
}

int2048 minus(int2048 a, const int2048 &b) { return a.minus(b); }

// Integer2 operators
int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const {
  int2048 r(*this);
  if (!r.d.empty()) r.neg = !r.neg;
  return r;
}

int2048 &int2048::operator=(const int2048 &) = default;

int2048 &int2048::operator+=(const int2048 &b) { return add(b); }
int2048 operator+(int2048 a, const int2048 &b) { return a += b; }

int2048 &int2048::operator-=(const int2048 &b) { return minus(b); }
int2048 operator-(int2048 a, const int2048 &b) { return a -= b; }

int2048 &int2048::operator*=(const int2048 &b) {
  if (d.empty() || b.d.empty()) { d.clear(); neg = false; return *this; }
  // choose algorithm
  std::vector<unsigned int> prod;
  size_t n = d.size(), m = b.d.size();
  if (n < 32 || m < 32) prod = mul_naive(d, b.d);
  else prod = mul_fft(d, b.d);
  d.swap(prod);
  neg = (neg != b.neg) && !d.empty();
  normalize();
  return *this;
}

int2048 operator*(int2048 a, const int2048 &b) { return a *= b; }

int2048 &int2048::operator/=(const int2048 &b) {
  // floor division
  if (b.d.empty()) return *this; // undefined, but won't happen per problem
  if (d.empty()) return *this;
  int2048 A = *this; A.neg = false;
  int2048 B = b; B.neg = false;
  int2048 q, r;
  divmod_abs(A, B, q, r);
  bool signs_diff = (neg != b.neg);
  if (signs_diff && !r.d.empty()) {
    // floor adjustment: q = q + (-1)
    // since q is non-negative here
    // q += 1 with sign negative
    // Equivalent to: q = q + 1; then negate after assign sign
    q.add(int2048(1));
  }
  q.neg = signs_diff && !r.d.empty();
  q.normalize();
  *this = q;
  return *this;
}

int2048 operator/(int2048 a, const int2048 &b) { return a /= b; }

  int2048 &int2048::operator%=(const int2048 &b) {
    if (b.d.empty()) return *this; // undefined
    if (d.empty()) return *this;
    int2048 A = *this; A.neg = false;
    int2048 B = b; B.neg = false;
    int2048 q, r;
    divmod_abs(A, B, q, r);
  bool signs_diff = (neg != b.neg);
  if (signs_diff && !r.d.empty()) {
    // rem = r + b, with b positive magnitude B and sign of b.
    // Since r < B and B > 0, r + (-B) = -(B - r)
    int2048 tmpB = B;
    // tmpB >= r always here
    sub_abs_to(tmpB.d, r.d);
    r.d.swap(tmpB.d);
    r.neg = b.neg; // sign of divisor
  } else {
    r.neg = false; // same sign case -> remainder non-negative
  }
    *this = r;
    return *this;
  }

int2048 operator%(int2048 a, const int2048 &b) { return a %= b; }

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s; is >> s; x.read(s); return is;
}
std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.d.empty()) { os << 0; return os; }
  if (x.neg) os << '-';
  os << x.d.back();
  char buf[16];
  for (size_t i = x.d.size() - 1; i-- > 0;) {
    std::snprintf(buf, sizeof(buf), "%09u", x.d[i]);
    os << buf;
  }
  return os;
}

bool operator==(const int2048 &a, const int2048 &b) { return a.neg == b.neg && a.d == b.d; }
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) {
  if (a.neg != b.neg) return a.neg;
  int c = int2048::cmp_abs(a, b);
  return a.neg ? (c > 0) : (c < 0);
}
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
