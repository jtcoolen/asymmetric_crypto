#include <quadratic_residues.h>

int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

int64_t pow_mod(int64_t base, int64_t power, int64_t mod) {
  int64_t res = 1;
  while (power > 0) {
    if ((power & 1) > 0) {
      res = modulo(res * base, mod);
    }
    power >>= 1;
    base = modulo(base * base, mod);
  }
  return res;
}

int64_t modular_inverse(int64_t a, int64_t b) {
  int64_t t, nt, r, nr, q, tmp;
  if (b < 0) {
    b = -b;
  }
  if (a < 0) {
    a = b - (-a % b);
  }
  t = 0;
  nt = 1;
  r = b;
  nr = a % b;

  while (nr != 0) {
    q = r / nr;

    tmp = nt;
    nt = t - q * nt;
    t = tmp;

    tmp = nr;
    nr = r - q * nr;
    r = tmp;
  }
  if (r > 1) {
    return -1;
  }
  if (t < 0) {
    t += b;
  }
  return t;
}

int64_t jacobi(int64_t n, int64_t k) {
  if (k > 0 && k % 2 != 1) {
    return -1;
  }

  int64_t tmp;

  n %= k;
  int64_t t = 1;

  while (n != 0) {
    while (n % 2 == 0) {
      n >>= 1;
      int64_t r = k % 8;
      if (r == 3 || r == 5) {
        t = -t;
      }
    }
    tmp = n;
    n = k;
    k = tmp;
    if (n % 4 == 3 && k % 4 == 3) {
      t = -t;
    }
    n %= k;
  }

  if (k == 1) {
    return t;
  } else {
    return 0;
  }
}

int64_t shanks_tonelli(int64_t a, int64_t p) {
  if (a == 0) {
    return 0;
  }

  if (jacobi(a, p) != 1) {
    return -1;
  }

  int64_t n = 2;
  for (; n < p; n++) {
    if (jacobi(n, p) == -1) {
      break;
    }
  }

  int64_t q = p - 1;
  int64_t s = 0;
  while (q % 2 == 0) {
    q >>= 1;
    s++;
  }

  int64_t r = pow_mod(a, (q + 1) >> 1, p);
  int64_t y = (((r * r) % p) * (modular_inverse(a, p))) % p;
  int64_t b = pow_mod(n, q, p);
  int64_t j = 0;
  int64_t b_pow = 0;

  for (int64_t k = 0; k < s; k++) {
    b_pow = pow_mod(((pow_mod(b, j << 1, p) * y) % p), 1 << (s - 2 - k), p);
    if (b_pow != 1) {
      j ^= (1 << k);
    }
  }

  return ((pow_mod(b, j, p) * r) % p);
}
