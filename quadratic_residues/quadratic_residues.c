#include <quadratic_residues.h>
#include <assert.h>

// Modulo sur entiers signés
int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

// exponentiation square & multiply
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


// Euclide étendu
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

// Symbole de Legendre
short Legendre(uint64_t a, uint64_t p) {
  assert(a > 0);
  if (a == 1) {
    return 1;
  }
  // ici on a défini a non signé donc positif
  //if (a == -1) {
  //  return (((p - 1) >> 2) & 1) == 0 ? 1 : -1; // Legendre(-1,p) = (-1)^((p-1)/2) donc on se base sur la parité de l'exposant
  //}
  if (a % 2 == 0) {
    short pow = ((p * p - 1) / 8) % 2 == 0 ? 1 : -1;
    return Legendre(a >> 1, p) * pow;
  }
  // Loi de réciprocité quadratique
  short pow = ((a - 1) * (p - 1) / 4) % 2 == 0 ? 1 : -1;
  return Legendre(p % a, a) * pow;
}

// p premier
int64_t shanks_tonelli(int64_t a, int64_t p) {
  if (a == 0) {
    return 0;
  }

  // véification des conditions de l'algorithme
  // On attend que a soit un résidu quadratique modulo p, sinon échec
  if (Legendre(a, p) != 1) {
    return -1;
  }

  int64_t n = 2; // calcul de n non résidu quadratique mod p
  for (; n < p; n++) {
    if (Legendre(n, p) == -1) {
      break;
    }
  }

  // p - 1 = 2^s * q
  int64_t q = p - 1;
  int64_t s = 0;
  while (q % 2 == 0) {
    q >>= 1;
    s++;
  }

  int64_t r = pow_mod(a, (q + 1) >> 1, p); // O(log^3 p)
  int64_t y = (((r * r) % p) * (modular_inverse(a, p))) % p; // O(log^2 p)
  int64_t b = pow_mod(n, q, p); // b^q mod p racine 2^s-ième de l'unité, O(log^3 p)
  int64_t j = 0;
  int64_t b_pow = 0;

  // Calcul de la puissance j de b telle que b^(2*j)*r^2/a = 1 mod p.
  // Pour k>1, on suppose connu j=j_0+j_1*2,...,j_(k-1)*2^(k-1) tel que
  // ((b^j)^2 * r^2 / a)^(2^(s-2-(k-1))) = 1 mod p.
  // D'où j_k = 0 si ((b^j)^2 * r^2 / a)^(2^(s-2-k)) = 1 mod p, j_k=1 sinon.
  // Le cas k=0 correspond au calcul de (r^2/a)^(2^(s-2)), et s'il n'est pas congru à 1 mod p alors j=1.
  // Complexité : O(log p log^3 p) soit O(log^4 p)
  for (int64_t k = 0; k < s; k++) {
    b_pow = pow_mod(((pow_mod(b, j << 1, p) * y) % p), 1 << (s - 2 - k), p); // O(log^3 p)
    if (b_pow != 1) {
      j ^= (1 << k);
    }
  }
  // b^(2*j)*r^2/a = 1 mod p => (b^j*r)^2 = a mod p.
  return ((pow_mod(b, j, p) * r) % p); // O(log p log^2 p)
  // Donc la complexité totale est O(log^4 p)
}
