#include <math.h>
#include <quadratic_sieve.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Détermine la taille de la base de nombres premiers
size_t basis_len(unsigned long long P, mpz_t n) {
  size_t i = 0;
  mpz_t p;
  mpz_init(p);
  mpz_nextprime(p, p);
  i++;
  while (mpz_cmp_ui(p, P) <= 0) {
    mpz_nextprime(p, p);
    if (mpz_legendre(n, p) == 1) {
      i++;
    }
  }
  return i - 1;
}

// Nombres de la forme t^2 - n
struct smooth_candidate {
  mpz_t remaining_factor; // vaut 1 si t^2-n est B friable (B-smooth)
  mpz_t t;
  mpz_t number;
  char *basis_factors_parity;
};


// Macro pour accéder au coefficient (i,j) de la matrice M
// (les coefficients de M sont stockés dans un tableau unidimensionnel
#define MatCoeff(M, i, j, ncol) M[i * (ncol) + j]

void swap_rows(char *matrix, size_t nrow, size_t ncol, size_t row1,
               size_t row2) {
  for (size_t i = 0; i < ncol; i++) {
    size_t tmp = MatCoeff(matrix, row1, i, ncol);
    MatCoeff(matrix, row1, i, ncol) = MatCoeff(matrix, row2, i, ncol);
    MatCoeff(matrix, row2, i, ncol) = tmp;
  }
}

// Pivot de Gauss sur une matrice à coefficients dans F_2
void gaussian_elimination(char *matrix, size_t nrow, size_t ncol) {
  for (int k = 0; k < nrow; k++) {
    if (MatCoeff(matrix, k, k, ncol) ==
        0) { // on permute la colonne courante avec une colonne pivot (avec un 1 sur la diagonale)
      for (int l = k; l < nrow; l++) {
        if (MatCoeff(matrix, l, k, ncol) == 1) {
          swap_rows(matrix, nrow, ncol, l, k);
          break;
        }
      }
    }
    // pour les lignes en dessous du pivot, on soustrait chaque ligne par la ligne du pivot
    for (int i = k + 1; i < nrow; i++) {
      if (MatCoeff(matrix, i, k, ncol)) {
        for (int j = 0; j < ncol; j++)
          MatCoeff(matrix, i, j, ncol) ^= MatCoeff(matrix, k, j, ncol);
      }
    }
  }
}

void transpose(const char *matrix, size_t nrow, size_t ncol, char *transpose) {
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      MatCoeff(transpose, j, i, nrow) = MatCoeff(matrix, i, j, ncol);
    }
  }
}

void quadratic_sieve(mpz_t n, unsigned long long P, unsigned long long A,
                     mpz_t fac1, mpz_t fac2) {
  size_t B_len = basis_len(P, n);
  // printf("len=%ld\n", B_len);
  mpz_t *B = malloc(B_len * sizeof(mpz_t));
  if (B == NULL) {
    return;
  }

  // Initialisation de la base de facteurs premiers <= P
  mpz_t p;
  mpz_init(p);
  size_t i = 0;
  mpz_nextprime(p, p);
  mpz_init(B[i]);
  mpz_set(B[i], p);
  i++;
  for (; i < B_len;) {
    mpz_nextprime(p, p);
    if (mpz_legendre(n, p) == 1) {
      mpz_init(B[i]);
      mpz_set(B[i], p);
      i++;
    }
  }

  struct smooth_candidate *S = malloc(A * sizeof(struct smooth_candidate));
  if (S == NULL) {
    return;
  }

  mpz_t sqrt_n;
  mpz_init(sqrt_n);
  mpz_sqrt(sqrt_n, n);
  mpz_t *fptr;


  // Initialisation de l'ensemble S (sieving interval) du crible
  for (size_t i = 0; i < A; i++) {
    fptr = &S[i].remaining_factor;
    S[i].basis_factors_parity = malloc(B_len * sizeof(char));

    if (S[i].basis_factors_parity == NULL) {
      return;
    }
    memset(S[i].basis_factors_parity, 0, B_len * sizeof(char));
    mpz_init(*fptr);
    mpz_init(S[i].t);
    mpz_init(S[i].number);
    mpz_add_ui(*fptr, sqrt_n, i + 1);
    mpz_set(S[i].t, *fptr);
    mpz_mul(*fptr, *fptr, *fptr);
    mpz_sub(*fptr, *fptr, n);
    mpz_set(S[i].number, *fptr);
  }

  // Crible basique
  for (size_t i = 0; i < B_len; i++) {
    for (size_t j = 0; j < A; j++) {
      while (mpz_divisible_p(S[j].remaining_factor, B[i])) {
        mpz_divexact(S[j].remaining_factor, S[j].remaining_factor, B[i]);
        S[j].basis_factors_parity[i] ^= 1;
      }
    }
  }

  char *smooth_num_vectors = malloc((B_len + 1) * B_len * sizeof(char));
  if (smooth_num_vectors == NULL) {
    return;
  }
  int *smooth_num_indices = malloc((B_len + 1) * sizeof(*smooth_num_indices));
  if (smooth_num_indices == NULL) {
    return;
  }

  mpz_t u, v, diff;
  mpz_inits(u, v, diff, NULL);

  int done = 0;

  // On considère |B|+1 vecteurs de parité d'un nombre B friable consécutifs dans la liste
  // à partir du premier vecteur ayant cette propriété, puis du deuxième, etc... tant que cela
  // ne permet pas de trouver un facteur non trivial de N. C'est le rôle de la variable offset.
  size_t offset = 0;
  while (!done) {
    size_t j = 0;
    for (; offset < A && j < B_len + 1; offset++) {
      if (mpz_cmp_ui(S[offset].remaining_factor, 1) == 0) {
        memcpy(&smooth_num_vectors[j * B_len], S[offset].basis_factors_parity,
               B_len * sizeof(char));
        smooth_num_indices[j] = offset;
        j++;
      }
    }

    // Eviter de sortir de l'intervalle du crible
    if (offset == A - 1 - B_len) {
      return;
    }

    char *smooth_num_vectors2 = malloc((B_len + 1) * B_len * sizeof(char));
    if (smooth_num_vectors2 == NULL) {
      return;
    }

    // Pivot de Gauss sur la matrice base de vecteurs de parité
    transpose(smooth_num_vectors, B_len + 1, B_len, smooth_num_vectors2);
    gaussian_elimination(smooth_num_vectors2, B_len, B_len + 1);

    int f;
    for (f = 0; f < B_len; f++) {
      if (MatCoeff(smooth_num_vectors2, f, f, B_len + 1) != 1) {
        break;
      }
    }

    // "Back substitution"
    for (int k = f - 1; k >= 0; k--) {
      for (int i = k - 1; i >= 0; i--) {
        if (MatCoeff(smooth_num_vectors2, i, k, B_len + 1)) {
          for (int j = 0; j < B_len; j++)
            MatCoeff(smooth_num_vectors2, i, j, (B_len + 1)) ^=
                MatCoeff(smooth_num_vectors2, k, j, (B_len + 1));
        }
      }
    }

    // Initialisation du noyau de la matrice base de vecteurs de parité
    char *nullspace = malloc((B_len + 1) * sizeof(*nullspace));
    memset(nullspace, 0, B_len + 1);
    nullspace[f] = 1;

    for (int i = 0; i < f; i++) {
      nullspace[i] = MatCoeff(smooth_num_vectors2, i, f, (B_len + 1));
    }

    // Le noyau calculé donne une relation de dépendance linéaire entre les vecteurs de parité
    // On obtient u^2 = (t0 * t1 * ... * t|B|)^2 = (t0^2-n) * ... * (t|B|^2-n) mod n = v^2 mod n
    // avec u!=v mod n d'où u-v!=0 mod n et le pgcd de u-v et n peut donner un facteur non trivial de n
    mpz_set_ui(u, 1);
    mpz_set_ui(v, 1);
    for (size_t i = 0; i < B_len + 1; i++) {
      if (nullspace[i] == 1) {
        mpz_mul(u, u, S[smooth_num_indices[i]].t);
        mpz_mul(v, v, S[smooth_num_indices[i]].number);
      }
    }
    mpz_sqrt(v, v);

    mpz_mod(v, v, n);
    mpz_mod(u, u, n);

    mpz_sub(diff, v, u);
    mpz_gcd(fac1, diff, n);
    if (mpz_cmp(n, fac1) != 0) {
      mpz_add(diff, v, u);
      mpz_gcd(fac2, diff, n);
      done = 1;
    }
  }
}
