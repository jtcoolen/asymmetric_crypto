#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <stdint.h>


// On représent un polynome à coefficients dans le corps fini Fp par
// - le tableau de ses coefficients
// - sa capacité (degré maximal autorisé)
// - la caractéristique p
// - son degré
struct polynomial_in_Fp {
  int64_t *coefficients;
  uint64_t capacity;
  uint64_t characteristic;
  int64_t degree;
};

// Un corps fini de degré q=p^d (q éléments) est construit via le quotient de l'anneau principal Fp[x] par l'idéal maximal engendré par un polynome P de degré d irréductible sur Fp
// Fq=Fp[x]/(P(x))
struct finite_field_element {
  struct polynomial_in_Fp *poly; // écriture d'un élément dans la Fp-base (1,a,...,a^(d-1)) de Fp[x]/(P) avec a racine de P
  struct polynomial_in_Fp *mod; // polynome irréductible
  uint64_t multiplicative_group_order; // ordre du groupe multiplicatif (Fq)*
};

struct polynomial_in_Fp *polynomial_in_Fp_new(uint64_t capacity,
                                                    uint64_t characteristic);

void polynomial_in_Fp_free(struct polynomial_in_Fp *pfe);

int64_t polynomial_in_Fp_degree(const struct polynomial_in_Fp *pfe);

struct polynomial_in_Fp *
polynomial_in_Fp_from_array(int64_t *coefficients, uint64_t capacity,
                               uint64_t characteristic);

struct polynomial_in_Fp *
polynomial_in_Fp_copy(const struct polynomial_in_Fp *pfe);

void polynomial_in_Fp_print(const struct polynomial_in_Fp *pfe);

struct polynomial_in_Fp *
polynomial_in_Fp_add(const struct polynomial_in_Fp *pfe1,
                        const struct polynomial_in_Fp *pfe2);

struct polynomial_in_Fp *
polynomial_in_Fp_subtract(const struct polynomial_in_Fp *pfe1,
                             const struct polynomial_in_Fp *pfe2);

struct polynomial_in_Fp *
polynomial_in_Fp_multiplication(const struct polynomial_in_Fp *pfe1,
                                   const struct polynomial_in_Fp *pfe2);

void polynomial_in_Fp_remainder_mut(
    struct polynomial_in_Fp *pfe, const struct polynomial_in_Fp *divisor);

void polynomial_in_Fp_division(const struct polynomial_in_Fp *pfe1,
                                  const struct polynomial_in_Fp *pfe2,
                                  struct polynomial_in_Fp **remainder,
                                  struct polynomial_in_Fp **quotient);

struct polynomial_in_Fp *
polynomial_in_Fp_gcd(const struct polynomial_in_Fp *pfe1,
                        const struct polynomial_in_Fp *pfe2);

struct polynomial_in_Fp *prime_field_element_gcd_extended(
    const struct polynomial_in_Fp *pfe1,
    const struct polynomial_in_Fp *pfe2,
    struct polynomial_in_Fp **bezout_coefficient_pfe1,
    struct polynomial_in_Fp **bezout_coefficient_pfe2);

void finite_field_element_free(struct finite_field_element *ffe);

struct finite_field_element *
finite_field_element_copy(const struct finite_field_element *ffe);

struct finite_field_element *
finite_field_element_multiplication(const struct finite_field_element *ffe1,
                                    const struct finite_field_element *ffe2);

struct polynomial_in_Fp *prime_field_element_scalar_multiplication(
    int64_t scalar, const struct polynomial_in_Fp *pfe);

struct finite_field_element *
finite_field_element_inverse(const struct finite_field_element *ffe);

struct finite_field_element *
finite_field_element_division(const struct finite_field_element *ffe1,
                              const struct finite_field_element *ffe2);

int finite_field_element_is_generator(const struct finite_field_element *ffe);
#endif
