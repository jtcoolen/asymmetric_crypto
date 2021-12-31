#ifndef DSA_H
#define DSA_H

#include <gmp.h>
#include <stdlib.h>

struct DSA_public_key {
  mpz_t *p;
  mpz_t *q;
  mpz_t *generator_cyclic_subgroup_order_q;
  mpz_t *generator_pow;
};

struct DSA_private_key {
  mpz_t *x;
};

struct DSA_signature {
  mpz_t *r;
  mpz_t *s;
};

void DSA_keypair(struct DSA_private_key *private_key,
                 struct DSA_public_key *public_key);

void DSA_sign(const void *plaintext, size_t plaintext_length,
              const struct DSA_public_key *public_key,
              const struct DSA_private_key *private_key,
              struct DSA_signature *signature);

int DSA_check_signature(const void *plaintext, size_t plaintext_length,
                        const struct DSA_signature *signature,
                        const struct DSA_public_key *public_key);

#endif