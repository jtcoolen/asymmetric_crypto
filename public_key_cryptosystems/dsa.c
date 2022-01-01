#include <blake2b.h>
#include <dsa.h>
#include <helpers.h>
#include <stdio.h>
#include <string.h>

void DSA_keypair(struct DSA_private_key *private_key,
                 struct DSA_public_key *public_key) {
  mpz_t p, q, mod;
  mpz_inits(p, q, mod, NULL);

  gmp_randstate_t rd_state;
  gmp_randinit_mt(rd_state);

  do {
    mpz_urandomb(q, rd_state, qLEN);
    mpz_nextprime(q, q);
  } while (!(mpz_sizeinbase(q, 2) < qLEN));

  do {
    mpz_urandomb(p, rd_state, pLEN);
    mpz_mod(mod, p, q);
    mpz_sub(p, p, mod);
    mpz_add_ui(p, p, 1);
  } while (!mpz_probab_prime_p(p, 50));

  mpz_t g, g0, pdec, pow;
  mpz_inits(g, g0, pdec, pow, NULL);
  do {
    mpz_urandomb(g0, rd_state,
                 pLEN); // TODO: for now 0 <= g0 <= p - 1, ensure 2 < g0 < p - 2
    mpz_sub_ui(pdec, p, 1);
    mpz_divexact(pow, pdec, q);
    mpz_powm(g, g0, pow, p);
  } while (mpz_cmp_ui(g, 1) == 0);

  mpz_clears(g0, pdec, pow, mod, NULL);

  mpz_t x, y;
  mpz_inits(x, y, NULL);
  mpz_urandomb(x, rd_state, qLEN);
  mpz_powm(y, g, x, p);

  gmp_randclear(rd_state);

  mpz_inits(public_key->p, public_key->q,
            public_key->generator_cyclic_subgroup_order_q,
            public_key->generator_pow, private_key->x, NULL);

  mpz_set(public_key->p, p);
  mpz_set(public_key->q, q);
  mpz_set(public_key->generator_cyclic_subgroup_order_q, g);
  mpz_set(public_key->generator_pow, y);

  mpz_set(private_key->x, x);
  mpz_clears(p, q, g, y, x, NULL);
}

void DSA_sign(const void *plaintext, size_t plaintext_length,
              const struct DSA_public_key *public_key,
              const struct DSA_private_key *private_key,
              struct DSA_signature *signature) {
  mpz_t k, r, s, mpz_hash, xr, kinv;
  mpz_inits(k, r, s, mpz_hash, xr, kinv, NULL);

  gmp_randstate_t rd_state;
  gmp_randinit_mt(rd_state);

  do {
    mpz_urandomb(k, rd_state, qLEN);
    mpz_powm(r, public_key->generator_cyclic_subgroup_order_q, k,
             public_key->p);
    mpz_mod(r, r, public_key->q);

    signed long hash;
    char h[sizeof hash];
    blake2b(h, sizeof hash, NULL, 0, // optional secret key
            plaintext, plaintext_length);

    memcpy(&hash, h, sizeof hash);
    mpz_set_si(mpz_hash, hash);

    mpz_invert(kinv, k, public_key->q);
    mpz_mul(xr, private_key->x, r);
    mpz_add(xr, mpz_hash, xr);
    mpz_mul(s, xr, kinv);

    mpz_mod(s, s, public_key->q);
  } while (mpz_cmp_ui(r, 0) == 0 || mpz_cmp_ui(s, 0) == 0);

  mpz_clears(k, mpz_hash, kinv, xr, NULL);
  gmp_randclear(rd_state);

  mpz_set(signature->r, r);
  mpz_set(signature->s, s);
}

int DSA_check_signature(const void *plaintext, size_t plaintext_length,
                        const struct DSA_signature *signature,
                        const struct DSA_public_key *public_key) {

  if (!(mpz_sgn(signature->r) == 1 &&
        mpz_cmp(signature->r, public_key->q) < 0 &&
        mpz_sgn(signature->s) == 1 &&
        mpz_cmp(signature->s, public_key->q) < 0)) {
    return 0;
  }

  mpz_t w, u1, u2, v, g, y, mpz_hash;
  mpz_inits(w, u1, u2, v, g, y, mpz_hash, NULL);

  mpz_invert(w, signature->s, public_key->q);

  signed long hash;
  char h[sizeof hash];
  blake2b(h, sizeof hash, NULL, 0, // optional secret key
          plaintext, plaintext_length);

  memcpy(&hash, h, sizeof hash);
  mpz_set_si(mpz_hash, hash);

  mpz_mul(u1, mpz_hash, w);
  mpz_mod(u1, u1, public_key->q);
  mpz_mul(u2, signature->r, w);
  mpz_mod(u2, u2, public_key->q);

  mpz_powm(g, public_key->generator_cyclic_subgroup_order_q, u1, public_key->p);
  mpz_powm(y, public_key->generator_pow, u2, public_key->p);
  mpz_mul(v, g, y);
  mpz_mod(v, v, public_key->p);
  mpz_mod(v, v, public_key->q);

  return mpz_cmp(signature->r, v) == 0;
}