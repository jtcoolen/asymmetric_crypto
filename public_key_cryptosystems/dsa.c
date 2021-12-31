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
  mpz_urandomb(q, rd_state, 160);

  mpz_urandomb(p, rd_state, 1024);
  mpz_mod(mod, p, q);
  mpz_sub(p, p, mod);
  mpz_add_ui(p, p, 1);

  mpz_t g, g0, pdec, pow;
  mpz_inits(g, g0, pdec, pow, NULL);
  do {
    mpz_urandomb(g0, rd_state,
                 1024); // TODO: for now 0 <= g0 <= p - 1, ensure 2 < g0 < p - 2
    mpz_sub_ui(pdec, p, 1);
    mpz_divexact(pow, p, q);
    mpz_powm(g, g0, pow, p);
  } while (mpz_cmp_ui(g, 1) == 0);

  mpz_clears(g0, pdec, pow, mod, NULL);

  mpz_t x, y;
  mpz_inits(x, y, NULL);
  mpz_urandomb(x, rd_state, 160);
  mpz_powm(y, g, x, p);
  mpz_out_str(stdout, 10, *&p);

  gmp_randclear(rd_state);

  public_key->p = &p;
  public_key->q = &q;
  public_key->generator_cyclic_subgroup_order_q = &g;
  public_key->generator_pow = &y;

  private_key->x = &x;
}

void DSA_sign(const void *plaintext, size_t plaintext_length,
              const struct DSA_public_key *public_key,
              const struct DSA_private_key *private_key,
              struct DSA_signature *signature) {
  mpz_t k, r, s, rtmp, mpz_hash, tmp2, tmp3, tmp4;
  mpz_inits(k, r, s, rtmp, mpz_hash, tmp2, tmp3, tmp4, NULL);

  gmp_randstate_t rd_state;
  gmp_randinit_mt(rd_state);

  do {
    mpz_urandomb(k, rd_state, 160);

    mpz_powm(rtmp, k, *public_key->generator_cyclic_subgroup_order_q,
             *public_key->p);
    mpz_mod(r, rtmp, *public_key->p);

    char h[4];
    blake2b(h, 4, NULL, 0, // optional secret key
            plaintext, plaintext_length);
    unsigned long hash;
    memcpy(&hash, h, 4);
    mpz_set_ui(mpz_hash, hash);

    mpz_invert(rtmp, k, *public_key->q);
    mpz_mul(tmp3, *private_key->x, r);
    mpz_add(tmp2, mpz_hash, tmp3);
    mpz_mul(tmp4, tmp2, tmp3);
    mpz_mod(s, tmp4, *public_key->q);
    printf("\n\ns=");
    mpz_out_str(stdout, 10, s);
    printf("\nr=");
    mpz_out_str(stdout, 10, r);
  } while (mpz_cmp_ui(r, 0) == 0 || mpz_cmp_ui(s, 0) == 0);

  signature = malloc(sizeof(struct DSA_signature));
  signature->r = &r;
  signature->s = &s;
}

int DSA_check_signature(const void *plaintext, size_t plaintext_length,
                        const struct DSA_signature *signature,
                        const struct DSA_public_key *public_key) {
  /*if (!((mpz_sgn(&signature->r) == 1 &&
         mpz_cmp(signature->r, public_key->q) < 0) &&
        (0 < signature->s < public_key->q))) {
    return 0;
  }*/

  mpz_t w, u1, u2, v, tmp1, tmp2, tmp3, mpz_hash;

  mpz_invert(w, *signature->s, *public_key->q);

  char h[4];
  blake2b(h, 4, NULL, 0, // optional secret key
          plaintext, plaintext_length);
  uint32_t hash;
  memcpy(&hash, h, 4);
  mpz_set_ui(mpz_hash, hash);

  mpz_mul(tmp1, mpz_hash, w);
  mpz_mod(u1, tmp1, *public_key->q);
  mpz_mul(tmp2, *signature->r, w);
  mpz_mod(u2, tmp2, *public_key->q);

  mpz_powm(tmp1, *public_key->generator_cyclic_subgroup_order_q, u1,
           *public_key->p);
  mpz_powm(tmp2, *public_key->generator_pow, u2, *public_key->p);
  mpz_mul(tmp3, tmp1, tmp2);
  mpz_mod(v, tmp3, *public_key->q);

  return mpz_cmp(*signature->r, v) == 0;
}