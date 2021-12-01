
#include <helpers.h>
#include <rsa.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void RSA_keypair(int64_t p, int64_t q, struct RSA_private_key *private_key,
                 struct RSA_public_key *public_key) {
  private_key->p = p;
  private_key->q = q;
  public_key->n = p * q;
  int64_t multiplicative_order_mod_n = (p - 1) * (q - 1);
  int64_t e;
  do {
    e = randrange(2, multiplicative_order_mod_n - 1);
  } while (gcd(e, multiplicative_order_mod_n) != 1);
  public_key->e = e;
  private_key->d = modular_inverse(public_key->e, multiplicative_order_mod_n);
}

uint64_t RSA_encrypt(uint64_t plaintext,
                     const struct RSA_public_key *public_key) {
  return pow_mod(plaintext, public_key->e, public_key->n);
}

uint64_t RSA_decrypt(uint64_t ciphertext,
                     const struct RSA_public_key *public_key,
                     const struct RSA_private_key *private_key) {
  return pow_mod(ciphertext, private_key->d, public_key->n);
}