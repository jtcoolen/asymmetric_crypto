#ifndef RSA_H
#define RSA_H

#include <stdint.h>

struct RSA_private_key {
  int64_t p;
  int64_t q;
  int64_t d;
};

struct RSA_public_key {
  int64_t n; // p*q
  int64_t e;
};

void RSA_keypair(int64_t p, int64_t q, struct RSA_private_key *private_key,
                 struct RSA_public_key *public_key);

uint64_t RSA_encrypt(uint64_t plaintext,
                     const struct RSA_public_key *public_key);

uint64_t RSA_decrypt(uint64_t ciphertext,
                     const struct RSA_public_key *public_key,
                     const struct RSA_private_key *private_key);

#endif