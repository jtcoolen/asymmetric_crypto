#ifndef ELGAMAL_H
#define ELGAMAL_H

#include <stddef.h>
#include <stdint.h>

struct ElGamal_private_key {
  uint64_t pow; // a
};

struct ElGamal_public_key {
  struct finite_field_element *generator;     // g
  struct finite_field_element *generator_pow; // g^a
};

struct ElGamal_ciphertext {
  struct finite_field_element *ffe1; // g^k
  struct finite_field_element *ffe2; // m(g^a)^k
};

struct ElGamal_signature_public_key {
  int64_t prime_field_characteristic;
  int64_t generator;     // g
  int64_t generator_pow; // g^a
};

struct ElGamal_signature_private_key {
  int64_t pow; // a
};

struct ElGamal_signature {
  int64_t ipfe1; // g^k
  int64_t ipfe2; // (H(m)-ar)k^{-1}
};

void ElGamal_keypair(struct finite_field_element *ffe,
                     struct ElGamal_private_key *private_key,
                     struct ElGamal_public_key *public_key);

struct ElGamal_ciphertext *
ElGamal_encrypt(const struct finite_field_element *plaintext,
                const struct ElGamal_public_key *public_key);

struct finite_field_element *
ElGamal_decrypt(const struct ElGamal_ciphertext *ciphertext,
                const struct ElGamal_private_key *private_key);

void ElGamal_signature_keypair(
    int64_t prime_field_characteristic, int64_t Zp_generator,
    struct ElGamal_signature_private_key *private_key,
    struct ElGamal_signature_public_key *public_key);

struct ElGamal_signature *
ElGamal_sign(const void *plaintext, size_t plaintext_length,
             const struct ElGamal_signature_public_key *public_key,
             const struct ElGamal_signature_private_key *private_key);

int ElGamal_check_signature(
    const void *plaintext, size_t plaintext_length,
    const struct ElGamal_signature *signature,
    const struct ElGamal_signature_public_key *public_key);

#endif