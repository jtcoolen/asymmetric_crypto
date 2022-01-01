#include <blake2b.h>
#include <elgamal.h>
#include <helpers.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ElGamal_keypair(struct finite_field_element *ffe,
                     struct ElGamal_private_key *private_key,
                     struct ElGamal_public_key *public_key) {

  public_key->generator = finite_field_element_copy(ffe);

  int64_t ffe_order =
      fast_pow(ffe->mod->characteristic, polynomial_in_Fp_degree(ffe->mod));
  private_key->pow = randrange(2, ffe_order - 1);

  public_key->generator_pow = finite_field_element_power(ffe, private_key->pow);
}

struct ElGamal_ciphertext *
ElGamal_encrypt(const struct finite_field_element *plaintext,
                const struct ElGamal_public_key *public_key) {

  struct ElGamal_ciphertext *ciphertext =
      malloc(sizeof(struct ElGamal_ciphertext));

  int64_t k = randrange(1, plaintext->multiplicative_group_order);

  ciphertext->ffe1 = finite_field_element_power(public_key->generator, k);

  struct finite_field_element *pow =
      finite_field_element_power(public_key->generator_pow, k);
  ciphertext->ffe2 = finite_field_element_multiplication(plaintext, pow);
  finite_field_element_free(pow);

  return ciphertext;
}

struct finite_field_element *
ElGamal_decrypt(const struct ElGamal_ciphertext *ciphertext,
                const struct ElGamal_private_key *private_key) {
  int64_t ffe_order = fast_pow(ciphertext->ffe1->mod->characteristic,
                               polynomial_in_Fp_degree(ciphertext->ffe1->mod));

  struct finite_field_element *pow = finite_field_element_power(
      ciphertext->ffe1, ffe_order - private_key->pow - 1);
  struct finite_field_element *plaintext =
      finite_field_element_multiplication(pow, ciphertext->ffe2);
  finite_field_element_free(pow);
  return plaintext;
}

void ElGamal_signature_keypair(
    int64_t prime_field_characteristic, int64_t Zp_generator,
    struct ElGamal_signature_private_key *private_key,
    struct ElGamal_signature_public_key *public_key) {

  public_key->prime_field_characteristic = prime_field_characteristic;
  public_key->generator = Zp_generator;
  private_key->pow = randrange(2, prime_field_characteristic - 1);

  public_key->generator_pow =
      pow_mod(Zp_generator, private_key->pow, prime_field_characteristic);
}

struct ElGamal_signature *
ElGamal_sign(const void *plaintext, size_t plaintext_length,
             const struct ElGamal_signature_public_key *public_key,
             const struct ElGamal_signature_private_key *private_key) {
  int64_t k;
  struct ElGamal_signature *signature =
      malloc(sizeof(struct ElGamal_signature));

  do {
    do {
      k = randrange(2, public_key->prime_field_characteristic - 2);
    } while (gcd(k, public_key->prime_field_characteristic - 1) != 1);

    signature->ipfe1 = pow_mod(public_key->generator, k,
                               public_key->prime_field_characteristic);

    uint32_t hash;
    char h[sizeof hash];
    blake2b(h, sizeof hash, NULL, 0, // optional secret key
            plaintext, plaintext_length);

    memcpy(&hash, h, sizeof hash);

    signature->ipfe2 = modulo(
        (hash - private_key->pow * signature->ipfe1) *
            modular_inverse(k, public_key->prime_field_characteristic - 1),
        public_key->prime_field_characteristic - 1);
  } while (signature->ipfe2 == 0);

  return signature;
}

int ElGamal_check_signature(
    const void *plaintext, size_t plaintext_length,
    const struct ElGamal_signature *signature,
    const struct ElGamal_signature_public_key *public_key) {
  if (signature->ipfe1 <= 0 ||
      signature->ipfe1 >= public_key->prime_field_characteristic) {
    return 0;
  }
  if (signature->ipfe2 <= 0 ||
      signature->ipfe2 >= public_key->prime_field_characteristic - 1) {
    return 0;
  }

  uint32_t hash;
  char h[sizeof hash];
  blake2b(h, sizeof hash, NULL, 0, // optional secret key
          plaintext, plaintext_length);

  memcpy(&hash, h, sizeof hash);

  memcpy(&hash, h, sizeof hash);

  return pow_mod(public_key->generator, hash,
                 public_key->prime_field_characteristic) ==
         modulo(pow_mod(public_key->generator_pow, signature->ipfe1,
                        public_key->prime_field_characteristic) *
                    pow_mod(signature->ipfe1, signature->ipfe2,
                            public_key->prime_field_characteristic),
                public_key->prime_field_characteristic);
}