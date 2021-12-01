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
  int64_t ffe_order =
      fast_pow(ciphertext->ffe1->mod->characteristic,
               polynomial_in_Fp_degree(ciphertext->ffe1->mod));

  struct finite_field_element *pow = finite_field_element_power(
      ciphertext->ffe1, ffe_order - private_key->pow - 1);
  struct finite_field_element *plaintext =
      finite_field_element_multiplication(pow, ciphertext->ffe2);
  finite_field_element_free(pow);
  return plaintext;
}

struct ElGamal_signature *
ElGamal_sign(const void *plaintext, size_t plaintext_length,
             const struct ElGamal_signature_public_key *public_key,
             const struct ElGaml_signature_private_key *private_key) {
  int64_t k;
  do {
    k = randrange(2, public_key->prime_field_characteristic - 2);
  } while (gcd(k, public_key->prime_field_characteristic - 1) != 1);

  struct ElGamal_signature *signature =
      malloc(sizeof(struct ElGamal_signature));
  signature->ipfe1 =
      pow_mod(public_key->generator, k, public_key->prime_field_characteristic);

  // be careful of signed representation(one's complement)
  uint8_t h[4];
  blake2b(h, 4, NULL, 0, // optional secret key
          plaintext, plaintext_length);
  int64_t hash;
  memcpy(&hash, h, 4);

  signature->ipfe2 =
      modulo((hash - private_key->pow * signature->ipfe1) *
                 modular_inverse(k, public_key->prime_field_characteristic - 1),
             public_key->prime_field_characteristic - 1);

  return signature;
}

int ElGamal_check_signature(
    const void *plaintext, size_t plaintext_length,
    const struct ElGamal_signature *signature,
    const struct ElGamal_signature_public_key *public_key) {
  if (signature->ipfe1 <= 0 ||
      signature->ipfe1 >= public_key->prime_field_characteristic) {
    fprintf(stderr, "wrong signature");
    exit(EXIT_FAILURE);
  }
  if (signature->ipfe2 <= 0 ||
      signature->ipfe2 >= public_key->prime_field_characteristic - 1) {
    fprintf(stderr, "wrong signature");
    exit(EXIT_FAILURE);
  }

  uint8_t h[4];
  blake2b(h, 4, NULL, 0, // optional secret key
          plaintext, plaintext_length);
  int64_t hash;
  memcpy(&hash, h, 4);

  return pow_mod(public_key->generator, hash,
                 public_key->prime_field_characteristic) ==
         modulo(pow_mod(public_key->generator_pow, signature->ipfe1,
                        public_key->prime_field_characteristic) *
                    pow_mod(signature->ipfe1, signature->ipfe2,
                            public_key->prime_field_characteristic),
                public_key->prime_field_characteristic);
}