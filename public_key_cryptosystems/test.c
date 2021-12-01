#include <elgamal.h>
#include <finite_field.h>
#include <rsa.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {

  printf("RSA:\n");
  struct RSA_private_key privk;
  struct RSA_public_key pubk;
  uint64_t p = 443;
  uint64_t q = 743;
  RSA_keypair(p, q, &privk, &pubk);
  uint64_t plaintext = 0xabcf;
  printf("plaintext=%ld, %ld, %ld\n", plaintext, privk.p, p * q);
  const uint64_t ciphertext = RSA_encrypt(plaintext, &pubk);
  printf("ciphertext=%ld\n", ciphertext);
  printf("recovered plaintext=%ld\n", RSA_decrypt(ciphertext, &pubk, &privk));

  printf("ElGamal:\n");
  // ElGamal
  // polisirreducible(('x^6+'x^2+'x+2)*Mod(1,719)) = 1

  // GF(719^6), of multiplicative order 719^6-1 =
  // 138157142546011680

  int64_t prime = 719; // 719 Sophie Germain prime
  struct polynomial_in_Fp *poly =
      polynomial_in_Fp_from_array((int64_t[]){0, 1, 1}, 3, prime);
  struct polynomial_in_Fp *mod = polynomial_in_Fp_from_array(
      (int64_t[]){1, 1, 1, 0, 0, 0, 1}, 7, prime);
  struct finite_field_element ffe;
  ffe.poly = polynomial_in_Fp_copy(poly);
  ffe.mod = polynomial_in_Fp_copy(mod);
  ffe.multiplicative_group_order = 138157142546011680;

  struct ElGamal_private_key eg_privk;
  struct ElGamal_public_key eg_pubk;
  ElGamal_keypair(&ffe, &eg_privk, &eg_pubk);
  struct ElGamal_ciphertext *eg_ciphertext = ElGamal_encrypt(&ffe, &eg_pubk);

  printf("ciphertext: \n");
  polynomial_in_Fp_print(eg_ciphertext->ffe2->poly);
  
  struct finite_field_element *eg_decrypt =
      ElGamal_decrypt(eg_ciphertext, &eg_privk);
      
  printf("\nrecovered plaintext: \n");
  polynomial_in_Fp_print(eg_decrypt->poly);
  printf("\n");

  /*printf("El Gamal signature:\n");

  struct ElGamal_signature *eg_sig = ElGamal_sign("hello world!", 13, &eg_pubk, &eg_privk);
  printf("Checking signature : %d", ElGamal_check_signature("hello world!", 13, eg_sig, &eg_pubk));
  */

  return 0;
}