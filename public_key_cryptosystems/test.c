#include <dsa.h>
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

  // F_q, q=719^6, of multiplicative order 719^6-1 =
  // 138157142546011680

  int64_t prime = 719; // 719 Sophie Germain prime
  struct polynomial_in_Fp *poly =
      polynomial_in_Fp_from_array((int64_t[]){0, 1, 1}, 3, prime);
  struct polynomial_in_Fp *mod =
      polynomial_in_Fp_from_array((int64_t[]){1, 1, 1, 0, 0, 0, 1}, 7, prime);
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

  printf("El Gamal signature:\n");
  struct ElGamal_signature_public_key eg_sig_pubk;
  struct ElGamal_signature_private_key eg_sig_privk;
  ElGamal_signature_keypair(719, 11, &eg_sig_privk,
                            &eg_sig_pubk); // 11 generator of Z/719Z

  char *msg = "hello world!";
  size_t msg_len = 13;
  struct ElGamal_signature *eg_sig =
      ElGamal_sign(msg, msg_len, &eg_sig_pubk, &eg_sig_privk);

  printf("(r,s)=(%ld, %ld)\n", eg_sig->ipfe1, eg_sig->ipfe2);
  printf("Checking signature : %d",
         ElGamal_check_signature(msg, msg_len, eg_sig, &eg_sig_pubk));

  printf("\nDSA signature:\n");
  struct DSA_public_key dsa_pubk;
  struct DSA_private_key dsa_privk;

  DSA_keypair(&dsa_privk, &dsa_pubk);
  printf("\n");
  mpz_out_str(stdout, 10, *dsa_pubk.p);

  struct DSA_signature sig;
  DSA_sign(msg, msg_len, &dsa_pubk, &dsa_privk, &sig);

  return 0;
}