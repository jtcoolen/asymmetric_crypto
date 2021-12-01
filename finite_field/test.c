#include <finite_field.h>
#include <stdlib.h>
#include <stdio.h>

int main(void) {
  int64_t b[4] = {1,2,3,4};
  int64_t b1[3] = {1, 2};
  struct prime_field_element *p = prime_field_element_from_array(b, 4, 3);
  struct prime_field_element *q = prime_field_element_from_array(b1, 2, 3);
  prime_field_element_print(p);
  printf("\n");
  prime_field_element_print(q);
  struct prime_field_element *pq = prime_field_element_multiplication(p, q);
  printf("\n");
  prime_field_element_print(pq);
  printf("\n");
  struct prime_field_element *subpq, *rem, *quot;
  prime_field_element_division(p, q, &rem, &quot);
  prime_field_element_print(rem);
  printf("\n");
  prime_field_element_print(quot);
  printf("\n\np-q=");
  subpq = prime_field_element_subtract(p, q);
  prime_field_element_print(subpq);
  printf("\n\n");
  int64_t b2[14] = {0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  struct prime_field_element *r = prime_field_element_from_array(b2, 14, 3);
  prime_field_element_print(r);
  printf("\n");
  int64_t b3[13] = {0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  struct prime_field_element *s = prime_field_element_from_array(b3, 13, 3);
  prime_field_element_print(s);
  printf("\n");
  struct prime_field_element *gcd_pr = prime_field_element_gcd(p, r);
  prime_field_element_print(gcd_pr);
  printf("\n");
  struct finite_field_element *rp = malloc(sizeof(struct finite_field_element));
  rp->poly = prime_field_element_copy(r);
  rp->mod = prime_field_element_copy(p);

  struct finite_field_element *sp = malloc(sizeof(struct finite_field_element));
  sp->poly = prime_field_element_copy(s);
  sp->mod = prime_field_element_copy(p);

  printf("\nrp_inv=");
  struct finite_field_element *rp_inv = finite_field_element_inverse(rp);
  prime_field_element_print(rp_inv->poly);
  printf("\nsp_inv=");
  struct finite_field_element *sp_inv = finite_field_element_inverse(sp);
  prime_field_element_print(sp_inv->poly);
  printf("\n");

  struct finite_field_element *res = finite_field_element_multiplication(rp, sp);
  prime_field_element_print(res->poly);
  printf("\n\ngcd(p,r)");
  struct prime_field_element *u, *v;
  struct prime_field_element *gcd = prime_field_element_gcd_extended(p, r, &u, &v);
  printf("\ngcd=");
  prime_field_element_print(gcd);
  printf("\n");
  prime_field_element_print(u);
  printf("\n");
  prime_field_element_print(v);

  printf("\n");
  struct prime_field_element *prmul = prime_field_element_multiplication(p, r);
  prime_field_element_print(prmul);
  printf("\n");
  struct prime_field_element *gcd_2 = prime_field_element_gcd(prmul, s);
  prime_field_element_print(gcd_2);

  struct prime_field_element *gcd3 = prime_field_element_gcd_extended(prmul, s, &u, &v);
  printf("\ngcd=");
  prime_field_element_print(gcd3);
  printf("\n");
  prime_field_element_print(u);
  printf("\n");
  prime_field_element_print(v);
  printf("\n");
  
  prime_field_element_free(p);
  prime_field_element_free(q);
  prime_field_element_free(pq);
  prime_field_element_free(rem);
  prime_field_element_free(quot);
  prime_field_element_free(r);
  prime_field_element_free(s);
  prime_field_element_free(gcd_pr);
  finite_field_element_free(rp);
  finite_field_element_free(sp);
  finite_field_element_free(res);
  prime_field_element_free(u);
  prime_field_element_free(v);
  prime_field_element_free(gcd);
  finite_field_element_free(rp_inv);
  finite_field_element_free(sp_inv);
  prime_field_element_free(subpq);


  return 0;
}