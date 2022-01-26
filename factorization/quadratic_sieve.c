#include <quadratic_sieve.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


size_t basis_len(unsigned long long P, mpz_t n) {
  size_t i = 0;
  mpz_t p;
  mpz_init(p);
  mpz_nextprime(p, p);
  i++;
  while (mpz_cmp_ui(p, P) <= 0) {
    mpz_nextprime(p, p);
    if (mpz_legendre(n, p) == 1) {
      i++;
    }
  }
  return i - 1;
}

// numbers of the form t^2 - n
struct smooth_candidate {
  mpz_t remaining_factor; // equals to 1 if smooth
  mpz_t t;
  mpz_t number;
  char *basis_factors_parity;
};

#define matrix(i,j) matrix[i * ncol + j]
#define MatCoeff(M, i, j, ncol) M[i * (ncol) + j]

void swap_rows(char *matrix, size_t nrow, size_t ncol, size_t row1, size_t row2) {
  for (size_t i = 0; i < ncol; i++) {
    size_t tmp = matrix(row1, i);
    matrix(row1, i) = matrix(row2, i);
    matrix(row2, i) = tmp;
  }
}

size_t get_pivot(char *matrix, size_t nrow, size_t ncol, size_t l) {
  size_t m = abs(matrix(l, 0));
  if (m == 1) return 0;
  size_t pivot = nrow;
  for (size_t i = 1; i < ncol; i++) {
    if (abs(matrix(l, i)) > m) {
      m = abs(matrix(l, i));
      pivot = i;
      break;
    }
  }
  return pivot;
}

size_t get_pivot2(char *matrix, size_t nrow, size_t ncol, size_t h, size_t k) {
  size_t m = abs(matrix(h, k));
  size_t pivot = 0;
  for (size_t i = h; i < nrow; i++) {
    if (abs(matrix(i, k)) > m) {
      m = abs(matrix(i, k));
      pivot = i;
      break;
    }
  }
  return pivot;
}

// Gaussian elimination over GF(2)
void gaussian_elimination(char *matrix, size_t nrow, size_t ncol) {
  /*size_t h = 1; // pivot ligne
  size_t k = 1; // pivot colonne

  while (h < nrow && k < ncol) {
    size_t pivot = get_pivot2(matrix, nrow, ncol, h,k);
    if (matrix(pivot, k) == 0) {
      k++;
    } else {
      swap_rows(matrix, nrow, ncol, pivot, h);
      for (size_t i = h + 1; i < nrow; i++) {
	matrix(i, k) = 0;
	for (size_t j = k + 1; j < ncol; j++) {
	  matrix(i, j) ^= matrix(h, j);
	}
      }
      h++;
      k++;
    }
    }*/
  for(int k=0; k<nrow; k++) {
    if (matrix(k,k) == 0) { // Swap with pivot if current diagonal is 0
      for(int l = k; l < nrow; l++) {
	if (matrix(l,k) == 1) {
	  swap_rows(matrix, nrow, ncol, l, k);
	  break;
	}
      }
    }
    // For rows below pivot
    for(int i=k+1; i<nrow; i++) {
      // If row can be subtracted, subtract every element (using xor)
      if (matrix(i,k)) {
	for(int j=0; j<ncol; j++)
	  matrix(i,j) ^= matrix(k,j);
	//for(size_t i=0; i<A.size(); i++)
	//   print_range(A[i]);
	//std::cout << '\n';
      }
    }
    }
  
  
}

void transpose(const char *matrix, size_t nrow, size_t ncol, char *transpose) {
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      MatCoeff(transpose, j, i, nrow) = MatCoeff(matrix, i, j, ncol);
    }
  }
}

void right_kernel(char *matrix, size_t n, size_t m) {
  
}

void quadratic_sieve(mpz_t n, unsigned long long P, unsigned long long A) {
  size_t B_len = basis_len(P, n);
  //printf("len=%ld\n", B_len);
  mpz_t *B = malloc(B_len * sizeof(mpz_t));
  if (B == NULL) {
    return;
  }

  mpz_t p;
  mpz_init(p);
  size_t i = 0;
  mpz_nextprime(p, p);
  mpz_init(B[i]);
  mpz_set(B[i], p);
  i++;
  for (; i < B_len;) {
    mpz_nextprime(p, p);
    if (mpz_legendre(n, p) == 1) {
      mpz_init(B[i]);
      mpz_set(B[i], p);
      i++;
    }
  }
  /*for (size_t i = 0; i < B_len; i++) {
    gmp_printf("i=%Zu ", B[i]);
  }
  printf("\n");*/

  struct smooth_candidate *S = malloc(A * sizeof(struct smooth_candidate));
  if (S == NULL) {
    return;
  }

  mpz_t sqrt_n;
  mpz_init(sqrt_n);
  mpz_sqrt(sqrt_n, n);
  mpz_t *fptr;

  for (size_t i = 0; i < A; i++) {
    fptr = &S[i].remaining_factor;
    S[i].basis_factors_parity = malloc(B_len * sizeof(char));

    if (S[i].basis_factors_parity == NULL) {
      return;
    }
    memset(S[i].basis_factors_parity, 0, B_len * sizeof(char));
    mpz_init(*fptr);
    mpz_init(S[i].t);
    mpz_init(S[i].number);
    mpz_add_ui(*fptr, sqrt_n, i + 1);
    mpz_set(S[i].t, *fptr); 
    mpz_mul(*fptr, *fptr, *fptr);
    mpz_sub(*fptr, *fptr, n);
    mpz_set(S[i].number, *fptr);
  }

  // Sieving smooth numbers in the S range
  // TODO: apply the optimisations seen in class
  for (size_t i = 0; i < B_len; i++) {
    for (size_t j = 0; j < A; j++) {
      while (mpz_divisible_p(S[j].remaining_factor, B[i])) {
	mpz_divexact(S[j].remaining_factor, S[j].remaining_factor, B[i]);
	S[j].basis_factors_parity[i] ^= 1; 
      }
    }
  }

  /*for (size_t i = 0; i < A; i++) {
    if (mpz_cmp_ui(S[i].remaining_factor, 1) == 0) {
      gmp_printf("i=%ld, %Zu, t=%Zu \n", i, S[i].number, S[i].t);
    }
    }*/

  char *smooth_num_vectors = malloc((B_len +1) * B_len * sizeof(char));
  if (smooth_num_vectors == NULL) {
    return;
  }
  int *smooth_num_indices = malloc((B_len + 1) * sizeof(*smooth_num_indices));
  if (smooth_num_indices == NULL) {
    return;
  }


  
  mpz_t u, v, fac, diff;
  mpz_inits(u, v, fac, diff, NULL);
  
  int done = 0;

  size_t offset = 0;
  while(!done) {
    size_t j = 0;
    for (; offset < A && j < B_len + 1; offset++) {
      if (mpz_cmp_ui(S[offset].remaining_factor, 1) == 0) {
	//gmp_printf("i=%ld, %Zu, t=%Zu \n", i, S[i].remaining_factor, S[i].t);
	/*for (size_t j = 0; j < B_len; j++) {
	  printf("%d", (int)S[offset].basis_factors_parity[j]); 
	  }*/
	//gmp_printf("   . %Zu\n", S[offset].number);
	memcpy(&smooth_num_vectors[j * B_len], S[offset].basis_factors_parity, B_len * sizeof(char));
	smooth_num_indices[j] = offset;
	j++;
      }
    }

    if (offset == A - 1) {
      return;
    }


  
    /*printf("\n");

    for (size_t i = 0; i < B_len + 1; i++) {
      for (size_t j = 0; j < B_len; j++) {
	printf("%d", smooth_num_vectors[i * B_len + j]); 
      }
      printf("\n");
    }
    printf("\n");*/

    char *smooth_num_vectors2 = malloc((B_len + 1) * B_len * sizeof(char));
    if (smooth_num_vectors2 == NULL) {
      return;
    }

    //printf("\nBlen=%d\n", B_len);

    transpose(smooth_num_vectors, B_len + 1, B_len, smooth_num_vectors2);

    /*printf("\n");
    for (size_t i = 0; i < B_len; i++) {
      for (size_t j = 0; j < B_len+1; j++) {
	printf("%d", smooth_num_vectors2[i *( B_len+1) + j]); 
      }
      printf("\n");
      }*/
  
    /*transpose(smooth_num_vectors, B_len + 1, B_len, smooth_num_vectors2);

      printf("\n");
      for (size_t i = 0; i < B_len ; i++) {
      for (size_t j = 0; j < 2*B_len+1; j++) {
      printf("%d", smooth_num_vectors2[i *( B_len+1) + j]); 
      }
      printf("\n");
      }*/

    gaussian_elimination(smooth_num_vectors2, B_len, B_len+1);

    /*printf("\n");
    for (size_t i = 0; i < B_len; i++) {
      for (size_t j = 0; j < B_len+1; j++) {
	printf("%d", smooth_num_vectors2[i *( B_len+1) + j]); 
      }
      printf("\n");
      }*/
  
    int f;
    for(f=0; f<B_len; f++) {
      if (MatCoeff(smooth_num_vectors2,f,f,B_len+1) != 1) {
	//printf("f=%d\n", f);
	break;
      }
    }

    // Back substitution on upper triangular matrix
    for(int k=f-1; k>=0; k--)
      {
	for(int i=k-1; i>=0; i--)
	  {
	    if (MatCoeff(smooth_num_vectors2, i, k, B_len+1))
	      {
		for(int j=0; j<B_len; j++)
		  MatCoeff(smooth_num_vectors2, i, j, (B_len+1)) ^= MatCoeff(smooth_num_vectors2, k, j, (B_len+1));
	      }
	  }
      }

    char *nullspace = malloc((B_len+1) * sizeof(*nullspace));
    memset(nullspace, 0, B_len+1);
    // Subject to change
    // First free variable is 1, rest are 0
    nullspace[f] = 1;

    for(int i=0; i<f; i++){
      nullspace[i] = MatCoeff(smooth_num_vectors2, i, f, (B_len+1));
    }
    //printf("\nf=%d\n", f);

    /*for(size_t i=0; i<B_len+1; i++) {
      printf(" %d-", nullspace[i]);
    }
    printf("\n");
    */

  
    mpz_set_ui(u, 1);
    mpz_set_ui(v, 1);
    for(size_t i=0; i<B_len+1; i++) {
      if (nullspace[i] == 1)
	{
	  //printf("idx=%d i=%d ", smooth_num_indices[i], i);
	  //gmp_printf("t=%Zu", S[smooth_num_indices[i]].t);
	  mpz_mul(u, u, S[smooth_num_indices[i]].t);
	  mpz_mul(v, v, S[smooth_num_indices[i]].number);
	}
    }
    mpz_sqrt(v, v);

    mpz_mod(v, v, n);
    mpz_mod(u, u, n);

    mpz_sub(diff, v, u);
    mpz_gcd(fac, diff, n);
    if (mpz_cmp(n, fac) != 0) {
      done = 1;
    }
  }

  gmp_printf("u=%Zu, v=%Zu\nGCD=%Zu\n\n", u, v, fac);

}
