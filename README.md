# TPs du cours de cryptographie asymétrique

## Utilisation

En premier lieu, il faut cloner le sous-module contenant une implémentation de table de hachage via `git submodule update --recursive --remote --init`.

## Arborescence du dossier

- TP 1: dossier `tp1`
  - Calcul de racine carré entière en O(log^3 n) opérations élémentaires (dans les faits en 0(log^2 n) car l'on opère sur des entiers de 64 bits où les opérations bits à bits s'effectuent en temps constant).
  - Algorithme d'Euclide avec le plus petit reste
  - Algorithme de Karatsuba

- TP 2: dossier `finite_field`
Bibliothèque d'opérations sur des corps premiers et finis.

- TP 3: dossier `quadratic_residues`
Symbole de Jacobi et algorithme de Shanks-Tonelli.

- TP 4: dossier `public_key_cryptosystems`
Cryptosystèmes à clé publique: chiffrements et signatures RSA, ElGamal et DSA (exploite le TP1 et la librairie GMP).

- TP 5: dossier `discrete_log`
Algo Baby Step-Giant Step de Shanks, Rho de Pollard (GMP)

- TP 6: dossier `primality_tests`
Tests de Fermat et de Solovay-Strassen.

- TP 7: dossier `factorization`
Rho-Pollard version factorisation et crible quadratique.

- TP 8: dossier `elliptic_curves`
Cardinalité d'une courbe elliptique (algos naif et utilisant Legendre). Multiplication binaire et 2^r-aire. Codage et décodage entier<->point sur courbe. Cryptosystème ElGamal. Chiffrement d'une chaine de caractères via le codage, déchiffrement d'un ensemble de points et retour d'une chaine de caractères.

- dossier `exam_p_plus_1_williams`: algorithme en O(logR log^2 n) (parcours des bits de R, multiplications,additions et inversions mod n)
Le calcul de la suite de Lucas produit un facteur de N=451889: 139
