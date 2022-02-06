# TPs du cours de cryptographie asymétrique

## Utilisation

En premier lieu, il faut cloner le sous-module contenant une implémentation de table de hachage via `git submodule update --recursive --remote --init`.

## Arborescence du dossier

[ x ] TP 1: dossiers <integer_sqrt> et <karatsuba>
- Calcul de racine carré entière en O(log^3 n) opérations élémentaires (dans les faits en 0(log^2 n) car l'on opère sur des entiers de 64 bits où les opérations bits à bits s'effectuent en temps O(1)).
- Algorithme de Karatsuba

[ x ] TP 2: dossier <finite_field>
Bibliothèque d'opérations sur des corps premiers et finis.

[ x ] TP 3: dossier <quadratic residues>
Symbole de Jacobi et algorithme de Shanks-Tonelli.

[ x ] TP 4: dossier <public_key_cryptosystems>
Cryptosystèmes à clé publique: chiffrements et signatures RSA, ElGamal et DSA (exploite le TP1 et la librairie GMP).

[ x ] TP 5: dossier <discrete_log>
Algo Baby Step-Giant Step de Shanks, Rho de Pollard (GMP)

[ x ] TP 6: dossier <primality_tests>
Tests de Fermat et de Solovay-Strassen.

[ x ] TP 7: dossier <factorization>
Rho-Pollard version factorisation et crible quadratique.

[ x ] TP 8: dossier <elliptic_curves>
Cardinalité d'une courbe elliptique (algos naif et utilisant Legendre). Multiplication binaire et 2^r-aire. Codage et décodage entier<->point sur courbe. Cryptosystème ElGamal.