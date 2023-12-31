\\ Version d'ElGamal qui opère sur des courbes elliptiques définies sur des corps premiers

\\ Le point à l'infini est représenté par oo dans nos calculs
\\ On travaille en coordonnées affines
\\ Si l'on emploie des coordonnées projectives, alors un représentant du point à l'infini est (0,1,0)

ellpointneg(E, P) = {
  if(P == oo, return(P));
  [P[1], -E.a1 * P[1] - E.a3 - P[2]]
}

\\ addition en caractéristique différente de 2 et 3
ellpointadd(E, P1, P2) = {
  if(P1 == oo, return(P2));
  if(P2 == oo, return(P1));
  if(E.p == 2 || E.p == 3,
    s = elladd(E, P1, P2);
    if(s == [0], s = oo);
    return(s));
  P1 = Mod(P1,	Es.p);
  P2 = Mod(P2, Es.p);
  if(P1 == ellpointneg(E, P2), return(oo));
  if(P1 == P2,
    r = (3 * P1[1]^2 + E.a4) / (2 * P1[2]);
    x = r^2 - 2 * P1[1];
    return ([x, r * (P1[1] - x) - P1[2]]);
  );
  r = ((P2[2] - P1[2]) / (P2[1] - P1[1]));
  x = r^2 - P1[1] - P2[1];
  return ([x, r * (P1[1] - x) - P1[2]]);
}

ellpointsub(E, P1, P2) = ellpointadd(E, P1, ellpointneg(E, P2))

\\ duplication en caractéristique quelconque
ellpointdup(E, P) = {
  if(P == oo, return(P));
  if(P[2] == 0, return(ellpointadd(E, P, P))); \\ car il y a un risque que le dénominateur du gradient ci-dessous ne soit pas inversible
  \\ Gradient de la tangente au point P:
  tangentgradient = (3 * P[1]^2 + 2 * E.a2 * P[1] - E.a1 * P[2] + E.a4) / (2 * P[2] + E.a1 * P[1] + E.a3);
  x = tangentgradient^2 + tangentgradient * E.a1 - E.a2 - 2 * P[1];
  [x, -E.a1 * x - E.a3 - tangentgradient * x + tangentgradient * P[1] - P[2]]
}

\\ multiplication binaire d'un point
ellpointmulbin(E, P, n) = {
  my(Q = oo);
  my(k = logint(n, 2));
  while(k != -1, \\ pour visiter le cas k = 0
    Q = ellpointdup(E, Q);
    if((n >> k) % 2 == 1, Q = ellpointadd(E, P, Q));
    k--;
  );
  Q
}

decomp(a) = {
  my(s = 0);
  my(h = a);
  while(h % 2 == 0, h = (h >> 1); s++);
  [s, h]
}


ellpointmul2raryPRECALC(E, P, n, r) = {
  my(multiplesP = vector(2^r+1, i, [0, 0]));
  multiplesP[1] = P;
  multiplesP[2] = ellpointdup(E, P);

  \\ précalcul, j commence à 1 pour remplir toutes les cases d'indice impair
  for(j = 1, 2^(r - 1), multiplesP[2 * j + 1] = ellpointadd(E, multiplesP[2 * j - 1], multiplesP[2]));
  multiplesP
}

\\ multiplication 2^r-aire modifiée d'un point
ellpointmul2rary(E, P, n, r, multiplesP) = {
  my(i = logint(n, 2^r));
  my(z, c, s, h);
  my(Q = oo);
  while(i != -1,
  
    c = (n >> (i * r)) % (2^r);
    if(c != 0, [s, h] = decomp(c));
    
    if(c != 0,
      Q = ellpointmulbin(E, Q, 2^(r - s));
      Q = ellpointadd(E, Q, multiplesP[h]));
      
    if(c == 0, s = r);
    
    Q = ellpointmulbin(E, Q, 2^s);
    i--;
  );
  Q
}

ellisoncurve2(E, P) = (P[2]^2 + E[1] * P[1] * P[2] + E[3] * P[2]) == (P[1]^3 + E[2] * P[1]^2 + E[4] * P[1] + E[5]);

ellcard_naive(E) = {
  my(c = 1, i, j); \\ point à l'infini
  for(i = 0, E.p - 1, for(j = 0, E.p - 1, if(ellisoncurve2(E,[i,j]), c++))); c
}

\\ ne fonctionne qu'en caractéristique différente de 2 et 3 !
ellcard_legendre(E) = {
  my(c = 1, i); \\ point à l'infini
  for(i = 0, E.p - 1, c += kronecker(lift(i^3 + E.a4 * i + E.a6), E.p) + 1); c
}

ellcard2(E) = {
  if(E.p == 2 || E.p == 3, return(ellcard_naive(E)));
  ellcard_legendre(E)
}

\\ Q: public key
\\ d: private key
ElGamal_gen_keys(E, P) = {
  my(d = random(ellcard(E) - 2) + 2);
  [ellpointmulbin(E, P, d), d]
}

ElGamal_encrypt(E, Q, P, msg) = {
  my(k = random(ellcard2(E) - 2) + 2);
  [ellpointmulbin(E, P, k), ellpointadd(E, msg, ellpointmulbin(E, Q, k))]
}

ElGamal_decrypt(E, d, ciphertext) = {
  ellpointsub(E, ciphertext[2], ellpointmul2rary(E, ciphertext[1], d, 4))
}

p = randomprime(2^1000);
a = Mod(2, p);
Es = ellinit([a^4, a^6], a);

time(E,P,n,i,multiplesP)={
  my(start=getabstime());
  ellpointmul2rary(E,P,n,i,multiplesP);
  getabstime()-start;
}

\\ temps en nanosecondes
\\ mesure du temps d'exécution de la multiplication 2^r-aire (on omet le précalcul)
P = random(Es);
for(n=199, 199, for(r = 1, 12, multiplesP = ellpointmul2raryPRECALC(Es, P, n, r); print("r = ", r, " ; temps (en ns) = ", time(Es,P,n,r,multiplesP))));
\*
AVEC PRECALCUL:
1
0
1
1
1
2
3
6
11
16
29
48
93
183
363
727
1456

SANS COMPTER LE PRECALCUL:
r = 1 ; temps (en ns) = 0
r = 2 ; temps (en ns) = 1
r = 3 ; temps (en ns) = 0
r = 4 ; temps (en ns) = 0
r = 5 ; temps (en ns) = 0
r = 6 ; temps (en ns) = 0
r = 7 ; temps (en ns) = 0
r = 8 ; temps (en ns) = 0
r = 9 ; temps (en ns) = 0
r = 10 ; temps (en ns) = 0
r = 11 ; temps (en ns) = 0
r = 12 ; temps (en ns) = 0

On en conclut que le précalcul domine très largement le reste de la multiplication 2^r-aire
Donc ce dernier algo est préférable lorsque l'on est amené à multiplier de nombreuses fois une courbe
*/
\\ mesure du temps d'exécution de la multiplication binaire
time2(E,P,n)={
  my(start=getabstime());
  ellpointmulbin(E,P,n);
  getabstime()-start;
}
P=random(Es); for(n=199, 209, print(time2(Es,P,n)));
\*0
0
1
0
0
0
1
0
0
0
0*/

test() = {
  success = 1;
  while(success,
    p = randomprime(2^10);
    a = Mod(2, p);
    Es = ellinit([a^4, a^6], a);
    if(Es != [],
      print("Es = ", Es);
      P = [0]; \\ point à l'infini selon PARI/GP
      while(P == [0], P = random(Es));
      print("P = ", P);
      
      [Q, d] = ElGamal_gen_keys(Es, P);
      msg = [0];
      while(msg == [0], msg = random(Es));
      print("plaintext = ", msg);
      ciphertext = ElGamal_encrypt(Es, Q, P, msg);
      print("ciphertext = ", ciphertext);
      print("decrypted = ", ElGamal_decrypt(Es, d, ciphertext));
      success = (msg == ElGamal_decrypt(Es, d, ciphertext));
      print("Success: ", success)));
}

\\ boucle infinie! Control+C pour la quitter!
test()