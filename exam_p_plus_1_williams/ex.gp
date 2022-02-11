\\ calcul du terme de rang n de la suite en O(R log^2 n)
U(n,N)= t=Mod(0,N); u0=Mod(0,N); u1=Mod(1,N); for(i=2,n, t=u1; u1=6*u1-u0; u0=t); u1

N = 451889;
R = 2520;
print("facteur = ", gcd(lift(U(R,N)),N));


\\ calcul du terme de rang R de la suite en O(log R log^2 n)
\\ (parcours des bits de R, analogue à la méthode square & multiply)
\\ On part du principe que u_{2^i} est inversible modulo N pour tout 0<i<log R
U_logR(R,N) = {
  u = Mod(0,N); // U_j,  j=0
  v = Mod(2, N); // V_j, j=0
  u2i=Mod(1,N); // U_{2^i}, i=0
  v2i=Mod(6,N); // V_{2^i}, i=0
  // U_R = U_{a_0 + a_1*2 + ... + a_{log(R)}*2^log(R)}
  // connaissant U_j, j=a_0+a_1 * 2 + ... + a_{i-1} * 2^{i-1}, on calcule U_{j+2^i} si a_i=1
  while(R > 0,
    g=gcd(lift(u2i),N); if(g!=1,return(g)); // Si U_{2^i} n'a pas d'inverse mod N, on retourne le pgcd (peut etre un facteur non trivial)

    // Si a_i == 1 alors on applique les formules U_{n+m} pour n=j et m=2^i
    // et V_{m+n} pour m=j et n=2^i (découle de la formule pour 2 Q^n U_{m+n} avec m=m+n et n=n).
    if(R % 2 == 1, t=u; u=(u * v2i + u2i * v)/2; v = (v2i * u - 2 * t)/u2i); // partie "multiply"
    // partie "square":  U_{2^{i+1}}=U_{2 * 2^i}, idem pour v2i
    u2i = u2i*v2i;
    v2i=v2i^2-2;
    // On considère le bit de R suivant: a_(i+1):
    R >>= 1);
  u
}

print("facteur = ", gcd(lift(U_logR(R,N)),N));