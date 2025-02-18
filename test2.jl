using Oscar

#Define the prime p for the p-adic valuation
p = 2
nu_p = tropical_semiring_map(QQ, p)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3) = Oscar.valued_ring(nu_p)[:t, :x1 , :x2, :x3]

#Define the polynomials
g1 = x1^2 + t*x2^2 - t^2*x3^2
g2 = t*x1^2 + x2^2 + (t + t^2)*x3^2
monomial_ordering(R, :lex)

sort([x1, 2, 3], lt=cmp(lex([t, x1, x2, x3])) == 1)
#sort([1, 2, 3], rev=true)

sort([1,2,3], lt= <)