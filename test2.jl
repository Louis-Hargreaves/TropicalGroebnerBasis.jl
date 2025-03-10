R, = residue_ring(ZZ, 23)
S, x = polynomial_ring(R, "x")

f = x^2 + 3x + 1

g = x^3 + 3x + 1


R = factor(f*g)

#S = factor_squarefree(f*g)

#T = factor_distinct_deg((x + 1)*g*(x^5+x^3+x+1))
