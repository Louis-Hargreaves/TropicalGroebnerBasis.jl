using Oscar

#Write something which takes a list of polynomials and a performs Oscar.tighten_simulation on them

#Define the prime p for the p-adic valuation
p = 2
nu_p = tropical_semiring_map(QQ, p)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3], internal_ordering=:deglex)

#Define the polynomials
g1 = x1^2 + t*x2^2 - t^2*x3^2
g2 = t*x1^2 + x2^2 + (t + t^2)*x3^2
g3 = t^4*x1^2 + (t^4 + t^5)*x2^2 + t^3*x3^2 
G = [g1, g2, g3]

#Write out the function which performs Oscar.tighten_simulation

G_red = tighten_simulation.(G, Ref(nu_p))

#Define the ordering on homogeneous functions
M =[-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
o1 = matrix_ordering(Rtx, M)
leading_monomial.(G_red, ordering=o1)

#Sort by leading monomial
G_red = sort(G_red, by = x -> leading_monomial(x, ordering=o1), rev=true)
leading_monomial(G_red[1], ordering=o1)