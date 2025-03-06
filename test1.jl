using Oscar
include("functions.jl")

#Write something which takes a list of polynomials and a performs Oscar.tighten_simulation on them

#Define the prime p for the p-adic valuation
p = 2
nu_p = tropical_semiring_map(QQ, p)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3])

#Define the polynomials
g1 = 4*x1 + 8*x2
tighten_simulation(g1, nu_p)