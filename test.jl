using Oscar
include("functions.jl")

#Define the prime p for the p-adic valuation
p = 2
nu_p = tropical_semiring_map(QQ, p)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3])

#Define the polynomials
g1 = x1^2 + t*x2^2 - t^2*x3^2
g2 = t*x1^2 + x2^2 + (t + t^2)*x3^2
g3 = t^4*x1^2 + (t^4 + t^5)*x2^2 + t^3*x3^2 
G = [g1, g2, g3]

#Write out the function which performs Oscar.tighten_simulation
G_red = tighten_simulation.(G, Ref(nu_p))

#Define the ordering on homogeneous functions
M1 = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
o1 = matrix_ordering(Rtx, M1)

G_sort = sort(G_red; lt=((m1, m2) -> cmp_monomial(o1, m1, m2) < 0), rev=true)

#Define an ordering to compare monomials
M2 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -1 0 0 0]
o2 = matrix_ordering(Rtx, M2)

for i in 1:(length(G_sort)-1)
    for j in i+1:length(G_sort)
        g_ia = leading_term(G_sort[i], ordering=o2)
        g_ja = leading_term(G_sort[j], ordering=o2)

        #Ensure that g_ia divides g_ja
        hcf = gcd(g_ia, g_ja)
        G_sort[j] = G_sort[j] * (g_ia/hcf)
        g_ja = leading_term(G_sort[j], ordering=o2)

        #Find the multiple difference and subtract
        t_bi = g_ja/g_ia

        G_sort[j] = G_sort[j] - G_sort[i] * t_bi
        G_sort[j] = tighten_simulation.(G_sort[j], Ref(nu_p))
    end
end
G_sort