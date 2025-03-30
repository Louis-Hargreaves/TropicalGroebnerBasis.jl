using Oscar
p = 2
nu_p = tropical_semiring_map(QQ, p)
M = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
o1 = matrix_ordering(Rtx, M)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3, x4, x5) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3])

#Define the polynomials
g1 = x1^2 + t*x2^2 - t^2*x3^2
g2 = t*x1^2 + x2^2 + (t + t^2)*x3^2
g3 = t^4*x1^2 + (t^4 + t^5)*x2^2 + t^3*x3^2 

h1 = t*x1^2*x2 + 5*t^2*x2^2*x3 + t^4*x4^2*x1
h2 = x1^3 + t^7*x2^3 + x5^3

#g4 = t^4*x1^2 + (t^4)*x2^2 + (t+t^3)*x3^2 
G = [g1, g3, g2]
#G = [g3, g2, g1, g4]

function p_minus_t_reduce(G, nu_p, o1)
    G_red = tighten_simulation.(G, Ref(nu_p))
    G_sort = sort(G_red; 
            by=(x->leading_monomial(x, ordering=o1)), 
            lt=((m1, m2) -> cmp(m1, m2) < 0),
            rev=true)
    #Perform first set of iterations
    for i in 1:(length(G_sort)-1)
        g_i = G_sort[i]
        LT_g_i = leading_term(g_i, ordering=o1)
        a_i = leading_exponent(LT_g_i)[2:end]
        t_bi = gen(Rtx, 1) ^ leading_exponent(LT_g_i)[1]
        g_i_ai = LT_g_i / prod(gens(Rtx)[2:end].^a_i)
        for j in i+1:length(G_sort)
            g_j = G_sort[j]
            coeffs_j = collect(coefficients(g_j))
            #Search for a term in g_j with the same a_i
            A_j = collect(exponents(g_j))
            trimmed_A_j = [a[2:end] for a in A_j]
            index = findfirst(x -> x == a_i, trimmed_A_j)
            g_j_ai = coeffs_j[index] * (gen(Rtx, 1) ^ A_j[index][1])
            
            G_sort[j] = ((g_i_ai * g_j)/t_bi) - ((g_j_ai * g_i)/t_bi)
            G_sort[j] = tighten_simulation.(G_sort[j], Ref(nu_p))
        end
    end

    #Perform second set of iterations
    for i in 1:(length(G_sort)-1)
        #g_i = G_sort[i]
        for j in i+1:length(G_sort)
            g_i = G_sort[i]
            g_j = G_sort[j]
            LT_g_j = leading_term(g_j, ordering=o1)
            a_j = leading_exponent(LT_g_j)[2:end]
            t_bj = gen(Rtx, 1) ^ leading_exponent(LT_g_j)[1]
            g_j_aj = LT_g_j / prod(gens(Rtx)[2:end].^a_j)

            #Search for a term in g_i with the same a_j
            A_i = collect(exponents(g_i))
            trimmed_A_i = [a[2:end] for a in A_i]
            index = findfirst(x -> x == a_j, trimmed_A_i)
            g_i_aj = collect(coefficients(g_i))[index] * (gen(Rtx, 1) ^ A_i[index][1])

            #Check if t_bj divides g_i_aj

            toReduce, q = divides(g_i_aj, t_bj)
            if toReduce
                G_sort[i] = ((g_j_aj * g_i)/t_bj) - q * g_j
                G_sort[i] = tighten_simulation(G_sort[i], nu_p)
            end
        end
    end
    return G_sort
end

G_sort = p_minus_t_reduce(G, nu_p, o1)