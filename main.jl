using Oscar

#Start with the polynomial ring
Rtx, (t, x1, x2, x3, x4) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3, :x4])

#Define the tropical semiring map / valuation
p = 3
nu_p = tropical_semiring_map(QQ, p)

#Define the ordering
M = [-1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]
o1 = matrix_ordering(Rtx, M)

#Define the polynomials
g1 = x1^2 + t*x2^2 - t^2*x3^2
g2 = t*x1^2 + x2^2 + (t + t^2)*x3^2
g3 = t^4*x1^2 + (t^4 + t^5)*x2^2 + t^3*x3^2 

h1 = x1*x2*x3 + 5*t^7*x1^3 + 17*x2^2*x3

G = [g1, g2, g3]
H = [h1]

function reduction_same_degree(G, nu_p, o1, sort_req=true)
    if sort_req
        G_red = tighten_simulation.(G, Ref(nu_p))
    end
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
            if !isnothing(index)
                g_j_ai = coeffs_j[index] * (gen(Rtx, 1) ^ A_j[index][1])
                G_sort[j] = ((g_i_ai * g_j)/t_bi) - ((g_j_ai * g_i)/t_bi)
                if !iszero(G_sort[j])
                    #To do: Fix the tighten_simulation function
                    G_sort[j] = tighten_simulation.(G_sort[j], Ref(nu_p))
                end
            end
        end
    end
    #Perform second set of iterations
    for i in 1:(length(G_sort)-1)
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
            if !isnothing(index)
                g_i_aj = collect(coefficients(g_i))[index] * (gen(Rtx, 1) ^ A_i[index][1])
                #Check if t_bj divides g_i_aj
                toReduce, q = divides(g_i_aj, t_bj)
                if toReduce
                    G_sort[i] = ((g_j_aj * g_i)/t_bj) - q * g_j
                    G_sort[i] = tighten_simulation(G_sort[i], nu_p)
                end
            end
        end
    end
    return G_sort
end

#AI tool used for this function
function integer_solutions(n, k)
    if k == 1
        return [[n]]  # Only one variable left, so it must be n
    end

    solutions = []
    for x in 0:n
        for rest in integer_solutions(n - x, k - 1)
            push!(solutions, [x; rest])  # Combine x with the rest of the solution
        end
    end
    return solutions
end

function reduction_at_once(G, H, nu_p, o1)
    E = []
    #Record the highest degree
    leading_exponents = leading_exponent.(G, ordering=o1)
    d = sum(collect(exponents(H[1]))[1][2:end])
    #Loop through all alpha = degree d 
    nvars_in_x = length(leading_exponents[1]) - 1
    #Generates all |α| = d
    alpha = integer_solutions(d, nvars_in_x)
    for a in alpha
        #Pick g
        for (g, expv) in zip(G, leading_exponents)
            #Checking if t^b x^a is in LT_G
            if all(a .>=  expv[2:end])
                tb_xa = gen(Rtx, 1) ^ expv[1] * prod(gens(Rtx)[2:end].^a)
                push!(E, g * (tb_xa / prod(gens(Rtx).^expv)) )
                break
            end
        end
    end
    #Perform the reduction
    H_U_E = vcat(H, E)
    H_U_E_tightened = tighten_simulation.(H_U_E, Ref(nu_p))
    sigma = sortperm(H_U_E_tightened; 
            by=(x->leading_monomial(x, ordering=o1)), 
            lt=((m1, m2) -> cmp(m1, m2) < 0),
            rev=true)
    H_U_E_red = reduction_same_degree(H_U_E_tightened, nu_p, o1)
    return [H_U_E_red[i] for i in sigma if sigma[i] in 1:length(H)]
end

reduction_same_degree(G, nu_p, o1)