using Oscar
p = 2
nu_p = tropical_semiring_map(QQ, p)

#Start with the polynomial ring
Rtx, (t, x1, x2, x3) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3])
function cmp_monomial(o1, m1, m2)
    m1 = leading_monomial(m1, ordering=o1)
    m2 = leading_monomial(m2, ordering=o1)
    return cmp(o1, m1, m2)
end

function p_minus_t_reduce(f, p)
    nu_p = tropical_semiring_map(QQ, p)
    t_list = [i[1] for i in collect(exponents(f))]
    val_ca = minimum(t_list)
    new_coefficients = [Int(2^t_list[i] / p^(val_ca)) for i in 1:length(t_list)]

    new_terms = collect(terms(f))
    for i in 1:length(new_coefficients)
        new_terms[i] = new_coefficients[i] * new_terms[i]
        new_terms[i] = new_terms[i] / (gen(Rtx, 1) ^ t_list[i])
    end
    new_terms = (gen(Rtx, 1) ^ val_ca) .* new_terms
    #Bin t and replace with p
    return sum(new_terms)
end
