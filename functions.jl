function cmp_monomial(o1, m1, m2)
    m1 = leading_monomial(m1, ordering=o1)
    m2 = leading_monomial(m2, ordering=o1)
    return cmp(o1, m1, m2)
end
