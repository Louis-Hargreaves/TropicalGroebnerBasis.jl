using Oscar

R,(z,y) = polynomial_ring(ZZ, [:z, :y], internal_ordering=:degrevlex)
p1 = z*y + z^2 + y^2 + z + y + 1
print(leading_monomial(p1))
