
using Oscar
R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

I = ideal([x^5+x^3, x+y^100, y*z^2 + z^3])
f = standard_basis(I, ordering = lex(R))

I = ideal([x^5+x^3, x+y^100])
f = standard_basis(I, ordering = lex(R))
(x+y^100)^5