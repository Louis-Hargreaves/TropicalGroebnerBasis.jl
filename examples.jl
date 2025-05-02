using Oscar

#Define the polynomial ring
Rtx, (t, x1, x2, x3) = polynomial_ring(ZZ, [:t, :x1, :x2, :x3])

#Example 3.8
R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x^5 + x^3, x + y^100, y*z^2 + z^3])
println("3.7: ", standard_basis(I, ordering=lex(R))[1])
println("3.7: ", standard_basis(I, ordering=lex(R))[2])
println("3.7: ", standard_basis(I, ordering=lex(R))[3])
println("3.7: ", standard_basis(I, ordering=lex(R))[4])

#Example 3.14
h1 = 10*x + 5*y^2
h2 = 3*x - z^4
I = ideal([h1, h2])
println("3.10 (not reduced): ", standard_basis(I, ordering=lex(R))[1])
println("3.10 (not reduced): ", standard_basis(I, ordering=lex(R))[2])
println("3.10 (reduced): ", standard_basis(I, ordering=lex(R), complete_reduction=true)[1])
println("3.10 (reduced): ", standard_basis(I, ordering=lex(R), complete_reduction=true)[2])

#Example 3.20
#Define the polynomial ring and valuation (tropical semiring map)
R,(t,x1,x2,x3) = QQ[:t, :x1, :x2, :x3];
nu_0 = tropical_semiring_map(QQ, max);

#Define the polynomials
g1 = (1 - t^3)*x1 + (t^3 - t^4)*x3
g2 = (1 - t^3)*x2 + (t^2 - t^4)*x3

#Define the ideal and weight vectors to calculate the initial ideal
I = ideal([g1, g2])
w = [[-1, 3, 3, 3] , [-1, 1, 4, 4], [-1, 1, 4, 5],
    [-1, 2, -1, 1], [-1, 1, 2, 4], [-1, 1, 3, 5],
    [-1, 1, 1, 10], [-1, 1, -6, 3], [-1, 1, -5, 5]]

#Calculate the initial ideal for each weight vector
for i in 1:9
    println("3.20: ", standard_basis(initial(I, nu_0, w[i])))
end

#Example 4.7
R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x^5 + x^3*y  + 1, x + y^100, y*z^2 + z^3])
println("4.7: ", standard_basis(I, ordering=lex(R))[1])
println("4.7: ", standard_basis(I, ordering=lex(R))[2])
println("4.7: ", standard_basis(I, ordering=lex(R))[3])
println("4.7: ", standard_basis(I, ordering=lex(R))[4])
