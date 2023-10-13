load("hecke.sage")

# 3 strands
H = BraidHeckeAlgebra(3, sqrt(-3))
print("We consider the", H)
print("specifically, powers of [1, 2]")
h = H.from_braid([1, 2])
print("It is", h^7 - h^4 == h^4 - h, "that h^7 - h^4 == h^4 - h")
print("This means that the coefficients of h^(3k+1) are h + k * (h^4 - h)")
print("This leads to the homfly polynomial:")
var('k')
h_diff = h^4 - h
coeffs = [h[i] + k * h_diff[i] for i in range(H.dimension())]
print(H(coeffs, (3*k + 1) * h.degree).homfly_polynomial())
print()

# 4 strands
H = BraidHeckeAlgebra(4, sqrt(-2))
print("Next, we consider the", H)
print("specifically, powers of [1, 2, 3]")
h = H.from_braid([1, 2, 3])
print("It is", h^17 - h^9 == h^9 - h, "that h^17 - h^9 == h^9 - h")
print("This means that the coefficients of h^(8k+1) are h + k * (h^9 - h)")
print("This leads to the homfly polynomial:")
var('k')
h_diff = h^9 - h
coeffs = [h[i] + k * h_diff[i] for i in range(H.dimension())]
print(H(coeffs, (8*k + 1) * h.degree).homfly_polynomial())
print()

# 6 strands
H = BraidHeckeAlgebra(6, sqrt(-1))
print("Finally, we consider the", H)
print("specifically, powers of [1, 2, 3, 4, 5]")
h = H.from_braid([1, 2, 3, 4, 5])
h37 = h^37
h73 = h37^2 * h.inverse()
print("It is", h73 - h37 == h37 - h, "that h^73 - h^37 == h^37 - h")
print("This means that the coefficients of h^(36k+1) are h + k * (h^37 - h)")
print("This leads to the homfly polynomial:")
var('k')
h_diff = h37 - h
coeffs = [h[i] + k * h_diff[i] for i in range(H.dimension())]
print(H(coeffs, (36*k + 1) * h.degree).homfly_polynomial())
print()

# Outputs (after ~30 minutes):
# We consider Braid Hecke algebra on 3 strands with parameter sqrt(-3) specifically powers of [1, 2]
# It is True that h^7 - h^4 == h^4 - h
# This means that the coefficients of h^(3k+1) are h + k * (h^4 - h)
# This leads to the homfly polynomial:
# (k*v^4 + k*v^2 + k + 1)*v^(6*k)

# We consider Braid Hecke algebra on 4 strands with parameter sqrt(-2) specifically powers of [1, 2, 3]
# It is True that h^9 - h^5 == h^5 - h
# This means that the coefficients of h^(8k+1) are h + k * (h^9 - h)
# This leads to the homfly polynomial:
# (2*k*v^6 + 2*k*v^4 + 2*k*v^2 + 2*k + 1)*v^(24*k)

# We consider Braid Hecke algebra on 6 strands with parameter I specifically powers of [1, 2, 3, 4, 5]
# It is True that h^73 - h^37 == h^37 - h
# This means that the coefficients of h^(36k+1) are h + k * (h^37 - h)
# This leads to the homfly polynomial:
# (6*k*v^10 + 6*k*v^8 + 6*k*v^6 + 6*k*v^4 + 6*k*v^2 + 6*k + 1)*v^(180*k)
