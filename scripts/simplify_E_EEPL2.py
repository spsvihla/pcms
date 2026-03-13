import sympy as sp
import pickle

# Define symbols
gamma, zeta_2, zeta_3, zeta_4, Theta_n, b_0, ln2 = sp.symbols('gamma zeta_2 zeta_3 zeta_4 Theta_n b_0 ln2')
ln_n = sp.symbols('ln_n')

# Expansions
E = {2: 1/(2*zeta_2), 1: (gamma*zeta_2 + zeta_3)/(zeta_2**2), 0: b_0}
D = {
    4: 1/(4*zeta_2**2),
    3: (gamma*zeta_2 + zeta_3)/zeta_2**3,
    2: -sp.Rational(9, 10)/zeta_2 + (3*gamma**2 + 4*zeta_3)/(2*zeta_2**2) + (3*gamma*zeta_3)/(zeta_2**3) + (2*zeta_3**2)/(zeta_2**4),
    1: 1 - (9*gamma + 20*zeta_3)/(5*zeta_2) + (5*gamma**3 + 20*zeta_3*gamma + 21*zeta_3 - 30*zeta_4)/(5*zeta_2**2) + (3*zeta_3*gamma**2 + 4*zeta_3**2)/(zeta_2**3) + (4*zeta_3**2*gamma)/(zeta_2**4) + (2*zeta_3**3)/(zeta_2**5),
    0: -b_0**2
}
M = {1: 1/ln2, 0: Theta_n}

# Build polynomials
poly_D = sum(D[i] * ln_n**i for i in D)
poly_E = sum(E[i] * ln_n**i for i in E)
poly_M = sum(M[i] * ln_n**i for i in M)

# Final calculation
final_expr = sp.simplify(poly_D - 2*poly_M*poly_E + poly_M**2)
poly_final = sp.Poly(final_expr, ln_n)

# Dictionary for pickle
coeff_dict = {i: poly_final.coeff_monomial(ln_n**i) for i in range(poly_final.degree() + 1)}

# Output LaTeX
print("--- LaTeX Coefficients ---")
for deg in range(poly_final.degree(), -1, -1):
    coeff = coeff_dict.get(deg, 0)
    print(f"Coefficient of ln(n)^{{{deg}}}:")
    print(f"$$ {sp.latex(sp.cancel(coeff))} $$\n")

# Save to pickle
with open('coefficients.pkl', 'wb') as f:
    pickle.dump(coeff_dict, f)

print("Coefficients saved to coefficients.pkl")