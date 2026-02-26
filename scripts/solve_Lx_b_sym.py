import sympy as sp

gamma, z2, z3, z4, b0, theta_n = sp.symbols('gamma zeta_2 zeta_3 zeta_4 b_0 Theta_n')

M = sp.Matrix([
    [2*gamma - 4*z2, 0, 0],
    [2*gamma**2 + 4*gamma*z2 - 8*gamma + 12*(z3 - 1), 2*gamma - 3*z2, 0],
    [2*gamma**3 - 4*gamma**2*z2 - 12*gamma*z3 + 24*gamma - 24*z4, 3*gamma*z2 - 6*gamma + 6*z3, 2*gamma - 2*z2]
])

Y = sp.Matrix([
    1 / (2*z2),
    (gamma*z2 + z3) / (z2**2) - 1,
    b0 - theta_n
])

X = M.LUsolve(Y)

A_sol = sp.simplify(X[0])
B_sol = sp.simplify(X[1])
C_sol = sp.simplify(X[2])

print(f"A = {A_sol}")
print(f"B = {B_sol}")
print(f"C = {C_sol}")