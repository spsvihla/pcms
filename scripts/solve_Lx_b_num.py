import sympy as sp

gamma_val = 0.5772156649
z2_val = 1.6449340668  # pi^2 / 6
z3_val = 1.2020569032
z4_val = 1.0823232337

gamma, z2, z3, z4, b0, theta_n = sp.symbols('gamma zeta_2 zeta_3 zeta_4 b_0 Theta_n')

M = sp.Matrix([
    [2*gamma - 4*z2, 0, 0],
    [2*gamma**2 + 4*gamma*z2 - 8*gamma + 12*(z3 - 1), 2*gamma - 3*z2, 0],
    [2*gamma**3 - 4*gamma**2*z2 - 12*gamma*z3 + 24*gamma - 24*z4, 3*gamma*z2 - 6*gamma + 6*z3, 2*gamma - 2*z2]
])

b0_val = (gamma**2 / (2*z2) + gamma*z3 / (z2**2) + z3**2 / z2**3 + 0.1).subs({
    gamma: gamma_val, z2: z2_val, z3: z3_val
})

Y = sp.Matrix([
    1 / (2*z2),
    (gamma*z2 + z3) / (z2**2) - 1,
    b0_val - theta_n
])

M_num = M.subs({gamma: gamma_val, z2: z2_val, z3: z3_val, z4: z4_val})
Y_num = Y.subs({gamma: gamma_val, z2: z2_val, z3: z3_val})

solution = M_num.LUsolve(Y_num)

A_val = solution[0].evalf()
B_val = solution[1].evalf()
C_linear = solution[2].evalf()

print(f"A = {A_val:.10f}")
print(f"B = {B_val:.10f}")
print(f"C = {C_linear}")