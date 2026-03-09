import sympy as sp

# Initialize symbols
k = 4
gamma, theta_n, b0, n = sp.symbols('gamma Theta_n b_0 n')
zeta = {i: sp.Symbol(f'zeta_{i}') for i in range(2, k + 2)}

# Define d2(n) and d1(n)
d2_n = (1 / (10 * zeta[2]**4)) * (
    10 * zeta[2]**4 - 29 * zeta[2]**3 - 20 * gamma * zeta[2]**3 - 
    10 * theta_n * zeta[2]**3 + 15 * gamma**2 * zeta[2]**2 + 
    20 * zeta[3] * zeta[2]**2 + 30 * gamma * zeta[3] * zeta[2] + 20 * zeta[3]**3
)

d1_n = (1 / (5 * zeta[2]**5)) * (
    4 * zeta[2]**5 - (14 * gamma + 20 * zeta[3]) * zeta[2]**4 + 
    (10 * theta_n + 5 * gamma**3 + 20 * gamma * zeta[3] + 21 * zeta[3] - 30 * zeta[4] - 5 * gamma**2) * zeta[2]**3 - 
    10 * gamma * zeta[3] * zeta[2]**2 - 10 * zeta[3]**2 * zeta[2] - 
    10 * theta_n * (gamma * zeta[2]**2 + zeta[3] * zeta[2]) + 
    20 * gamma * zeta[3]**2 * zeta[2] + zeta[3]**3
)

# Define vector c
c = sp.Matrix([
    theta_n**2 - 2 * theta_n * b0 - b0**2,
    d1_n,
    d2_n,
    (gamma * zeta[2] + zeta[3]) * (1 + zeta[3]) / (zeta[2]**3),
    1 / (4 * zeta[2]**2)
])

def get_B_val(k_idx, m_idx):
    if k_idx >= 2 and 0 <= m_idx <= k_idx - 2:
        s = 0
        for i in range(1, k_idx - m_idx):
            num = (gamma**(i-1)) * ((-1)**(k_idx - m_idx - 1)) * (zeta[k_idx - m_idx - i + 1] - 1) * sp.factorial(k_idx)
            den = sp.factorial(i + m_idx)
            s += num / den
        return s - (-gamma)**(k_idx - m_idx - 1)
    return 0

# Construct M' (5x5)
B_full = sp.Matrix(k + 2, k + 1, lambda i, j: get_B_val(i, j))
I_minus = sp.Matrix(k + 2, k + 1, lambda i, j: 1 if i == j + 1 else 0)
M = (B_full - I_minus).T
M_prime = M[0:k+1, 1:k+2]

# Solve M'C = -c
C_sol = M_prime.LUsolve(-c)

# Clean up results using simplification
C_sol_simplified = sp.simplify(C_sol)

# Output results
for i in range(k+1):
    print(f"C_{i+1} =")
    sp.pprint(C_sol_simplified[i])
    print("-" * 40)