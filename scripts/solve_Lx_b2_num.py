import sympy as sp

# Define symbols
k = 4
gamma, theta_n, n = sp.symbols('gamma Theta_n n')
zeta = {i: sp.Symbol(f'zeta_{i}') for i in range(2, k + 2)}

# 1. Define the specific value for b0
b0 = (gamma**2 / (2 * zeta[2]) + 
      (gamma * zeta[3]) / (zeta[2]**2) + 
      (zeta[3]**2) / (zeta[2]**3) + sp.Rational(1, 10))

# 2. Define d2(n) and d1(n)
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

# 3. Define vector c
c = sp.Matrix([
    theta_n**2 - 2 * theta_n * b0 - b0**2,
    d1_n,
    d2_n,
    (gamma * zeta[2] + zeta[3]) * (1 + zeta[3]) / (zeta[2]**3),
    1 / (4 * zeta[2]**2)
])

# 4. Construct Matrix M'
def get_B_val(k_idx, m_idx):
    if k_idx >= 2 and 0 <= m_idx <= k_idx - 2:
        s = 0
        for i in range(1, k_idx - m_idx):
            num = (gamma**(i-1)) * ((-1)**(k_idx - m_idx - 1)) * (zeta[k_idx - m_idx - i + 1] - 1) * sp.factorial(k_idx)
            den = sp.factorial(i + m_idx)
            s += num / den
        return s - (-gamma)**(k_idx - m_idx - 1)
    return 0

B_full = sp.Matrix(k + 2, k + 1, lambda i, j: get_B_val(i, j))
I_minus = sp.Matrix(k + 2, k + 1, lambda i, j: 1 if i == j + 1 else 0)
M_prime = (B_full - I_minus).T[0:k+1, 1:k+2]

# 5. Numerical Substitution
# Replace these values with your actual numerical parameters
numerical_values = {
    gamma: 0.5,
    zeta[2]: 1.2,
    zeta[3]: 1.1,
    zeta[4]: 1.0,
    zeta[5]: 0.9
}

M_num = M_prime.subs(numerical_values)
c_num = c.subs(numerical_values)

# Solve for C in terms of Theta_n
# Since c_num contains Theta_n, C_sol will be a vector of polynomials/expressions in Theta_n
C_sol = M_num.LUsolve(-c_num)

print(f"Solutions for C_1 through C_5 in terms of Theta_n (with gamma={numerical_values[gamma]}):")
for i in range(k+1):
    # Expand to see the polynomial structure clearly
    expr = sp.expand(C_sol[i])
    print(f"C_{i+1} = {expr}")