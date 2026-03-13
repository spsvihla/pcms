import sympy as sp
import numpy as np
import scipy.special as sps
import scipy.constants as spc
import pickle

# Load coefficients
with open('coefficients.pkl', 'rb') as f:
    d = pickle.load(f)

k = 4
gamma, Theta_n, b_0, ln2 = sp.symbols('gamma Theta_n b_0 ln2')
zeta = {i: sp.Symbol(f'zeta_{i}') for i in range(2, k + 6)}

# Compute b0 numerically for the substitution dictionary
z2_val = sps.zeta(2)
z3_val = sps.zeta(3)
gamma_val = 0.57721566490153286060
b0_val = (gamma_val**2 / (2 * z2_val)) + (gamma_val * z3_val) / (z2_val**2) + (z3_val**2) / (z2_val**3) + 0.1

# Setup substitution dictionary
subs_dict = {
    gamma: gamma_val,
    ln2: np.log(2),
    b_0: b0_val,
    **{zeta[i]: sps.zeta(i) for i in range(2, k + 6)}
}

c = sp.Matrix([d[i] for i in range(k + 1)])

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

# Numerical Solve
M_num = M_prime.subs(subs_dict)
c_num = c.subs(subs_dict)

C_sol = M_num.LUsolve(-c_num)

# Output
print(f"Solutions for C_1 through C_{k+1} (using computed b0 = {b0_val:.6f}):")
for i in range(k + 1):
    print(f"C_{i+1} =")
    sp.pprint(sp.expand(C_sol[i]))
    print("-" * 40)