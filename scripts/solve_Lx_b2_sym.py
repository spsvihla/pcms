import sympy as sp
import pickle

# Load d from a pkl file
# Ensure 'coefficients.pkl' contains a dictionary mapping integers 0 to k to expressions
with open('coefficients.pkl', 'rb') as f:
    d = pickle.load(f)

k = 4
gamma, Theta_n, b_0, ln2 = sp.symbols('gamma Theta_n b_0 ln_2')
zeta = {i: sp.Symbol(f'zeta_{i}') for i in range(2, k + 6)}

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
M = (B_full - I_minus).T
M_prime = M[0:k+1, 1:k+2]

C_sol = M_prime.LUsolve(-c)
C_sol_simplified = sp.simplify(C_sol)

for i in range(k + 1):
    print(f"C_{i+1} =")
    sp.pprint(C_sol_simplified[i])
    print("-" * 40)