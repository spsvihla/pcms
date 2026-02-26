import sympy as sp

def simplify_matrix_terms():
    gamma = sp.Symbol('gamma')
    z2, z3, z4 = sp.symbols('zeta_2 zeta_3 zeta_4')

    m11 = gamma - 4 + gamma + 4*(1 - z2)
    m21 = gamma**2 - 4*gamma + 12*(z3 - 1) + 4*gamma*(z2 - 1) + gamma**2
    m22 = gamma - 3 + 3*(1 - z2) + gamma
    m31 = (gamma**3 - 4*gamma**2 + 12*gamma - 24 + 24*(1 - z4) + 
           12*gamma*(1 - z3) + 4*gamma**2*(1 - z2) + gamma**3)
    m32 = gamma**2 - 3*gamma + 6 + 6*(z3 - 1) + 3*gamma*(z2 - 1) - gamma**2
    m33 = gamma - 2 + gamma + 2*(1 - z2)

    expressions = [
        [m11, 0, 0, 0],
        [m21, m22, 0, 0],
        [m31, m32, m33]
    ]

    for row in expressions:
        simplified_row = [sp.simplify(expr) for expr in row]
        print(simplified_row)

if __name__ == "__main__":
    simplify_matrix_terms()