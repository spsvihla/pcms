import sympy as sp

def compute_series_expansion(k_val):
    n, gamma = sp.symbols('n gamma')
    m = sp.symbols('m', integer=True)
    i = sp.symbols('i', integer=True)
    log_n = sp.log(n)

    total_sum = 0
    
    for m_val in range(k_val):
        inner_sum = 0
        limit = k_val - m_val
        
        for i_val in range(1, limit + 1):
            term = (gamma**(i_val - 1) * (-1)**(k_val - m_val - i_val) * sp.factorial(k_val)) / \
                   sp.factorial(i_val + m_val)
            inner_sum += term
            
        total_sum += inner_sum * (log_n**m_val)

    return sp.simplify(total_sum)

k_input = 4
result = compute_series_expansion(k_input)
print(f"Expansion for k={k_input}:")
sp.pprint(result)