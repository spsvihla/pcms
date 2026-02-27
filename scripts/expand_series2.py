import sympy as sp

def get_expansion(k):
    L = sp.Symbol('log(n)')
    gamma = sp.Symbol('gamma')
    
    def zeta_sym(s):
        return sp.Symbol(f'zeta({s})')

    expansion = L**k - L**(k-1)
    
    for m in range(k - 1):
        inner_contribution = 0
        for i in range(1, k - m):
            term1_num = ((-1)**(k - m - 1)) * (zeta_sym(k - m - i + 1) - 1) * (gamma**(i - 1))
            term1_den = sp.factorial(i + m)
            term1 = sp.factorial(k) * (term1_num / term1_den)
            inner_contribution += term1
        
        term2 = -((-gamma)**(k - m - 1))
        expansion += (inner_contribution + term2) * (L**m)
        
    return expansion

def print_ordered_results(expr, k):
    L = sp.Symbol('log(n)')
    
    # Convert to a polynomial in L to force ordering
    poly_expr = sp.Poly(expr, L)
    terms = poly_expr.all_terms() # Returns list of ((pow,), coeff)
    
    # Sort terms by power descending
    sorted_terms = sorted(terms, key=lambda x: x[0][0], reverse=True)
    
    print(f"--- Expansion for k={k} (Ordered) ---")
    
    full_latex = []
    for (pow_tuple,), coeff in sorted_terms:
        p = pow_tuple
        # Pretty print to console
        print(f"Power log^{p}(n):")
        sp.pprint(coeff)
        print("-" * 20)
        
        # Build LaTeX string
        l_term = f"\\log^{{{p}}}(n)" if p > 1 else ("\\log(n)" if p == 1 else "")
        coeff_latex = sp.latex(sp.simplify(coeff))
        full_latex.append(f"\\left( {coeff_latex} \\right) {l_term}")

    print("\nFull LaTeX representation:")
    print(" + ".join(full_latex).replace("+ -", "- ") + r" + O(\log^{-1}(n))")

if __name__ == "__main__":
    k_val = 6
    result = get_expansion(k_val)
    print_ordered_results(result, k_val)