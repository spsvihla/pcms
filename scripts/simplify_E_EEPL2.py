import sympy as sp

gamma, z2, z3, z4, theta_n, b0 = sp.symbols('gamma zeta_2 zeta_3 zeta_4 Theta_n b_0')

# E constants
E2 = 1 / (2 * z2)
E1 = (gamma * z2 + z3) / (z2**2)
# E0 = gamma**2 / (2 * z2) + (gamma * z3) / (z2**2) + (z3**2) / (z2**3) + sp.Rational(1, 10)
E0 = b0

# D constants
D2 = (-sp.Rational(9, 10) / z2 + 
      (3 * gamma**2 + 4 * z3) / (2 * z2**2) + 
      (3 * gamma * z3) / (z2**3) + 
      (2 * z3**2) / (z2**4))

D1 = (1 - (9 * gamma + 20 * z3) / (5 * z2) + 
      (5 * gamma**3 + 20 * z3 * gamma + 21 * z3 - 30 * z4) / (5 * z2**2) + 
      (3 * z3 * gamma**2 + 4 * z3**2) / (z2**3) + 
      (4 * z3**2 * gamma) / (z2**4) + 
      (2 * z3**3) / (z2**5))

D0 = -E0**2

# Formulating Primes
D2_prime_raw = D2 - 2*E1 - 2*E2*theta_n + 1
D1_prime_raw = D1 - 2*E0 - 2*E1*theta_n + 2*theta_n
D0_prime_raw = D0 - 2*E0*theta_n + theta_n**2

# Using cancel() for best fraction simplification
D2_final = sp.cancel(D2_prime_raw)
D1_final = sp.cancel(D1_prime_raw)
D0_final = sp.cancel(D0_prime_raw)

print("--- D2 PRIME ---")
sp.pprint(D2_final)
print("\n--- D1 PRIME ---")
sp.pprint(D1_final)
print("\n--- D0 PRIME ---")
sp.pprint(D0_final)