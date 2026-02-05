from scipy.special import zeta

GAMMA = 0.57721566490153286060651209008240243

A0 = 1/(2*zeta(2))
B0 = (GAMMA*zeta(2) + zeta(3))/(zeta(2)**2)
b0 = GAMMA**2/(2*zeta(2)) + GAMMA*zeta(3)/(zeta(2)**2) + zeta(3)**2/zeta(2)**3 + 1/10

A1 = 1 / (2 * zeta(2))
B1 = (GAMMA * zeta(2) + zeta(3)) / zeta(2)**2

A2 = 1 / (4 * zeta(2)**2)
B2 = (GAMMA * zeta(2) + zeta(3)) / (zeta(2)**3) 
C2 = -9/(10*zeta(2)) + (3*GAMMA**2 + 4*zeta(3))/(2*zeta(2)**2) + 3*GAMMA*zeta(3)/(zeta(2)**3) + 2*zeta(3)**2/(zeta(2)**4)
D2 = 1 - (9*GAMMA + 20*zeta(3))/(5*zeta(2)) + (5*GAMMA**3 + 20*GAMMA*zeta(3) + 21*zeta(3) - 30*zeta(4))/(5*zeta(2)**2) + (3*zeta(3)*GAMMA**2 + 4*zeta(3)**2)/(zeta(2)**3) + (4*GAMMA*zeta(3)**2)/(zeta(2)**4) + (2*zeta(3)**3)/zeta(2)**5

A3 = (2*zeta(3)/zeta(2)**2 - 1/zeta(2))
B3 = 1 - (2*GAMMA + 4*zeta(3))/zeta(2) + (4*zeta(3)*GAMMA + 4*zeta(3) - 6*zeta(4))/zeta(2)**2 + 4*zeta(3)**2/zeta(2)**3