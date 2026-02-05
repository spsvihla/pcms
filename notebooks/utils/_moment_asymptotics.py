import numpy as np
from scipy.integrate import quad

import pcms.tree
from .constants import *

##
# Helper functions
##
def _integrand(x, n):
    return (1 - x**n) / (1 - x)

def _harmonic_number(n):
    hn, _ = quad(_integrand, 0, 1, args=(n,))
    return hn

##
# E[L_n]
##
ln_cache = {1: 0}
def exp_ln_exact(n):
    if n in ln_cache:
        return ln_cache[n]

    s = np.sum([exp_ln_exact(i) / (n-i) for i in range(1,n)])
    ln_cache[n] = 1 + s / _harmonic_number(n-1)
    return ln_cache[n]

def draw_ln(n_samples, n_leaves):
    depths = []
    for _ in range(n_samples):
        tree = pcms.tree.cbst(n_leaves=n_leaves, planted=False)
        _, depths_ = tree.find_leaves(return_depths=True)
        depths.append(np.random.choice(depths_))
    return np.array(depths)

def ELn(n): 
    return A0 * np.log(n)**2 + B0 * np.log(n) + b0

##
# E[EPL(Tn)]
##
def exp_epln_exact(n):
    return n * exp_ln_exact(n)

def draw_epl(n_samples, n_leaves):
    epls = []
    for _ in range(n_samples):
        tree = pcms.tree.cbst(n_leaves=n_leaves, planted=False)
        epls.append(tree.find_epl())
    return np.array(epls)

def EEPLn(n): 
    return A1 * n * np.log(n)**2 + B1 * n * np.log(n) + b0 * n

##
# E[EPL(Tn)^2]
##
zn_cache = {1: 0}
def exp_zn_exact(n):
    if n in zn_cache:
        return zn_cache[n]
    
    s1 = np.sum([(i-1)/(n-i)*exp_zn_exact(i) for i in range(1,n)])
    s2 = np.sum([exp_ln_exact(i)*exp_ln_exact(n-i) for i in range(1,n)])

    zn_cache[n] = 2*exp_ln_exact(n) - 1 + (s1+s2)/(n-1)/_harmonic_number(n-1)
    return zn_cache[n]

ln2_cache = {1: 0}
def exp_ln2_exact(n):
    if n in ln2_cache:
        return ln2_cache[n]
    
    s = np.sum([exp_ln2_exact(i)/(n-i) for i in range(1,n)])

    ln2_cache[n] = 2*exp_ln_exact(n) - 1 + s/_harmonic_number(n-1)
    return ln2_cache[n]

def exp_epln2_exact(n):
    return n*exp_ln2_exact(n) + n*(n-1)*exp_zn_exact(n)

def draw_epl2(n_samples, n_leaves):
    epl2s = draw_epl(n_samples, n_leaves)
    return np.pow(epl2s, 2)

def EEPLn2(n): 
    return A2 * n**2 * np.log(n)**4 + B2 * n**2 * np.log(n)**3 \
        + C2 * n**2 * np.log(n)**2 + D2 * n**2 * np.log(n)

##
# Var(EPL(Tn))
##
def var_epl_exact(n):
    return exp_epln2_exact(n) - exp_epln_exact(n)**2

def VEPLn(n):
    return A3 * n**2 * np.log(n)**2 + B3 * n**2 * np.log(n)