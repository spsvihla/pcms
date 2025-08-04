import numpy as np
import pcms.haar 

ys = np.linspace(-0.01, 0.01, 1000)

n_leaves = 10000
f = np.random.dirichlet(alpha=np.ones((n_leaves,), dtype=float))
g = np.random.dirichlet(alpha=np.ones((n_leaves,), dtype=float))
func = f - g

pcms.haar.cdf_proj_cbst(ys, func)