#!/usr/bin/env python3

## @file plot.py
#  @brief Functions for plotting sparse phylogenetic covariance matrices.
#  @author Sean Svihla
#

import numpy as np
import scipy
import matplotlib.pyplot as plt


## @brief Plot sparsity pattern of the dense matrix by downsampling.
#
#  Note that this function does not check that the decimation factor is
#  sufficient to store the decimated dense matrix in memory.
#
#  @param scipy.sparse.csr_matrix   C           The sparse phylogenetic 
#                                                   covariance matrix.
#  @param scipy.sparse.csr_matrix   Q           Wavelet basis matrix.
#  @param int                       decimation  Factor for downsampling.
#  @param str                       cmap        Name of a color map. 
#
def spy_dense(C, Q, decimation, cmap='hot'):
    n, m = C.shape
    T = scipy.sparse.csc_matrix(
            (np.ones((m,)), np.arange(m), 
                np.r_[np.arange(0, m, decimation), m]), 
            (m, (m-1)//decimation + 1)
        )
    S = scipy.sparse.csr_matrix(
            (np.ones((n,)), np.arange(n),
                np.r_[np.arange(0, n, decimation), n]),
            ((n-1)//decimation + 1, n)
        )
    QT = Q.transpose()
    CD = S @ Q @ C @ QT @ T
    plt.imshow(CD.todense(), cmap=cmap, interpolation='nearest')
    plt.colorbar()
