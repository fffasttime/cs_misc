# PCA compress and rebuild

import numpy as np

def pca(data, PCA_DIM):
    # cetter
    means = np.mean(data, axis = 0)
    data -= means
    
    data = data.T

    # # Eigen decomposition
    #cov = data @ data.T
    #diag, P = np.linalg.eig(cov)
    # # cov = P @ np.diag(diag) @ P.T , so diag = P.T @ cov @ P
    #P=P.T

    # svd is a easiser way
    U, diag, V = np.linalg.svd(data, full_matrices=False)
    P = U.T

    print("eig value: ", diag)
    P1 = P[:PCA_DIM,:]

    # compressed data
    output = P1 @ data

    # rebuild
    data1 = P1.T @ output
    data1 = data1.T
    data1 += means
    
    return data1