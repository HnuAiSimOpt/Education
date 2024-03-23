import numpy as np

def lagrange_basis(x, x i):

    return np.prod([(x - xj) / (x i - xj) for j, xj in enumerate(x) if xj != x i], axis=0)

def lagrange_interpolation(x, y, x i):

    return sum(yi * lagrange_basis(x i, xj) for xj, yi in zip(x, y))