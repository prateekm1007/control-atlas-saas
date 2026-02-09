import numpy as np

def scrambled_null(mi_func, sig_a, sig_b, n=10):
    vals = []
    for _ in range(n):
        vals.append(mi_func(sig_a, np.random.permutation(sig_b)))
    return np.mean(vals) + 2 * np.std(vals)
