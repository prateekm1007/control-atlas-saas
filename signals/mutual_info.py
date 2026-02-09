import numpy as np
from sklearn.metrics import mutual_info_score

def mutual_information(x, y, bins=20):
    x_d = np.digitize(x, np.histogram(x, bins)[1])
    y_d = np.digitize(y, np.histogram(y, bins)[1])
    return mutual_info_score(x_d, y_d)
