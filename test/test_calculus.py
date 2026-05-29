
from scipy.integrate import cumulative_trapezoid
import numpy as np

y = [0, 4, 16, 36, 64, 100]

cumint = cumulative_trapezoid(y, initial=0)
print(cumint)

cumint = np.cumsum(y)
print(cumint)
