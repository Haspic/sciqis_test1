
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.integrate as integrate
from math import factorial, sqrt, pi, exp


class GeneralizedLaguerrePolynomial:

    def __init_cache__(self, a, x, k):
        self.cache = np.zeros(k+1)

        self.cache[0] = 1
        if k > 0:
            self.cache[1] = 1 + a - x

    def __call__(self, a, x, k):
        """Init or reset cache for given a and x value"""
        self.__init_cache__(a, x, k)
        return self.polynomial(a, x, k)

    def polynomial(self, a, x, k):

        # if k < 2:
        #     return self.cache[k]

        for i in range(2, k+1):
            self.cache[i] = self._k_evaluate(a, x, i,
                                             self.cache[i-1],
                                             self.cache[i-2])

        return self.cache

    @staticmethod
    def _k_evaluate(a, x, k, L_1, L_2):
        return ((2*(k-1) + 1 + a - x)*L_1 - ((k-1) + a)*L_2) / k


def WignerFunction(x, p, m, n, RETURN="complex"):
    L = GeneralizedLaguerrePolynomial()

    p1 = 1/pi * exp(-x**2 - p**2)
    p2 = (-1)**n * (x -1j*p)**(m-n)
    p3 = sqrt(2**(m-n) * factorial(n)/factorial(m))
    p4 = L(a=m-n, x=2*x**2 + 2*p**2, k=n)[-1]

    match RETURN:

        case "complex":
            return p1 * p2 * p3 * p4

        case "real":
            return (p1 * p2 * p3 * p4).real

        case "imag":
            return (p1 * p2 * p3 * p4).imag

def WignerFunction_to_DensityMatrix(W, m, n):

    real_func = lambda x, p: WignerFunction(x, p, m, n, "real") * W(x, p)
    imag_func = lambda x, p: WignerFunction(x, p, m, n, "imag") * W(x, p)

    real_integral = integrate.dblquad(real_func, -np.inf, np.inf, -np.inf, np.inf)
    imag_integral = integrate.dblquad(imag_func, -np.inf, np.inf, -np.inf, np.inf)

    integral = real_integral[0] + 1j * imag_integral[0]

    return 2 * pi * integral

def DensityMatrix_to_WignerFunction(densityMaxtrix, x, p):

    N = len(densityMaxtrix)

    final_sum = 0

    for m in range(N):
        for n in range(N):

            final_sum += densityMaxtrix[m, n] * WignerFunction(x, p, m, n)

    return final_sum

""" ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### """

alpha = 1+1j
N = 2
density_matrix = np.zeros((N, N), dtype=np.complex128)

def W(x, p):
    return 1/pi * exp(-(x - sqrt(2)*alpha.real)**2 - (p - sqrt(2)*alpha.imag)**2)


for m in range(N):
    for n in range(N):
        print("m:", m,"n:", n)
        density_matrix[m, n] = WignerFunction_to_DensityMatrix(W, m, n)


print(density_matrix)

for X, P in ((0, 0), (1, 1), (1, 2), (2, 1), (2, 2)):
    print("X:", X, "P:", P)
    print("Real Wigner:", W(X, P))
    print("Inverted Wigner:", DensityMatrix_to_WignerFunction(density_matrix, X, P))
