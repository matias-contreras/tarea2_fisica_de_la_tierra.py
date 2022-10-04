
import numpy as np
from matplotlib import pyplot as plt


G = 6.67*10**(-11)


def altura(r, h, R):
    """
    encuentra la altura del modelo de montaña en un punto r
    :param r:distancia en el eje horizontal desde el origen
    :param h: altura maxima
    :param R: ancho máximo
    :return: altura topografica en r
    """
    h_r = h * np.cos(np.pi * r / (2 * R))
    return h_r


def U_cilindro(H1, R, t, rho):
    """

    :param H1: distancia desde la cara mas cercana de la placa hasta el punto de medicion de U
    :param R: radio del cilindro
    :param t: espesor del cilindro
    :param rho: densidad del cilindro
    :return: potencial producido por el cilindro en su eje a una distancia H_1
    """
    L1 = np.sqrt(H1**2 + R**2)
    H2 = H1 + t
    L2 = np.sqrt(H2**2 + R**2)

    U = -np.pi * G * rho * (H2 * L2 - H1 * L1 + H1**2 - H2**2 +
                            R**2 * np.log((t + np.sqrt(t**2 + R**2)) / R))

    return U
