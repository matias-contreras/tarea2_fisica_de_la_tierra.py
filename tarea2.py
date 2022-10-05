import numpy as np
import matplotlib.pyplot as plt


G = 6.67*10**(-11)


def altura(r, h, R):
    """
    encuentra la altura del modelo de montaña en un punto r
    :param r:distancia en el eje horizontal desde el origen [kilometros]
    :param h: altura maxima [kilometros]
    :param R: ancho máximo [kilometros]
    :return: altura topografica en r [kilometros]
    """
    h_r = h * np.cos(np.pi * r / (2 * R))
    return h_r

def radios(h, h_max, R):
    """
    es la inversa de la funcion altura(r, h, R)
    :param h: altura [Km]
    :param h_max: altura maxima[Km]
    :param R: radio maximo[Km]
    :return: entrega el radio correspondiente a cierta altura[Km]
    """
    r_h = 2 * R * np.arccos(h / h_max) / np.pi
    return r_h

def raiz(h, porcentaje, rho_c=2700, rho_m=3300):
    """
    encuentra el tamaño de la raiz
    :param h: altura
    :param porcentaje: porcentaje de isostacia completa
    :param rho_c: densidad corteza
    :param rho_m: densidad manto
    :return: profundidad de la raiz
    """
    raiz = rho_c * h * porcentaje/100 / (rho_m - rho_c)
    return raiz


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

def U_montaña(num_cilindros, h_topo, R, h_corteza,
              porcentaje_raiz, rho_c = 2700, rho_m=3300):
    """
    :param num_cilindros: cantidad de cilindros usados para el calculo
    :param h_topo: altitud montaña
    :param R: radio maximo montaña
    :param h_corteza: espesor corteza
    :param porcentaje_raiz: porcentaje de raiz isostatica que se usa
    :return: potencial en el centro de la montaña
    """

    delta_rho = rho_m - rho_c
    h = np.linspace(0, h_topo, num_cilindros + 1)
    radi = radios(h, h_topo, R)
    h_raiz = raiz(h, porcentaje_raiz, rho_c, rho_m)
    h1 = h_raiz + h_corteza
    t_topo = h[1]
    t_raiz = h_raiz[1]
    u_raiz = U_cilindro(h1, radi, t_raiz, delta_rho)
    u_raiz = np.delete(u_raiz, -1)
    u_raiz = np.sum(u_raiz)
    u_topo = U_cilindro(h, radi, t_topo, rho_c)
    u_topo = np.delete(u_topo, -1)
    u_topo = np.sum(u_topo)

    u_total = - u_raiz + u_topo
    return u_total

def ondulacion(num_cilindros, h_topo, R, h_corteza,
              porcentaje_raiz, rho_c = 2700, rho_m=3300, g0 = 9.8):

    potencial = U_montaña(num_cilindros, h_topo, R, h_corteza,
              porcentaje_raiz, rho_c, rho_m)
    n = -potencial / g0
    h_raiz = raiz(h_topo, porcentaje_raiz, rho_c, rho_m)
    return [n, potencial, h_raiz]


final = ondulacion(20, 4000, 300000, 33000, 80)
print('ondulacion = ' + str(final[0]) + ' metros.')
print('potencial = ' + str(final[1]))
print('raiz = ' + str(final[2]) + ' metros')
# porcentaje de raiz entre 65% y 80% para ondulacion de 30-40 metros
