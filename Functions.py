import numpy as np
from Vehicle import *
import pandas as pd
from scipy.interpolate import splprep, splev

g = 9.81


def max_velocity_rear(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, W):
    """Calculation of max velocity for rear axle"""

    m_rear = m * b / l

    c2 = (-CoPy * W * alfa_cl * r + CoPy * alfa_cl * r + CoPz * W * alfa_cd * r + CoPz * alfa_cd * r -
          W * alfa_cd * h * r + W * alfa_cl * l * r + 2 * W * h * l * m_rear + W * l * m_rear * w + alfa_cd * h * r
          + alfa_cd * r * w - alfa_cl * l * r + 2 * h * l * m_rear + l * m_rear * w) / (
                     2 * l * r)

    c1 = 0

    c0 = (W * a * g * m - a * g * m) / (2 * l)

    coeff = [c2, c1, c0]
    Y = np.roots(coeff)

    return Y[1]


def normal_force_rear_outer(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, v, d,acc):
    """Calculation of normal force on rear outer tire"""

    m_rear = m * b / l

    Normal_Force = ((a * m * g + (h + w / 2) * m*acc + alfa_cl * v ** 2 * (
            l - CoPy) + alfa_cd * v ** 2 * ((w / 2) + CoPz)) / (2 * l)) + m_rear * v ** 2 / (d * r) * (
                                        (2 * h + w) / 2)

    return Normal_Force


def normal_force_rear_inner(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, v, d, acc):
    """Calculation of normal force on rear inner tire"""

    m_rear = m * b / l

    Normal_Force = (((a * m * g + (h + w / 2) * m*acc + alfa_cl * v ** 2 * (
            l - CoPy) + alfa_cd * v ** 2 * ((w / 2) + CoPz)) / (2 * l)) - m_rear * v ** 2 / (d*r) * (
                                        (2 * h + w) / 2))

    return Normal_Force


def normal_force_front_outer(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, v, d, acc):
    """Calculation of normal force on front outer tire"""

    m_front = m * a / l

    normal_force = ((b * m * g - (h + w / 2) * m*acc + alfa_cl * v ** 2 * (
            l - CoPy) - alfa_cd * v ** 2 * ((w / 2) + CoPz)) / (2 * l)) + m_front * v ** 2 / (d*r) * (
                                         (2 * h + w) / 2)

    return normal_force


def normal_force_front_inner(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, v, d, acc):
    """Calculation of normal force on front inner tire"""

    m_front = m * a / l

    normal_force = ((b * m * g - (h + w / 2) * m*acc + alfa_cl * v ** 2 * (
            l - CoPy) - alfa_cd * v ** 2 * ((w / 2) + CoPz)) / (2 * l)) - m_front * v ** 2 / (d * r) * (
                           (2 * h + w) / 2)

    return normal_force


def add_new_points(name_of_track):
    """With approximation add extra points in track"""

    df = pd.read_csv(name_of_track, index_col=0)

    x = df["x"].to_numpy()
    y = df["y"].to_numpy()

    # fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=2
    # is needed in order to force the spline fit to pass through all the input points.
    tck, u = splprep([x, y], s=2)  # parameter s is important for smoothing

    # evaluate the spline fits for (as many as you like) evenly spaced distance values
    xi, yi = splev(np.linspace(0, 1, 10000), tck)

    new_df = pd.DataFrame()
    new_df["x"] = xi
    new_df["y"] = yi

    new_df.to_csv(name_of_track)


def radius_of_corner(name_of_track):
    """From three points of track we get a radius
    apply this for each point of the track"""

    df = pd.read_csv(name_of_track, index_col=0)

    for x in range(len(df.index) - 2):
        A = np.array([df.loc[x, 'x'], df.loc[x, 'y']])
        B = np.array([df.loc[x + 1, 'x'], df.loc[x + 1, 'y']])
        C = np.array([df.loc[x + 2, 'x'], df.loc[x + 2, 'y']])

        AB = B - A
        BC = C - B
        AC = C - A

        a = np.linalg.norm(AB)
        b = np.linalg.norm(BC)
        c = np.linalg.norm(AC)
        d = np.sqrt(2.0 * a ** 2 * b ** 2 +
                    2.0 * b ** 2 * c ** 2 +
                    2.0 * c ** 2 * a ** 2 -
                    a ** 4 - b ** 4 - c ** 4)

        A = np.array([df.loc[x, 'x'], df.loc[x, 'y']])
        B = np.array([df.loc[x + 1, 'x'], df.loc[x + 1, 'y']])
        s = np.sqrt((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2)  # calculation of total distance in x and y direction
        df.at[x, "s"] = s

        if d == 0:
            # if their is no difference between this and previous point, radius is given as 0
            r = 0
            df.at[x, "R"] = math.floor(r)

        else:
            r = (a * b * c) / np.sqrt(2.0 * a ** 2 * b ** 2 +
                                      2.0 * b ** 2 * c ** 2 +
                                      2.0 * c ** 2 * a ** 2 -
                                      a ** 4 - b ** 4 - c ** 4)
            # calculation of the radius given three points
            df.at[x, "R"] = math.floor(r)

    df.fillna(method="ffill",
              inplace=True)  # if error occurs, pandas gives NaN value. Automatic filled with previous value

    df.to_csv(name_of_track)


def max_velocity_front_inner(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, mu, K, d):

    m_rear = m * b / l
    m_front = m * a / l

    c4 = 0.25*(1.0*CoPy**2*alfa_cl**2*d**2*mu**2*r**2 + 2.0*CoPy*CoPz*alfa_cd*alfa_cl*d**2*mu**2*r**2 +
               2.0*CoPy*alfa_cd*alfa_cl*d**2*h*mu**2*r**2 + 2.0*CoPy*alfa_cd*alfa_cl*d**2*mu**2*r**2*w -
               2.0*CoPy*alfa_cl**2*d**2*l*mu**2*r**2 + 4.0*CoPy*alfa_cl*d*h*l*m_front*mu**2*r +
               2.0*CoPy*alfa_cl*d*l*m_front*mu**2*r*w + 1.0*CoPz**2*alfa_cd**2*d**2*mu**2*r**2 +
               2.0*CoPz*alfa_cd**2*d**2*h*mu**2*r**2 + 2.0*CoPz*alfa_cd**2*d**2*mu**2*r**2*w -
               2.0*CoPz*alfa_cd*alfa_cl*d**2*l*mu**2*r**2 + 4.0*CoPz*alfa_cd*d*h*l*m_front*mu**2*r +
               2.0*CoPz*alfa_cd*d*l*m_front*mu**2*r*w - 4.0*K**2*alfa_cd**2*d**2*l**2*r**2 +
               1.0*alfa_cd**2*d**2*h**2*mu**2*r**2 + 2.0*alfa_cd**2*d**2*h*mu**2*r**2*w +
               1.0*alfa_cd**2*d**2*mu**2*r**2*w**2 - 2.0*alfa_cd*alfa_cl*d**2*h*l*mu**2*r**2 -
               .0*alfa_cd*alfa_cl*d**2*l*mu**2*r**2*w + 4.0*alfa_cd*d*h**2*l*m_front*mu**2*r +
               6.0*alfa_cd*d*h*l*m_front*mu**2*r*w + 2.0*alfa_cd*d*l*m_front*mu**2*r*w**2 +
               1.0*alfa_cl**2*d**2*l**2*mu**2*r**2 - 4.0*alfa_cl*d*h*l**2*m_front*mu**2*r -
               2.0*alfa_cl*d*l**2*m_front*mu**2*r*w - 1.0*d**2*l**2*m_front**2 + 4.0*h**2*l**2*m_front**2*mu**2 +
               4.0*h*l**2*m_front**2*mu**2*w + 1.0*l**2*m_front**2*mu**2*w**2)/(d**2*l**2*r**2)

    c2 = 0.5*(-1.0*CoPy*alfa_cl*b*d*g*m*mu**2*r - 1.0*CoPz*alfa_cd*b*d*g*m*mu**2*r - 1.0*alfa_cd*b*d*g*h*m*mu**2*r -
              1.0*alfa_cd*b*d*g*m*mu**2*r*w + 1.0*alfa_cl*b*d*g*l*m*mu**2*r - 2.0*b*g*h*l*m*m_front*mu**2 -
              1.0*b*g*l*m*m_front*mu**2*w)/(d*l**2*r)

    c0 = 0.25*b**2*g**2*m**2*mu**2/l**2

    coeff = [c4, 0, c2, 0, c0]

    return coeff


def max_velocity_front_outer(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, mu, K,d):

    m_rear = m * b / l
    m_front = m * a / l

    c4 = 0.25*(1.0*CoPy**2*alfa_cl**2*d**2*mu**2*r**2 + 2.0*CoPy*CoPz*alfa_cd*alfa_cl*d**2*mu**2*r**2 +
               2.0*CoPy*alfa_cd*alfa_cl*d**2*h*mu**2*r**2 + 2.0*CoPy*alfa_cd*alfa_cl*d**2*mu**2*r**2*w -
               2.0*CoPy*alfa_cl**2*d**2*l*mu**2*r**2 - 4.0*CoPy*alfa_cl*d*h*l*m_front*mu**2*r -
               2.0*CoPy*alfa_cl*d*l*m_front*mu**2*r*w + 1.0*CoPz**2*alfa_cd**2*d**2*mu**2*r**2 +
               2.0*CoPz*alfa_cd**2*d**2*h*mu**2*r**2 + 2.0*CoPz*alfa_cd**2*d**2*mu**2*r**2*w -
               2.0*CoPz*alfa_cd*alfa_cl*d**2*l*mu**2*r**2 - 4.0*CoPz*alfa_cd*d*h*l*m_front*mu**2*r -
               2.0*CoPz*alfa_cd*d*l*m_front*mu**2*r*w - 4.0*K**2*alfa_cd**2*d**2*l**2*r**2 +
               1.0*alfa_cd**2*d**2*h**2*mu**2*r**2 + 2.0*alfa_cd**2*d**2*h*mu**2*r**2*w +
               1.0*alfa_cd**2*d**2*mu**2*r**2*w**2 - 2.0*alfa_cd*alfa_cl*d**2*h*l*mu**2*r**2 -
               2.0*alfa_cd*alfa_cl*d**2*l*mu**2*r**2*w - 4.0*alfa_cd*d*h**2*l*m_front*mu**2*r -
               6.0*alfa_cd*d*h*l*m_front*mu**2*r*w - 2.0*alfa_cd*d*l*m_front*mu**2*r*w**2 +
               1.0*alfa_cl**2*d**2*l**2*mu**2*r**2 + 4.0*alfa_cl*d*h*l**2*m_front*mu**2*r +
               2.0*alfa_cl*d*l**2*m_front*mu**2*r*w - 1.0*d**2*l**2*m_front**2 + 4.0*h**2*l**2*m_front**2*mu**2 +
               4.0*h*l**2*m_front**2*mu**2*w + 1.0*l**2*m_front**2*mu**2*w**2)/(d**2*l**2*r**2)

    c2 = 0.5*(-1.0*CoPy*alfa_cl*b*d*g*m*mu**2*r - 1.0*CoPz*alfa_cd*b*d*g*m*mu**2*r - 1.0*alfa_cd*b*d*g*h*m*mu**2*r -
              1.0*alfa_cd*b*d*g*m*mu**2*r*w + 1.0*alfa_cl*b*d*g*l*m*mu**2*r + 2.0*b*g*h*l*m*m_front*mu**2 +
              1.0*b*g*l*m*m_front*mu**2*w)/(d*l**2*r)


    c0 = 0.25*b**2*g**2*m**2*mu**2/l**2

    coeff = [c4, 0, c2, 0, c0]

    return coeff


def max_velocity_rear_outer(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, mu, K, d):

    m_rear = m * b / l
    m_front = m * a / l

    c4 = 0.25*(1.0*CoPy**2*alfa_cl**2*d**2*mu**2*r**2 - 2.0*CoPy*CoPz*alfa_cd*alfa_cl*d**2*mu**2*r**2 +
               2.0*CoPy*alfa_cd*alfa_cl*d**2*h*mu**2*r**2 - 2.0*CoPy*alfa_cl**2*d**2*l*mu**2*r**2 -
               4.0*CoPy*alfa_cl*d*h*l*m_rear*mu**2*r - 2.0*CoPy*alfa_cl*d*l*m_rear*mu**2*r*w +
               1.0*CoPz**2*alfa_cd**2*d**2*mu**2*r**2 - 2.0*CoPz*alfa_cd**2*d**2*h*mu**2*r**2 +
               2.0*CoPz*alfa_cd*alfa_cl*d**2*l*mu**2*r**2 + 4.0*CoPz*alfa_cd*d*h*l*m_rear*mu**2*r +
               2.0*CoPz*alfa_cd*d*l*m_rear*mu**2*r*w - 4.0*K**2*alfa_cd**2*d**2*l**2*r**2 +
               1.0*alfa_cd**2*d**2*h**2*mu**2*r**2 - 2.0*alfa_cd*alfa_cl*d**2*h*l*mu**2*r**2 -
               4.0*alfa_cd*d*h**2*l*m_rear*mu**2*r - 2.0*alfa_cd*d*h*l*m_rear*mu**2*r*w +
               1.0*alfa_cl**2*d**2*l**2*mu**2*r**2 + 4.0*alfa_cl*d*h*l**2*m_rear*mu**2*r +
               2.0*alfa_cl*d*l**2*m_rear*mu**2*r*w - 1.0*d**2*l**2*m_rear**2 + 4.0*h**2*l**2*m_rear**2*mu**2 +
               4.0*h*l**2*m_rear**2*mu**2*w + 1.0*l**2*m_rear**2*mu**2*w**2)/(d**2*l**2*r**2)

    c2 = 0.5*(-1.0*CoPy*a*alfa_cl*d*g*m*mu**2*r + 1.0*CoPz*a*alfa_cd*d*g*m*mu**2*r - 1.0*a*alfa_cd*d*g*h*m*mu**2*r +
              1.0*a*alfa_cl*d*g*l*m*mu**2*r + 2.0*a*g*h*l*m*m_rear*mu**2 + 1.0*a*g*l*m*m_rear*mu**2*w)/(d*l**2*r)

    c0 = 0.25*a**2*g**2*m**2*mu**2/l**2

    coeff = [c4, 0, c2, 0, c0]

    return coeff


def max_velocity_rear_inner(a, b, m, g, h, w, alfa_cl, l, CoPy, alfa_cd, CoPz, r, mu, K, d):

    m_rear = m * b / l
    m_front = m * a / l

    c4 = 0.25*(1.0*CoPy**2*alfa_cl**2*d**2*mu**2*r**2 + 2.0*CoPy*CoPz*alfa_cd*alfa_cl*d**2*mu**2*r**2 +
               2.0*CoPy*alfa_cd*alfa_cl*d**2*h*mu**2*r**2 + 2.0*CoPy*alfa_cd*alfa_cl*d**2*mu**2*r**2*w -
               2.0*CoPy*alfa_cl**2*d**2*l*mu**2*r**2 + 4.0*CoPy*alfa_cl*d*h*l*m_rear*mu**2*r +
               2.0*CoPy*alfa_cl*d*l*m_rear*mu**2*r*w + 1.0*CoPz**2*alfa_cd**2*d**2*mu**2*r**2 +
               2.0*CoPz*alfa_cd**2*d**2*h*mu**2*r**2 + 2.0*CoPz*alfa_cd**2*d**2*mu**2*r**2*w -
               2.0*CoPz*alfa_cd*alfa_cl*d**2*l*mu**2*r**2 + 4.0*CoPz*alfa_cd*d*h*l*m_rear*mu**2*r +
               2.0*CoPz*alfa_cd*d*l*m_rear*mu**2*r*w - 4.0*K**2*alfa_cd**2*d**2*l**2*r**2 +
               1.0*alfa_cd**2*d**2*h**2*mu**2*r**2 + 2.0*alfa_cd**2*d**2*h*mu**2*r**2*w +
               1.0*alfa_cd**2*d**2*mu**2*r**2*w**2 - 2.0*alfa_cd*alfa_cl*d**2*h*l*mu**2*r**2 -
               2.0*alfa_cd*alfa_cl*d**2*l*mu**2*r**2*w + 4.0*alfa_cd*d*h**2*l*m_rear*mu**2*r +
               6.0*alfa_cd*d*h*l*m_rear*mu**2*r*w + 2.0*alfa_cd*d*l*m_rear*mu**2*r*w**2 +
               1.0*alfa_cl**2*d**2*l**2*mu**2*r**2 - 4.0*alfa_cl*d*h*l**2*m_rear*mu**2*r -
               2.0*alfa_cl*d*l**2*m_rear*mu**2*r*w - 1.0*d**2*l**2*m_rear**2 + 4.0*h**2*l**2*m_rear**2*mu**2 +
               4.0*h*l**2*m_rear**2*mu**2*w + 1.0*l**2*m_rear**2*mu**2*w**2)/(d**2*l**2*r**2)

    c2 = 0.5*(-1.0*CoPy*a*alfa_cl*d*g*m*mu**2*r - 1.0*CoPz*a*alfa_cd*d*g*m*mu**2*r - 1.0*a*alfa_cd*d*g*h*m*mu**2*r -
              1.0*a*alfa_cd*d*g*m*mu**2*r*w + 1.0*a*alfa_cl*d*g*l*m*mu**2*r - 2.0*a*g*h*l*m*m_rear*mu**2 -
              1.0*a*g*l*m*m_rear*mu**2*w)/(d*l**2*r)

    c0 = 0.25*a**2*g**2*m**2*mu**2/l**2

    coeff = [c4, 0, c2, 0, c0]

    return coeff