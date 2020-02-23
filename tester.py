import numpy as np
import math
from Vehicle import *
from Functions import *

# pre equations
mass = construction()[0]
coef_friction = tires()[0]
alpha_Cl = alpha_downforce()[0]
alpha_Cd = alpha_drag()
CP = alpha_downforce()[1]
CG = construction()[4]
KF = 0  # which axis is driven, front equals zero
KR = 1
g = 9.81
a = construction()[1]
b = construction()[2]
height_CG = construction()[3]
wheelbase = construction()[5]
CoPy = CP / 100 * wheelbase
CoPz = alpha_downforce()[2]
track_width = construction()[7]
w = drive_train()[0]*2
W = construction()[8]
d = construction()[7]


radius = 8.125
vx_entry = 9.669

k = normal_force_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)/normal_force_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)

print(k)