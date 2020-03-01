from Code import *
import matplotlib.pyplot as plt
from scipy import optimize

#lap_time_simulation(track="Test_short", formula_data="Svarog data")
track = "Test_short"
formula_data = "Svarog data"
name_of_track = "Tracks/" + track + ".csv"
df = pd.read_csv(name_of_track, index_col=0, low_memory=False)
df_data = pd.read_csv(formula_data, index_col=0)


# basic parameters
mass = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Mass"].index[0]])
CG = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "CG in y"].index[0]])
height_CG = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "CG in z"].index[0]])
track_width = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Track width"].index[0]])
wheelbase = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Wheelbase"].index[0]])
W = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "W"].index[0]])

# suspension parameters
coef_friction = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Coefficient of Friction"].index[0]])

# aerodynamics parameters
air_density = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Air density"].index[0]])
frontal_area = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Frontal area"].index[0]])
coef_of_df = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Coefficient of DF"].index[0]])
coef_of_drag = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Coefficient of Drag"].index[0]])
alpha_Cl = air_density * frontal_area * coef_of_df/2
alpha_Cd = air_density * frontal_area * coef_of_drag/2
CP = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "CoP in y"].index[0]])
CoPy = CP / 100 * wheelbase
CoPz = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "CoP in z"].index[0]])

# drivetrain parameters
KF = 0  # which axis is driven, front equals zero
KR = 0.5
g = 9.81
a = wheelbase*CG/100
b = wheelbase*(100-CG)/100
w = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])*2
max_rpm = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max RPM"].index[0]])
tire_radius = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])
gear_ratio = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Gear ratio"].index[0]])
v_max_teo = max_rpm*math.pi*2*tire_radius/(60*gear_ratio)


"""for radius in range(1, 40, 4):
    max_vel_fi = max_velocity_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF, d=track_width)

    max_vel_fo = max_velocity_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF, d=track_width)

    max_vel_ri = max_velocity_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR, d=track_width)

    max_vel_ro = max_velocity_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR, d=track_width)

    def vel_fi(x):
        c4 = max_vel_fi[0]
        c2 = max_vel_fi[2]
        c0 = max_vel_fi[4]
        return (c4*x ** 4 + c2*x + c0)

    def vel_fo(x):
        c4 = max_vel_fo[0]
        c2 = max_vel_fo[2]
        c0 = max_vel_fo[4]
        return (c4*x ** 4 + c2*x + c0)

    def vel_ri(x):
        c4 = max_vel_ri[0]
        c2 = max_vel_ri[2]
        c0 = max_vel_ri[4]
        return (c4*x ** 4 + c2*x + c0)

    def vel_ro(x):
        c4 = max_vel_ro[0]
        c2 = max_vel_ro[2]
        c0 = max_vel_ro[4]
        return (c4*x ** 4 + c2*x + c0), [c4, c2, c0]


    root_fi = optimize.newton(vel_fi, radius, maxiter=1000)
    root_fo = optimize.newton(vel_fo, radius, maxiter=1000)
    root_ri = optimize.newton(vel_ri, radius, maxiter=1000)
    root_ro = optimize.newton(vel_ro, radius, maxiter=10000)
    plt.scatter(radius, root_fi, color="b", label="Front inner")
    plt.scatter(radius, root_fo, color="y", label="Front outer")
    plt.scatter(radius, root_ri, color="r", label="Rear inner")
    #plt.scatter(radius, root_ro, color="g", label="Rear outer")
    #print(root_ro, radius)"""

radius = 20
"""max_vel_fi = max_velocity_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF, d=track_width)

max_vel_fo = max_velocity_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF, d=track_width)

max_vel_ri = max_velocity_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR, d=track_width)

max_vel_ro = max_velocity_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR, d=track_width)
x = np.linspace(0, 20, num=100)

f_fi = []
f_fo = []
f_ri = []
f_ro = []

for i in range(len(x)):
    f_fi.append(max_vel_fi[0]*x[i]**4 + max_vel_fi[2]*x[i]**2 + max_vel_fi[4])
    f_fo.append(max_vel_fo[0] * x[i] ** 4 + max_vel_fo[2] * x[i] ** 2 + max_vel_fo[4])
    f_ri.append(max_vel_ri[0] * x[i] ** 4 + max_vel_ri[2] * x[i] ** 2 + max_vel_ri[4])
    f_ro.append(max_vel_ro[0] * x[i] ** 4 + max_vel_ro[2] * x[i] ** 2 + max_vel_ro[4])

plt.plot(x, f_fi, label="Front inner")
plt.plot(x, f_fo, label="Front outer")
plt.plot(x, f_ri, label="Rear inner")
plt.plot(x, f_ro, label="Rear outer")
plt.legend()
plt.show()"""


x = np.linspace(0, 25, num=400)
print(x)
for radius in range(5, 50, 5):
    force_outer = []
    for v in x:
        m_rear = mass * b / wheelbase
        Force_rear_outer = (((coef_friction * (((a * mass * g - (height_CG + w / 2) * alpha_Cd * v ** 2 + alpha_Cl * v ** 2 * (
            wheelbase - CoPy) + alpha_Cd * v ** 2 * ((w / 2) + CoPz)) / (2 * wheelbase)) + m_rear * v ** 2 / (track_width*radius) * (
                          (2 * height_CG + w)/2)))**2 - (0.7*m_rear * v ** 2 / radius) ** 2) ** (1 / 2)
            - KR * alpha_Cd * v ** 2)
        """Force_rear_inner = (((coef_friction * (((a * mass * g - (height_CG + w / 2) * alpha_Cd * v ** 2 + alpha_Cl * v ** 2 * (
            wheelbase - CoPy) + alpha_Cd * v ** 2 * ((w / 2) + CoPz)) / (2 * wheelbase)) - m_rear * v ** 2 / (track_width*radius) * (
                          (2 * height_CG + w)/2)))**2 - (0.475*m_rear * v ** 2 / radius) ** 2) ** (1 / 2)
            - KR * alpha_Cd * v ** 2)"""
        force_outer.append(Force_rear_outer)
    plt.plot(x, force_outer, label=radius)


plt.legend()
plt.show()