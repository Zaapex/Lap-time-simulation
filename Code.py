from Functions import *
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

name_of_track = "FSAE_track"
df = pd.read_csv(name_of_track, index_col=0, low_memory=False)

initial_conditions = {"vx_initial": 0, "t_initial": 0, "a_initial": 7.4, "sigma_initial": 0, "vy_initial": 0,
                      "ang_vel_initial": 0}


df.at[0, "vx_entry"] = initial_conditions["vx_initial"]
df.at[0, "time"] = initial_conditions["t_initial"]
df.at[0, "acceleration"] = initial_conditions["a_initial"]

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


for x in range(len(df.index)):
    """"Here we calculate max velocity in given corner"""

    radius = df.loc[x, 'R']  # get radius at this segment, between two points

    max_Velocity_cal_rear = max_velocity_rear(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                          alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, W=W)

    max_Velocity_force = drive_train()[2]


    df.at[x, "vx_max"] = min(max_Velocity_cal_rear, max_Velocity_force)

    df.fillna(method="ffill", inplace=True)

for x in range(len(df.index)):

    vx_entry = df.loc[x, 'vx_entry']  # entry velocity in this section
    acc_x = df.loc[x, "acceleration"]
    radius = df.loc[x, "R"]
    vx_max = df.loc[x, "vx_max"]
    s = df.loc[x, "s"]

    # normal force on each tire
    f_nor_r_o = normal_force_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_r_i = normal_force_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_f_o = normal_force_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_f_i = normal_force_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)

    # friction force for each tire
    f_fri_r_o = f_nor_r_o * tires()[0]
    f_fri_r_i = f_nor_r_i * tires()[0]
    f_fri_f_o = f_nor_f_o * tires()[0]
    f_fri_f_i = f_nor_f_i * tires()[0]

    # some total forces on car
    F_centripental = (mass * vx_entry ** 2 / radius)
    F_drag = alpha_drag() * vx_entry ** 2
    F_nor_total = f_nor_f_i + f_nor_f_o + f_nor_r_i + f_nor_r_o


    if f_fri_r_o**2 - (F_centripental * f_nor_r_o / F_nor_total)**2 < 0:
        print(x)
    else:
        F_acc_r_o = np.sqrt(f_fri_r_o ** 2 - (F_centripental * f_nor_r_o / F_nor_total) ** 2)
        F_acc_r_i = np.sqrt(f_fri_r_i ** 2 - (F_centripental * f_nor_r_i / F_nor_total) ** 2)
        F_acc_f_o = np.sqrt(f_fri_f_o ** 2 - (F_centripental * f_nor_f_o / F_nor_total) ** 2)
        F_acc_f_i = np.sqrt(f_fri_f_i ** 2 - (F_centripental * f_nor_f_i / F_nor_total) ** 2)


    #print(F_acc_r_o, F_acc_r_i, F_acc_f_o, F_acc_f_i)

    if vx_entry < vx_max:
        if (F_nor_total * tires()[0]) ** 2 - F_centripental ** 2 > 0:
            F_acceleration_1 = min(F_acc_r_o, F_acc_r_i)*2 - F_drag
            F_limit = drive_train()[4] - F_drag
            acc = min(F_acceleration_1, F_limit) / mass
            vx_exit = np.sqrt(vx_entry ** 2 + 2 * s * acc)
            df.at[x, "vx_exit"] = vx_exit
            df.at[x + 1, "vx_entry"] = vx_exit
            df.at[x + 1, "acceleration"] = acc

        else:
            df.at[x, "vx_exit"] = df.at[x, "vx_entry"]
            df.at[x + 1, "vx_entry"] = df.at[x, "vx_entry"]
            df.at[x + 1, "acceleration"] = df.at[x, "acceleration"]

    else:
        df.at[x, "vx_exit"] = vx_max
        df.at[x + 1, "vx_entry"] = vx_max
        df.at[x + 1, "acceleration"] = 0


df.fillna(method="ffill", inplace=True)

for x in range(len(df.index) - 2, 0, -1):

    vx_entry = df.loc[x, 'vx_entry']
    vx_exit = df.loc[x, "vx_exit"]
    radius = df.loc[x, "R"]
    acc_x = df.loc[x+1, "acceleration"]

    s = df.loc[x, "s"]

    F_centripental = (mass * vx_entry ** 2 / radius)
    F_drag = alpha_drag() * vx_entry ** 2

    # normal force on each tire
    f_nor_r_o = normal_force_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                        alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_r_i = normal_force_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                        alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_f_o = normal_force_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)
    f_nor_f_i = normal_force_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                         alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=d, v=vx_entry)

    # friction force for each tire
    f_fri_r_o = f_nor_r_o * tires()[0]
    f_fri_r_i = f_nor_r_i * tires()[0]
    f_fri_f_o = f_nor_f_o * tires()[0]
    f_fri_f_i = f_nor_f_i * tires()[0]

    F_friction = min(f_fri_r_o, f_fri_r_i)*2 + min(f_fri_f_o, f_fri_f_o)*2


    if F_friction ** 2 - F_centripental ** 2 < 0:
        max_deceleration = 0
        new_v_entry = np.sqrt(vx_exit ** 2 - 2 * max_deceleration * s)
        df.at[x, "vx_entry"] = new_v_entry
        df.at[x - 1, "vx_exit"] = new_v_entry
        df.at[x, "acceleration"] = max_deceleration
    else:
        F_braking = np.sqrt(F_friction ** 2 - F_centripental ** 2)
        F_deceleretion = F_drag + F_braking
        max_deceleration = -F_deceleretion / mass

        if vx_entry > vx_exit:
            new_v_entry = np.sqrt(vx_exit ** 2 - 2 * max_deceleration * s)
            df.at[x, "vx_entry"] = new_v_entry
            df.at[x - 1, "vx_exit"] = new_v_entry
            df.at[x, "acceleration"] = max_deceleration
        else:
            pass

plt.plot(list(range(len(df.index))), df["vx_max"][:], "g", label="V max")
plt.plot(list(range(len(df.index))), df["vx_entry"][:], "b", label="V entry")
plt.plot(list(range(len(df.index))), df["acceleration"][:], "r.", label="Acceleration")

"""plt.plot(list(range(4500)), df["vx_max"][:4500], "r")
plt.plot(list(range(4500)), df["vx_entry"][:4500], "b")
plt.plot(list(range(4500)), df["acceleration"][:4500], "g")"""

plt.legend()
plt.show()

df.fillna(method="ffill", inplace=True)


for x in range(len(df.index)):
    vx_entry = df.loc[x, 'vx_entry']
    vx_exit = df.loc[x, "vx_exit"]
    s = df.loc[x, "s"]
    acc = df.loc[x, "acceleration"]

    time = abs(s / ((vx_entry + vx_exit) / 2))

    df.at[x, "time"] = time
    df.at[x, "time_total"] = df["time"].sum()

    t_total = df["time"].sum()

print(df["time"].sum())