from Functions import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def lap_time_simulation(track, formula_data, simulations_name):
    name_of_track = "Tracks/" + track + ".csv"
    df = pd.read_csv(name_of_track, index_col=0, low_memory=False)
    df_data = pd.read_csv(formula_data, index_col=0)

    initial_conditions = {"vx_initial": 0, "t_initial": 0, "a_initial": 0}

    df.at[0, "vx_entry"] = initial_conditions["vx_initial"]
    df.at[0, "time"] = initial_conditions["t_initial"]
    df.at[0, "acceleration"] = initial_conditions["a_initial"]

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
    KR = 1
    g = 9.81
    a = wheelbase*CG/100
    b = wheelbase*(100-CG)/100
    mass_rear = mass * b / wheelbase
    mass_front = mass * a / wheelbase
    w = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])*2
    max_rpm = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max RPM"].index[0]])
    tire_radius = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])
    gear_ratio = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Gear ratio"].index[0]])
    v_max_teo = max_rpm*math.pi*2*tire_radius/(60*gear_ratio)
    max_torque = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max torque"].index[0]])
    max_rear_wheel_torque = max_torque*gear_ratio/tire_radius*2

    # braking system parameters
    brake_bias = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Brake bias"].index[0]])

    for x in range(len(df.index)):
        """"Here we calculate max velocity in given corner"""

        # track width added to the radius -> we need the center of the car
        radius = df.loc[x, 'R'] + track_width/2  # get radius at this segment, between two points

        max_vel_fi = max_velocity_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                              CoPy=CoPy,
                                              alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF,
                                              d=track_width)

        max_vel_ri = max_velocity_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                             CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR,
                                             d=track_width)

        def vel_fi(x):
            c4 = max_vel_fi[0]
            c2 = max_vel_fi[2]
            c0 = max_vel_fi[4]
            return (c4 * x ** 4 + c2 * x ** 2 + c0)

        def vel_ri(x):
            c4 = max_vel_ri[0]
            c2 = max_vel_ri[2]
            c0 = max_vel_ri[4]
            return (c4 * x ** 4 + c2 * x ** 2 + c0)

        if radius < 150:
            root_fi = optimize.newton(vel_fi, v_max_teo, maxiter=100000)
            root_ri = optimize.newton(vel_ri, v_max_teo, maxiter=100000)

            df.at[x, "vel_fi"] = root_fi
            df.at[x, "vel_ri"] = root_ri
            df.at[x, "vx_max"] = min(root_fi, root_ri)  # (root_fi + root_ri)/2
        else:
            df.at[x, "vel_fi"] = v_max_teo
            df.at[x, "vel_ri"] = v_max_teo
            df.at[x, "vx_max"] = v_max_teo

    df.fillna(method="ffill", inplace=True)

    for x in range(len(df.index)):

        vx_entry = df.loc[x, 'vx_entry']  # entry velocity in this section
        acc_x = df.loc[x, "acceleration"]
        radius = df.loc[x, "R"] + track_width/2
        vx_max = df.loc[x, "vx_max"]
        s = df.loc[x, "s"]

        # normal force on each tire
        f_nor_r_o = normal_force_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_entry, acc=acc_x)
        f_nor_r_i = normal_force_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_entry, acc=acc_x)
        f_nor_f_o = normal_force_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_entry, acc=acc_x)
        f_nor_f_i = normal_force_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_entry, acc=acc_x)

        # friction force for each tire
        f_fri_r_o = f_nor_r_o * coef_friction
        f_fri_r_i = f_nor_r_i * coef_friction
        f_fri_f_o = f_nor_f_o * coef_friction
        f_fri_f_i = f_nor_f_i * coef_friction

        # some total forces on car
        F_centripental = (mass * vx_entry ** 2 / radius)
        F_centripental_front = (mass_front*vx_entry**2/radius)
        F_centripental_rear = (mass_rear * vx_entry ** 2 / radius)
        F_drag = alpha_drag() * vx_entry ** 2
        F_nor_total = f_nor_f_i + f_nor_f_o + f_nor_r_i + f_nor_r_o


        if f_fri_r_i**2 - (F_centripental_rear/2)**2 < 0:
            F_acc_r_i = 0
            F_acc_r_o = 0
        else:
            F_acc_r_o = np.sqrt(f_fri_r_o ** 2 - (F_centripental_rear/2) ** 2)
            F_acc_r_i = np.sqrt(f_fri_r_i ** 2 - (F_centripental_rear/2) ** 2)
            #F_acc_f_o = np.sqrt(f_fri_f_o ** 2 - (F_centripental_front/2) ** 2)
            #F_acc_f_i = np.sqrt(f_fri_f_i ** 2 - (F_centripental_front/2) ** 2)

        if vx_entry < vx_max:
            F_acceleration_1 = min(F_acc_r_o, F_acc_r_i) * 2 - F_drag
            F_limit = max_rear_wheel_torque - F_drag
            acc = min(F_acceleration_1, F_limit) / mass
            vx_exit = np.sqrt(vx_entry ** 2 + 2 * s * acc)

            if vx_exit < vx_max:
                df.at[x, "vx_exit"] = vx_exit
                df.at[x + 1, "vx_entry"] = vx_exit
                df.at[x + 1, "acceleration"] = acc

            else:
                df.at[x, "vx_exit"] = vx_max
                df.at[x + 1, "vx_entry"] = vx_max
                df.at[x + 1, "acceleration"] = acc

        else:
            df.at[x, "vx_exit"] = vx_max
            df.at[x + 1, "vx_entry"] = vx_max
            df.at[x + 1, "acceleration"] = 0


    df.fillna(method="ffill", inplace=True)

    for x in range(len(df.index) - 2, 0, -1):
        vx_entry = df.loc[x, 'vx_entry']
        vx_exit = df.loc[x, "vx_exit"]
        radius = df.loc[x, "R"] + track_width/2
        acc_x = df.loc[x+1, "acceleration"]
        vx_max = df.loc[x, "vx_max"]

        s = df.loc[x, "s"]

        F_centripental_front = (mass_front * vx_exit ** 2 / radius)
        F_centripental_rear = (mass_rear * vx_exit ** 2 / radius)
        F_drag = alpha_Cd * vx_exit ** 2

        # normal force on each tire
        f_nor_r_o = normal_force_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                            alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_exit, acc=acc_x)
        f_nor_r_i = normal_force_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                            alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_exit, acc=acc_x)
        f_nor_f_o = normal_force_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_exit, acc=acc_x)
        f_nor_f_i = normal_force_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase, CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, d=track_width, v=vx_exit, acc=acc_x)

        # friction force for each tire
        f_fri_r_o = f_nor_r_o * coef_friction
        f_fri_r_i = f_nor_r_i * coef_friction
        f_fri_f_o = f_nor_f_o * coef_friction
        f_fri_f_i = f_nor_f_i * coef_friction

        # average of inner and outer tire friction
        F_front_friction = min(f_fri_f_o, f_fri_f_i)  # (f_fri_f_o + f_fri_f_i)/2
        F_rear_friction = min(f_fri_r_o, f_fri_r_i)  # (f_fri_r_o + f_fri_r_i)/2

        if vx_entry > vx_exit:
            F_brake_r_o = np.sqrt(F_rear_friction ** 2 - (F_centripental_rear / 2) ** 2)
            F_brake_f_o = np.sqrt(F_front_friction ** 2 - (F_centripental_front / 2) ** 2)

            """if F_front_friction ** 2 - (F_centripental_front / 2) ** 2 < 0:
                print(x, "front problem", F_front_friction, F_centripental_front)

            if F_rear_friction ** 2 - (F_centripental_rear / 2) ** 2 < 0:
                print(x, "rear problem", F_rear_friction, F_centripental_rear)"""

            # added brake bias
            if F_brake_f_o*brake_bias/100 < F_brake_r_o*(100-brake_bias)/100:
                F_braking = F_brake_f_o*4
                print(x, "front braking")
            else:
                F_braking = F_brake_r_o*4
                print(x, "Rear braking")

            F_deceleretion = F_drag + F_braking
            max_deceleration = - F_deceleretion / mass

            new_v_entry = np.sqrt(vx_exit ** 2 - 2 * max_deceleration * s)
            df.at[x, "vx_entry"] = new_v_entry
            df.at[x - 1, "vx_exit"] = new_v_entry
            df.at[x, "acceleration"] = max_deceleration

    plt.plot(list(range(len(df.index))), df["acceleration"][:], "r.", label="Pospešek")
    plt.plot(list(range(len(df.index))), df["vx_max"][:], "g", label="V max")
    plt.plot(list(range(len(df.index))), df["vx_entry"][:], "b", label="V vstopna")
    plt.title("Vstopna in maksimalna hitrost ter pojemek")
    plt.xlabel("Številka odseka")
    plt.ylim(-40, 40)
    plt.ylabel("Hitrost [m/s]" "\n"
               "Pojemek [m/s²]")
    #plt.plot(list(range(len(df.index))), df["vx_entry"][:], "b", label="V entry")
    #plt.plot(list(range(len(df.index))), df["acceleration"][:], "r.", label="Acceleration")

    """plt.plot(list(range(len(df.index))), df["vx_max"][:], "g.", label="V max")
    plt.plot(list(range(len(df.index))), df["vel_fi"][:], "b.", label="V front inner")
    plt.plot(list(range(len(df.index))), df["vel_ri"][:], "r.", label="V rear inner")"""

    """plt.plot(df["R"][:], df["vx_max"][:], "g.", label="V max")
    plt.plot(df["R"][:], df["vel_fi"][:], "b.", label="V front inner")
    plt.plot(df["R"][:], df["vel_ri"][:], "r.", label="V rear inner")"""


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

    df.at[0, "pot"] = df.loc[0, "s"]
    for x in range(len(df.index) - 1):
        pot = df.loc[x, "pot"] + df.loc[x + 1, "s"]
        df.at[x + 1, "pot"] = pot

    df.to_csv("Simulations/csv/" + simulations_name)

    return df["time"].sum()


author = "alex" #input("Author: ")
name_of_sim = "rep_test"  #input("Name of Simulation: ")
notes = "no" #input("Notes: ")
select_track = "Logatec" #input("Track name: ")
selected_settings = "Svarog data" #input("Data name: ")


df = pd.read_csv(selected_settings)

f = open("Simulations/txt/" + name_of_sim + ".txt", "w+")

f.write("Author: " + str(author) + "\n")
f.write("Name of Simulation: " + str(name_of_sim) + "\n")
f.write("Notes: " + str(notes) + "\n")

for line in range(len(df.index)):
    f.write(df["Parameter"][line] + ": " + df["Value"][line] + "\n")

lap_time = lap_time_simulation(select_track, selected_settings, name_of_sim)

f.write("Total time: " + str(lap_time) + "\n")
f.close()