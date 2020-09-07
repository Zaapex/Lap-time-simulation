import matplotlib.pyplot as plt
from scipy import integrate
from Functions import *
from Presentation import *
from scipy.interpolate import splprep, splev
import matplotlib.table
from scipy import optimize

df_data = pd.read_csv("Testiranja/2020-07-08_Logatec/autox_trc_on_all_runs_test.csv", index_col=0)
df_data_test = pd.read_csv("Simulations/csv/rep_test", index_col=0)



def change_data_columns(datoteka):
    df_data = pd.read_csv(datoteka)

    colum = list()
    for col in df_data.columns:
        colum.append(col)

    zanima = colum[1::3]

    for ele in colum[len(zanima):]:
        del df_data[ele]

    df_data.columns = zanima
    df_data = df_data[::-1]

    df_data.to_csv(datoteka[:-4:]+"_test.csv")
    df_data = pd.read_csv(datoteka[:-4:]+"_test.csv")

    del df_data["Unnamed: 0"]
    df_data.to_csv(datoteka[:-4:]+"_test.csv")

font = {'size'   : 24}

matplotlib.rc('font', **font)

variables = ['Acc_raw_1', 'Acc_raw_2', 'Acc_out_1', 'Acc_out_2', 'Brk_front', 'Brk_rear', 'TractTrq', 'RegTrq', 'PwrTrq',
             'StrAngl[deg]', 'TimeSinceRst[ms]', 'RPM_L', 'RPM_R', 'Gyro_X', 'Gyro_Y', 'Gyro_Z', 'Acc_X', 'Acc_Y', 'Acc_Z',
             'Vel_X', 'Vel_Y', 'Vel_Z', 'Shock_FL', 'Shock_FR', 'Shock_RL', 'Shock_RR', 'Acker_spd_RL[RPM]', 'Acker_spd_RR[RPM]',
             'Acker_spd_EST[RPM]', 'Acker_rad_FI[m]', 'Acker_rad_FO[m]', 'Acker_rad_RI[m]', 'Acker_rad_RO[m]',
             'Acker_rad_AVG[deg]', 'Acker_steer_state', 'slip RL', 'slip RR', 'RPM FL', 'RPM FR']


# Pospesevanje in zaviranje
"""plt.plot(list(range(len(df_data.index)))[120:1920], df_data['Acc_X'][0:1800], ".", label='Pospeševanje 1')
#plt.plot(list(range(len(df_data.index)))[:1800], df_data['Acc_Y'][:1800], ".", label='Acc_Y_1')
plt.plot(list(range(len(df_data.index)))[0:1800], df_data['Acc_X'][1800:3600], ".", label='Pospeševanje 2')
#plt.plot(list(range(len(df_data.index)))[:3800], df_data['Acc_X'][:3800], ".", label='Pospševanje 1')
#plt.plot(list(range(len(df_data.index)))[:1800], df_data['Acc_X'][1800:3600], ".", label='Pospeševanje 2')
plt.title("Pospeševanje in zaviranje")
plt.xlabel("Številka vzorca")
plt.ylabel("Pospešek [m/s²]")
plt.grid()
plt.legend(fontsize="medium")
plt.show()"""

# Autocross
#plt.plot(list(range(len(df_data.index)))[:], df_data['Brk pedal'][:], ".", label='Brk pedal')
#plt.plot(list(range(len(df_data.index)))[:], df_data['Brk_rear'][:], ".", label='Brk_rear')
#plt.plot(list(range(len(df_data.index)))[5800:7800], df_data['Acc_X'][5800:7800], ".", label='Acc_X_2')
#plt.plot(list(range(len(df_data.index)))[7800:9800], df_data['Acc_X'][7800:9800], ".", label='Acc_X_3')
#plt.plot(list(range(len(df_data.index)))[9800:11800], df_data['Acc_X'][9800:11800], ".", label='Acc_X_4')
#plt.plot(list(range(len(df_data.index)))[11800:13800], df_data['Acc_X'][11800:13800], ".", label='Acc_X_5')
#plt.legend()
#plt.show()

# RPM
# tole ni prav, prestavno razmerje?
"""for x in range(len(df_data.index)):

    rpm_1 = df_data.loc[x, 'RPM FL']*math.pi*2*0.228/(60*5.45)
    rpm_2 = df_data.loc[x, 'RPM FR']*math.pi*2*0.228/(60*5.45)
    df_data.at[x, "RPM vel"] = (rpm_1+rpm_2)/2"""

# vel x in vel y srednja vrednost
"""for x in range(len(df_data.index)):

    speed = (df_data.loc[x, 'Vel_X'] + df_data.loc[x, 'Vel_Y'])/2
    df_data.at[x, "speed"] = speed"""

# hitrost iz kotnega pospeška
"""for x in range(len(df_data.index)):

    radij = 6.1
    v = np.sqrt(abs(df_data.loc[x,'Acc_Y'])*radij)
    df_data.at[x, "vel_a"] = v"""

# pot
"""df_data_test.at[0, "pot"] = df_data_test.loc[0, "s"]
for x in range(len(df_data_test.index)-1):
    pot = df_data_test.loc[x, "pot"] + df_data_test.loc[x+1, "s"]
    df_data_test.at[x+1, "pot"] = pot"""

#df_data_test.to_csv("Testiranja/2020-07-08_Logatec/autox_sim")
"""df_data.at[5800, "pot"] = 0
for x in range(2000):
    pot = df_data.loc[x + 5800, "pot"] + 0.01*df_data.loc[x + 5800, "Vel_X"]
    print(pot, x)
    df_data.at[x+5801, "pot"] = pot"""

# Pospeški Autocross
"""plt.plot(np.linspace(0, 10000, 2000), df_data['Vel_X'][5800:7800], ".", label='Vel_X')
plt.plot(np.linspace(0, 10000, 2000), df_data["RPM vel"][5800:7800], ".", label="RPM vel")
plt.plot(list(range(len(df_data_test.index)))[::], df_data_test["vx_entry"][::], ".", label="vx_entry")
plt.plot(df_data['pot'][5800:7800], df_data['Vel_X'][5800:7800], ".", label="Vel_X")
plt.plot(df_data_test['pot'][:], df_data_test['vx_entry'][:], ".", label="vx_entry")
plt.plot(df_data['pot'][5800:7800], df_data['Acc_X'][5800:7800], ".", label="Pospešek meritev")
plt.plot(df_data_test['pot'][:], df_data_test['acceleration'][:], ".", label="Pospešek simulacija")
plt.grid()
plt.title("Primerjava pospeškov pri Autocross")
plt.xlabel("Pot [m]")
plt.ylabel("Pospešek [m/s]")
plt.legend()
plt.show()"""

#df_data.to_csv("Testiranja/2020-07-08_Logatec/autox_trc_on_all_runs_test.csv")

# Pospeševanje zaviranje primerjava
"""plt.plot(np.linspace(0, 500, 330), df_data['Acc_X'][1870:2200], ".", label='Acc_X')
#plt.plot(np.linspace(0, 5000, 300), df_data['Vel_X'][1900:2200], ".", label='Vel_X')
plt.plot(np.linspace(0, 500, 500), df_data_test['acceleration'][:500], ".", label='Acceleration')
#plt.plot(np.linspace(0, 5000, 501), df_data_test['vx_entry'][:], ".", label='Vx_entry')
plt.title("Pospeševanje in zaviranje primerjava pospeškov")
plt.xlabel("Številka vzorca")
plt.ylabel("Pospešeke [m/s²]")
plt.grid()
plt.legend()
plt.show()"""

# Vožnja v krogu
# basic parameters
"""df_data = pd.read_csv("Svarog data", index_col=0)
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
w = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])*2
max_rpm = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max RPM"].index[0]])
tire_radius = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])
gear_ratio = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Gear ratio"].index[0]])
v_max_teo = max_rpm*math.pi*2*tire_radius/(60*gear_ratio)
max_torque = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max torque"].index[0]])
max_rear_wheel_torque = max_torque*gear_ratio/tire_radius*2

radius = 6.1

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
    max_vel = (root_fi + root_ri)/2
    print(root_fi, root_ri, max_vel)"""

# možno odstopanje ker vzamem max hitrosti samo na prednji gumah, najbolje bi bilo da bi vzel od vseh štirih
# radij 8.75 m- skidpad je max hitrost 8.43 m/s
# radij 6.1 m - najmanjši radij je 7.01 m/s

#plt.plot(list(range(len(df_data.index)))[:2000], df_data['RPM FR'][:2000], ".", label='RPM FR')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['RPM FL'][:2000], ".", label='RPM FL')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['RPM_L'][:2000], ".", label='RPM_L')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['RPM_R'][:2000], ".", label='RPM_R')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['RPM vel'][:2000], ".", label='RPM vel')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['Acc_X'][:2000], ".", label='Acc_X')
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['Acc_Y'][:2000], ".", label='Centripentalni pospešek')
#plt.axhline(df_data['Acc_Y'][1000:2000].mean(), color="red", label="Povprečna vrednost pospeška")
#plt.plot(list(range(len(df_data.index)))[:2000], df_data["vel_a"][:2000], ".", label='Hitrost')
#plt.axhline(df_data["vel_a"][1000:2000].mean(), color="red", label="Povprečna vrednost hitrosti")
#plt.plot(list(range(len(df_data.index)))[:2000], df_data['speed'][:2000], ".", label='Speed')
#plt.grid()
#plt.title("Hitrost v ovinku radij 5 m")
#plt.xlabel("Številka vzorca")
#plt.ylabel("Hitrost [m/s]")
#plt.xlim(1000, 2000)
#plt.ylim(0, 16)
#plt.legend()
#plt.show()

# zelo visok bočni pospešek????
"""t = 5.2
s = 2*math.pi*8.75
u = s/t
print(u)

a = u**2/8.75
print(a)"""

# koeficient DF

"""for x in range(0, len(df_data.index)):
    kf = 62
    kr = 62
    force_FL = (1.747 - df_data.loc[x, 'Shock_FL'])*kf*25.4
    force_FR = (1.57 - df_data.loc[x, 'Shock_FR'])*kf*25.4
    force_RL = (1.334 - df_data.loc[x, 'Shock_RL']) * kr*25.4
    force_RR = (1.281 - df_data.loc[x, 'Shock_RR']) * kr*25.4

    df_data.at[x, 'force_FL'] = force_FL
    df_data.at[x, 'force_FR'] = force_FR
    df_data.at[x, 'force_RL'] = force_RL
    df_data.at[x, 'force_RR'] = force_RR
    df_data.at[x, "df_total"] = force_FR+force_RR+force_FL+force_RL

    a = 1.2
    density = 1.225
    hitrost = df_data.loc[x, 'Vel_X']
    hitrost_2 = df_data.loc[x, 'RPM_L']*math.pi*2*0.228/(60*5.45)
    df_data.at[x, "hitrost"] = hitrost_2

    if hitrost_2 == 0:
        cl = 0
    else:
        cl = 2 * (force_RL + force_FL + force_FR + force_RR) / (a*density*hitrost_2**2)
    df_data.at[x, "cl"] = cl"""

# aerodinamika test
"""x1 = 0
x2 = 2250
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['RPM_L'][x1:x2], ".", label='RPM_L')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['hitrost'][x1:x2], ".", label='Hitrost')
plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['Shock_FL'][x1:x2], ".", label='Shock_FL')
plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['Shock_FR'][x1:x2], ".", label='Shock_FR')
plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['Shock_RR'][x1:x2], ".", label='Shock_RR')
plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['Shock_RL'][x1:x2], ".", label='Shock_RL')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['Vel_X'][x1:x2], ".", label='Vel_X')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['cl'][x1:x2], ".", label='CL')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['force_RR'][x1:x2], ".", label='force_RR')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['force_FL'][x1:x2], ".", label='force_FL')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['force_FR'][x1:x2], ".", label='force_FR')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['force_RL'][x1:x2], ".", label='force_RL')
#plt.plot(list(range(len(df_data.index)))[x1:x2], df_data['df_total'][x1:x2], ".", label='df total')
#plt.ylim(-5, 5)
plt.legend()
plt.show()"""

# Sara acceleration, ni skidpada
"""plt.plot(list(range(len(df_data.index)))[:2000], df_data['Acc_Y'][:2000], ".", label='Acc_Y')
plt.plot(list(range(len(df_data.index)))[:2000], df_data['Acc_X'][:2000], ".", label='Acc_X')
plt.legend()
plt.show()"""

# Autocross
plt.plot(df_data['pot'][5800:7800], (df_data['Vel_X'][5800:7800]), ".", label='Vel_X')
plt.plot(df_data_test['pot'][:], (df_data_test['vx_entry'][:]), ".", label="vx_entry")
#plt.plot(df_data['pot'][5800:7800], df_data['Acc_X'][5800:7800], ".", label="Pospešek meritev")
#plt.plot(df_data_test['pot'][:], df_data_test['acceleration'][:], ".", label="Pospešek simulacija")
plt.grid()
plt.title("Primerjava hitrosti pri Autocross")
plt.xlabel("Pot [m]")
plt.ylabel("Hitrost [m/st]")
plt.legend()
plt.show()



"""fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(df_data['pot'][5800:7800], (df_data['Vel_X'][5800:7800]), ".", label='Vel_X')
ax1.plot(df_data_test['pot'][:], (df_data_test['vx_entry'][:]), ".", label="vx_entry")
#ax2.plot(df_data['pot'][5800:7800], (df_data['Vel_X'][5800:7800])*3.6, "b.", label='Vel_X')

data_r = [0,5,10,15,20,25,30,35]
data_r_2 = [x*3.6 for x in data_r]

ax1.set_title("Primerjava hitrosti simulacij in meritev")
ax1.set_xlabel("Pot [m]")
ax1.set_ylabel("Hitrost [m/s]")
ax2.set_ylabel("Hitrost [km/h]")
ax2.set_yticks(data_r)
ax2.set_yticklabels(data_r_2)

plt.show()"""

"""for x in range(2000):
    integral = integrate.simps(df_data['RPM vel'][5800:5801 + x], list(range(len(df_data.index)))[5800:5801 + x])
    df_data.at[x+5801, "pot"] = integral/100

df_data.to_csv("Testiranja/2020-07-08_Logatec/autox_trc_on_all_runs_test.csv")"""

"""for x in range(2000):
    pot = df_data.loc[x + 5801, "pot"]
    if pot > 160.94:
        print(pot, x)"""

#cas = 1531*0.01
"""plt.plot(df_data['pot'][5800:7800], df_data['RPM vel'][5800:7800], ".", label="Hitrost meritve")
plt.plot(df_data_test['pot'][:], df_data_test['vx_exit'][:], ".", label="Hitrost simulacija")
#plt.plot(df_data['pot'][5800:7800], df_data['Vel_X'][5800:7800], ".", label="Hitrost")
plt.grid()
plt.xlabel("Pot [m]")
plt.ylabel("Hitrost [m/s]")
plt.title("Primerjava hitrosti simulacije in meritev")
plt.legend()
plt.show()"""

# k poti dodamo x in y koordinate
"""pot_k = df_data_test["pot"]
for l in range(2000):
    pot = df_data.loc[l + 5800, "pot"] #+ 0.01*df_data.loc[x + 5800, "Vel_X"]
    count = len([i for i in pot_k if i > pot])
    try:
        df_data.at[l + 5801, "x"] = df_data_test.loc[10000 - count, "x"]
        df_data.at[l + 5801, "y"] = df_data_test.loc[10000 - count, "y"]
    except:
        continue"""

#plt.plot(df_data["x"][:], df_data["y"][:], ".")
#plt.show()
#df_data.to_csv("Testiranja/2020-07-08_Logatec/autox_trc_on_all_runs_test.csv")


#plt.plot(df_data['x'][5800:7800], df_data['Vel_X'][5800:7800], ".", label="Vel_X")
#plt.plot(df_data_test['x'][:], df_data_test['y'][:], df_data_test['vx_entry'][:], ".", label="vx_entry")
"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(df_data['x'][5800:7800], df_data['y'][5800:7800], df_data['Vel_X'][5800:7800], "r.", label="Vel_X")
plt.legend()
plt.show()"""

# Primerjava časov
"""df_data.at[5800, "time_total"] = -0.9
for x in range(2000):
    df_data.at[x + 5801, "time_total"] = df_data.loc[5800 + x, "time_total"] + 0.01"""
"""df_data_test.at[0, "time_total"] = 0
for x in range(len(df_data_test.index) - 1):
    df_data_test.at[x + 1, "time_total"] = df_data_test.loc[x, "time_total"] + df_data_test.loc[x, "time"]"""

"""x = df_data["pot"][5890:5800 + 1531]
y = df_data["time_total"][5890:5800+1531]
xvals = np.linspace(0, 160.94, 10001)
yinterp = np.interp(xvals, x, y)"""

#plt.plot(xvals, yinterp, '.', label="Interpoliran čas meritve")
#plt.plot(df_data["pot"][5890:5800 + 1531], df_data["time_total"][5890:5800+1531], ".", label="Čas meritve")
#plt.plot(df_data_test["pot"][:], df_data_test["time_total"][:], ".", label="Čas simulacije")
#time_delta = np.subtract(df_data_test["time_total"][:], yinterp)
#plt.plot(df_data_test["pot"][:], time_delta, "r.", label="Delta čas")
#plt.plot(df_data_test["pot"][:], time_delta/df_data_test["time_total"][:], "r.", label="Delta čas")

"""fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(df_data_test["pot"][:], time_delta, "r.", label="Časovna razlika med simulacijo in meritvijo")
ax2.plot(df_data_test["pot"][:], time_delta/df_data_test["time_total"][:]*100, "g.", label="Procentualna časovna razlika med simulacijo in meritvijo")

ax1.set_xlabel("Pot [m]")
ax1.set_ylabel("Čas [s]")
ax2.set_ylabel("Časovno odstopanje [%]")

ax1.legend()
ax2.legend()
plt.grid()
plt.title("Primerjava skupnih časov simulacije in meritve")
plt.show()"""

#plt.plot(df_data_test["pot"][:], df_data_test["vx_max"][:], ".", label="Maksimalna hitrost")
#plt.plot(df_data_test["pot"][:], df_data_test["vx_exit"][:], ".", label="Izstopna hitrost")
#plt.plot(df_data_test["pot"][:], df_data_test["R"][:], ".", label="Radij ovinka")

"""fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(df_data_test["pot"][:], df_data_test["vx_max"][:], ".", label="Maksimalna hitrost")
ax1.plot(df_data_test["pot"][:], df_data_test["vx_exit"][:], ".", label="Izstopna hitrost")

data_r = [0,5,10,15,20,25,30,35]
data_r_2 = [x*3.6 for x in data_r]

ax1.set_title("Primerjava maksimalne in izstopne hitrosti")
ax1.set_xlabel("Pot [m]")
ax1.set_ylabel("Hitrost [m/s]")
ax2.set_ylabel("Hitrost [km/h]")
ax2.set_yticks(data_r)
ax2.set_yticklabels(data_r_2)

ax1.legend(fontsize="large")
plt.show()"""

# Mavrar grafi
"""plt.subplot(311)
plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['Vel_X'][3800:5800], ".", label='Hitrost')
plt.title("Hitrost Autocross")

plt.ylabel("Hitrost [m/s]")
plt.legend()

plt.subplot(312)
plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['Shock_FL'][3800:5800], ".", label='Shock_FL')
plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['Shock_FR'][3800:5800], ".", label='Shock_FR')
plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['Shock_RL'][3800:5800], ".", label='Shock_RL')
# plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['Shock_RR'][3800:5800], ".", label='Shock_RR')
plt.title("Pomik vzmeti Autocross")

plt.ylabel("Lega [mm]")
plt.legend()

plt.subplot(313)
plt.plot(list(range(len(df_data.index)))[3800:5800], df_data['StrAngl[deg]'][3800:5800], ".", label='StrAngl[deg]')
plt.title("Zasuk volana Autocross")

plt.ylabel("Zasuk [deg]")
plt.legend()

plt.show()"""

# Skrt histogrami hitrosti vzmeti
"""shock_FL_vel = np.gradient(df_data['Shock_FL'][3900:5800])
shock_FR_vel = np.gradient(df_data['Shock_FR'][3900:5800])
shock_RL_vel = np.gradient(df_data['Shock_RL'][3900:5800])

n_bins = np.linspace(-0.1, 0.1, 50)

plt.subplot(221)
plt.hist(shock_FL_vel, bins=n_bins, weights=np.ones(len(shock_FL_vel))/len(shock_FL_vel), label="Hitrost Shock_FL")
plt.title("Porazdelitev hitrosti Shock_FL")
plt.xlabel("Hitrost [mm/s]")
plt.ylabel("Procenti")
plt.legend()

plt.subplot(222)
plt.hist(shock_FR_vel, bins=n_bins, weights=np.ones(len(shock_FR_vel))/len(shock_FR_vel), label="Hitrost Shock_FR")
plt.title("Porazdelitev hitrosti Shock_FR")
plt.xlabel("Hitrost [mm/s]")
plt.ylabel("Procenti")
plt.legend()

plt.subplot(223)
plt.hist(shock_RL_vel, bins=n_bins, weights=np.ones(len(shock_RL_vel))/len(shock_RL_vel), label="Hitrost Shock_RL")
plt.title("Porazdelitev hitrosti Shock_RL")
plt.xlabel("Hitrost [mm/s]")
plt.ylabel("Procenti")
plt.legend()

plt.show()"""

# change braking parameters

"""plt.plot(df_data['pot'][5800:7800], df_data['Vel_X'][5800:7800], ".", label='Vel_X')
#plt.plot(df_data_test['pot'][:], df_data_test["R"][:], ".", label="R")
plt.plot(df_data_test['pot'][:], df_data_test["vx_entry"][:], ".", label="vx_entry")
#plt.plot(df_data_test['pot'][:], df_data_test["vx_max"][:], ".", label="vx_max")
#plt.plot(df_data['pot'][5800:7800], df_data['Acc_X'][5800:7800], ".", label="Pospešek meritev")
#plt.plot(df_data_test['pot'][:], df_data_test['acceleration'][:], ".", label="Pospešek simulacija")
#plt.plot(df_data_test['pot'][:], df_data_test['vx_max'][:], ".", label="vx_max")
plt.grid()
plt.title("Primerjava hitrosti pri Autocross")
plt.xlabel("Pot [m]")
plt.ylabel("Hitrost [m/s]")
plt.ylim()
plt.legend()
plt.show()"""

"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(df_data_test["x"], df_data_test["y"], df_data_test["vx_entry"], "r.")
plt.show()"""

"""plt.plot(df_data["Tekmovanje"][2::7], df_data["Točke"][2::7], "o", label="Design report")
plt.ylim(0, 150)
plt.title("Design report rezultati")
plt.xlabel("Tekmovanje")
plt.ylabel("Dosežene točke")
plt.legend()
plt.show()
"""

# novi grafi za max hitrost v ovinku za radij
"""df_data = pd.read_csv("Svarog data", index_col=0)
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
w = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])*2
max_rpm = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max RPM"].index[0]])
tire_radius = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Tire radius"].index[0]])
gear_ratio = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Gear ratio"].index[0]])
v_max_teo = max_rpm*math.pi*2*tire_radius/(60*gear_ratio)
max_torque = float(df_data["Value"][df_data.loc[df_data['Parameter'] == "Max torque"].index[0]])
max_rear_wheel_torque = max_torque*gear_ratio/tire_radius*2

radiuses = (5, 10, 15, 20, 30, 40, 60, 100, 150, 200, 400)

for radius in radiuses:

    max_vel_fi = max_velocity_front_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                              CoPy=CoPy,
                                              alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF,
                                              d=track_width)

    max_vel_fo = max_velocity_front_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                              CoPy=CoPy,
                                              alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KF,
                                              d=track_width)

    max_vel_ri = max_velocity_rear_inner(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                             CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR,
                                             d=track_width)

    max_vel_ro = max_velocity_rear_outer(a, b, m=mass, g=g, h=height_CG, w=w, alfa_cl=alpha_Cl, l=wheelbase,
                                             CoPy=CoPy,
                                             alfa_cd=alpha_Cd, CoPz=CoPz, r=radius, mu=coef_friction, K=KR,
                                             d=track_width)

    x = np.linspace(0, 80, 150)
    v = x*3.6
    print(max_vel_ri)
    plt.plot(x, max_vel_fi[0] * x ** 4 + max_vel_fi[2] * x ** 2 + max_vel_fi[4], label="Polmer ovinka=" + str(radius))


fi = list()
fo = list()
ri = list()
ro = list()
plt.legend()
plt.ylabel("F_pospešek [N]")
plt.xlabel("Hitrost [m/s]")
plt.yscale('symlog', linthreshy=0.01)
plt.title("Sila pospeška, prednja notranja guma")
plt.xlim(0, 85)
plt.ylim(0, 10**7)
plt.show()"""

# novi grafi max hitrosti sim
"""plt.subplot(311)
plt.plot(list(range(len(df_data_test.index)))[:], df_data_test['vx_entry'][:],'b', label="V vstopna")
plt.ylabel("Hitrost [m/s]")
plt.xlabel("Procent proge [%]")
plt.ylim(0, 35)

plt.subplot(312)
plt.plot(list(range(len(df_data_test.index)))[:], df_data_test['vx_max'][:],'g', label="V max")
plt.ylabel("Hitrost [m/s]")
plt.xlabel("Procent proge [%]")
plt.ylim(0, 35)

plt.subplot(313)
plt.plot(list(range(len(df_data_test.index)))[:], df_data_test['acceleration'][:],'r.', label="Pospešek")
plt.ylabel("Pospešek [m/s²]")
plt.xlabel("Procent proge [%]")

plt.show()"""

# prikaz autocross steze
"""plt.plot(df_data_test["x"][:], df_data_test["y"][:], "o")
plt.xlabel("x smer [m]")
plt.ylabel("y smer [m]")
plt.title("Autocross steza")
plt.show()"""