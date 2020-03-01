from Code import *
import matplotlib.pyplot as plt
from scipy import optimize

lap_time_simulation(track="Test_short", formula_data="Svarog data")
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

