import pandas as pd

df = pd.DataFrame({"Parameter": pd.Series(["Basic", "Mass", "CG in y", "CG in z", "Track width", "Wheelbase", "W",
                                           "Suspension", "Coefficient of Friction",
                                           "Aerodynamics", "Air density", "Frontal area", "Coefficient of DF",
                                           "Coefficient of Drag", "CoP in y", "CoP in z",
                                           "Drivetrain", "Tire radius", "Max Power",
                                           "Max RPM", "Gear ratio", "Max torque"]),
                  "Value": pd.Series(["/", 290, 55, 0.1, 1, 1.250, 0.4, "/", 1.4, "/", 1.225, 1.2, 3.5, 1.7, 60, 0.1,
                                      "/",
                                      0.228, 80*100**3, 7600, 5.45, 90])})

selected_settings = ['C:/Users/Alex/Desktop/Sezona 19_20/Lap time simulator/Four wheels/Test data.csv']
df2 = pd.read_csv(selected_settings[0])

print(df2)
