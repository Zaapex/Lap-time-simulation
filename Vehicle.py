import math


def alpha_downforce():
        """Downforce parameters"""

        air_density = 1.225
        frontal_area = 1.2
        coeficient_of_DF = 3.5
        alpha_downforce = air_density * frontal_area * coeficient_of_DF/2
        cp_y = 60  # center of pressure in x direction, in percentage between wheels
        cp_z = 0.1  # center of pressure in z direction, given in m from wheels center

        return alpha_downforce, cp_y, cp_z

def alpha_drag():
        """Drag parameters"""

        air_density = 1.225
        frontal_area = 1.2
        coefficient_of_drag = 1.7
        alpha_drag = air_density * frontal_area * coefficient_of_drag/2

        return alpha_drag

def drive_train():
        """Power, gear ratio"""

        tire_radius = 0.228
        power = 80*100**3
        max_rpm = 7600
        gear_ratio = 5.45
        v_max = max_rpm*math.pi*2*tire_radius/(60*gear_ratio)
        constant_torque = 50
        max_torque = 90
        wheel_torque = constant_torque*gear_ratio
        max_wheel_torque = max_torque*gear_ratio
        max_wheel_force_rear = max_wheel_torque/tire_radius*2  # for two wheel it is duplicated
        max_wheel_force_front = 0

        return tire_radius, power, v_max, wheel_torque, max_wheel_force_rear, max_wheel_force_front

def construction():
        """General parameters"""

        mass = 290  # vehicle + driver
        CG = 55  # center of gravity given in percentage
        l = 1.250  # distance between wheels on different axis given in m
        a = l*CG/100  # distance from front wheels to CG
        b = l*(100-CG)/100  # distance from rear wheels to CG, in m
        h = 0.1  # height of CG given in mm from plane of both axis, in m
        Izz = 10*10**6  # moment of inertia around z axis- just assumption
        d = 1.0  # distance between wheels on the same axis  given in m
        W = 0.4  # 0.3147, 0.2244

        return mass, a, b, h, CG, l, Izz, d, W

def tires():
    """Tire data"""

    cor_stiff_front = 800  # cornering stiffness of front tires
    mu = 1.4  # coefficient of friction

    return mu, cor_stiff_front
