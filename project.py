import numpy as np
import sympy as sp
from sympy import solve, solveset
from sympy.abc import x, y, z
import matplotlib.pyplot as plt
from scipy import constants
import os
from pprint import pprint


class Cavity():
    def __init__(self, lambda0, R1, R2, d, ref1, ref2, g1=0, g2=0):
        self.lambda0 = lambda0
        self.R1 = R1     #Radius of Curvature for the 1st mirror
        self.R2 = R2     #Radius of Curvature for the 2nd mirror
        self.d = d       #Length of the Cavity
        self.g1 = g1
        self.g2 = g2
        self.ref1 = ref1 #Reflectivity of the 1st mirror
        self.ref2 = ref2 #Reflectivity of the 2nd mirror


    def g_parameters(self, d, R1, R2):
        g1 = 1 - (self.d/self.R1)
        g2 = 1 - (self.d/self.R2)
        self.g1 = g1
        self.g2 = g2
        return self.g1, self.g2

    @property
    def lambda0(self):
        return self._lambda0

    @lambda0.setter
    def lambda0(self, lambda0):
        self._lambda0 = lambda0


def main():
    print("\n")
    print("--------START PROGRAM--------\n")
    print("INSERT ALL MEASUREMENTS IN MILLIMETRES!!!\n")

    while True:
        lambda0 = input("Wavelength λ: ")
        try:
            lambda0 = float(lambda0)
            if lambda0 <= 0:
                print("Please input a positive number.")
                continue
            break
        except ValueError:
            print("Please input a valid number.")

    print("\n")

    while True:                         #Check for f or R
        select = input("Do you want to input the radius of curvature [0] of the focal lengths [1]? ")
        print("\n")
        if select == '0':
            while True:
                R1 = input("Radius of Curvature: R1 = ")
                try:
                    R1 = float(R1)
                    break
                except ValueError:
                    print("Radius of Curvature has to be a number.")
                    continue
            print("\n")
            while True:
                R2 = input("Radius of Curvature: R2 = ")
                try:
                    R2 = float(R2)
                    break
                except ValueError:
                    print("Radius of Curvature has to be a number.")
                    continue
            print("\n")
        elif select == '1':
            while True:
                f1 = input("f1 = ")
                try:
                    f1 = float(f1)
                    break
                except ValueError:
                    print("Focal Length has to be a number.")
                    continue
            while True:
                f2 = input("f2 = ")
                try:
                    f2 = float(f2)
                    break
                except ValueError:
                    print("Focal Length has to be a number.")
                    continue
            R1 = 2 * f1
            R2 = 2 * f2
            break
        else:
            print("Please input either 0 or 1")
            continue
        break

    while True:
        ref1 = input("Reflectivity 1 = ")
        try:
            ref1 = float(ref1)
            if ref1 <= 0 or ref1 >= 1:
                print("Reflectivity can take values from 0 to 1.")
                continue
            break
        except ValueError:
            print("Please insert a number from 0 to 1.")

    print("\n")

    while True:
        ref2 = input("Reflectivity 2 = ")
        try:
            ref2 = float(ref2)
            if ref2 <= 0 or ref2 >= 1:
                print("Reflectivity can take values from 0 to 1.")
                continue
            break
        except ValueError:
            print("Please insert a number from 0 to 1.")

    print("\n")

    while True:                         #Check the stability of the cavity
        while True:
            d = input("Cavity Length: d: ")
            try:
                d = float(d)
                if d <= 0:
                    print("Please input a positive number.")
                    continue
                break
            except ValueError:
                print("Please input a valid number.")
        print("\n")
        cavity = Cavity(lambda0, R1, R2, d, ref1, ref2)
        g1, g2 = cavity.g_parameters(d, R1, R2)
        condition = g1*g2
        if condition <= 0 or condition >= 1:
            print(f"The cavity is not stable (g1 * g2 = {condition}). Please try again.")
        else:
            print(f"The cavity is stable (g1 * g2 = {condition}).")
            break

    print("\n")

    pprint(matrices(cavity, y))
    print("\n")
    print(f"The minimum waist is at z0 = {minimum_waist_position(cavity):.06f} mm")
    print("\n")
    z0 = minimum_waist_position(cavity)
    print("\nThe matrix after finding the y value is: ")
    pprint(matrices(cavity, z0))
    print("\n")
    w0 = minimum_waist_1(cavity)
    print(f"The first method produces w0 = {w0:.04f} mm.")
    print("\n")
    w_0, z_R = minimum_waist_2(cavity)
    print(f"The second method produces w0 = {w0:.04f} mm and Rayleigh range z_R = {z_R:.04f} mm.")
    print("\n")

    if w0 - w_0 < 1E-6:
        print("The two methods yield the same w0.")
        print("\n")
    else:
        print("The two methods yield different results.")
        print("\n")
    graphs(d, z0, z_R, w0, lambda0)

    if not os.path.exists("data.txt"):
        with open("data.txt", "w") as f:
            pass
    with open("data.txt", "w") as f:
        f.write("------Properties of the cavity------\n")
        f.write(calculations(cavity))
        f.write("\n\n")
        f.write("------Input values------\n")
        f.write(f"λ = {lambda0} mm\n")
        f.write(f"Radii of Curvature -> ROC1 = {R1} mm, ROC2 = {R2} mm\n")
        f.write(f"Reflectivities -> R1 = {ref1}, R2 = {ref2}\n")
        f.write(f"Cavity Length -> L = {d} mm\n\n")
        f.write(f"------Calculations------\n")
        f.write(f"The minimum waist is at z0 = {minimum_waist_position(cavity):.06f} mm\n")
        f.write(f"The minimum waist is w0 = {w_0:.5f} mm and the Rayleigh distance is z_R = {z_R:.3f} mm")


def matrices(cavity, y):  #Define ABCD Matrices for the Round-Trip
    M1 = sp.Matrix([[1, y], [0, 1]])
    M2 = sp.Matrix([[1, 0], [-(2/cavity.R1), 1]])
    M3 = sp.Matrix([[1, cavity.d], [0, 1]])
    M4 = sp.Matrix([[1, 0], [-(2/cavity.R2), 1]])
    M5 = sp.Matrix([[1, (cavity.d)-(y)], [0, 1]])

    M_rt = sp.Matrix(M1 * M2 * M3 * M4 * M5) #Round trip matrix
    return M_rt


def minimum_waist_position(cavity):  #Find z0 where w(z) is w0 (minimum waist)
    M_rt = matrices(cavity, y)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    z0 = solve(D - A, y)
    z0 = float(z0[0])
    return z0


def minimum_waist_1(cavity):      #First Method to find w0, this uses the solution to
    z0 = minimum_waist_position(cavity)    #1/q=(C + D(1/q))/(A + B(1/q))
    M_rt = matrices(cavity, z0)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    w0 = np.sqrt(cavity.lambda0/ np.pi) * \
         np.sqrt(abs(float(B)) / (np.sqrt(1 - ((float(A) + float(D))/2)**2)))
    return w0


def minimum_waist_2(cavity):      #q=z+izR, z=0 at the start
    z0 = minimum_waist_position(cavity)
    M_rt = matrices(cavity, z0)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    q1_2 = solve(y - (float(A) * y + float(B))/(float(C) * y + float(D)), y)
    z_R = q1_2[0] * 1j
    z_R = float(z_R)
    w_0 = np.sqrt(cavity.lambda0 * z_R / np.pi)
    return w_0, z_R


def calculations(cavity):
    c = constants.speed_of_light
    fsr = c / (2.0 * (cavity.d * 1E-3))
    finesse = np.pi * np.sqrt(np.sqrt(cavity.ref1 * cavity.ref2)) / (1 - np.sqrt(cavity.ref1 * cavity.ref2))
    fwhm = fsr / finesse
    v0 = c / cavity.lambda0
    Q = v0 / fwhm
    photon_lifetime = 1 / (2 * np.pi * fwhm)
    results = f"FSR = {fsr:.3e} Hz \nFinesse = {finesse:.3e}\nPhoton lifetime = {photon_lifetime:.3e} s\nFWHM = {fwhm:.3e}\nQ = {Q:.3e} Hz"
    return results


def graphs(d, z0, z_R, w0, lambda0):
    fig, axes = plt.subplots(3, 2, figsize=(40,50))

    #Plots of w(z) and R(z)

    z = np.linspace(0, d, 100)
    w = w0 * (1 + ((z-z0)/z_R)**2)
    axes[0][0].plot(z, w, '-')
    w_max = w0 * (1 + ((d-z0)/z_R)**2)
    axes[0][0].axis([0, d, 0, w_max*3])
    axes[0][0].set_xlabel("z (mm)", fontsize=30, color='white')
    axes[0][0].set_ylabel("w(z) (mm)", fontsize=30, color='white')
    axes[0][0].set_title("Beam Waist Inside the Cavity", fontsize=35, color='white')
    axes[0][0].tick_params(axis="both", labelsize=22, colors="white")
    axes[0][0].grid(True)

    R = (z-z0+1E-6)*(1+(z_R/(z-z0+1E-6))**2)
    axes[0][1].plot(z, R, '-')
    axes[0][1].axis([0, d, -1E4, 1E4])
    axes[0][1].set_xlabel("z (mm)", fontsize=30, color='white')
    axes[0][1].set_ylabel("R(z) (mm)", fontsize=30, color='white')
    axes[0][1].set_title("Radius of Curvature Inside the Cavity", fontsize=35, color='white')
    axes[0][1].tick_params(axis="both", labelsize=22, colors="white")
    axes[0][1].grid(True)

    #Plot intensity at the first mirror and at the minimum beam waist position

    z_0 = 0

    w_z  = w0 * np.sqrt(1 + (z_0 / z_R)**2)

    if z_0 == 0:
        R_z = np.inf
    else:
        R_z = z_0 * (1 + (z_R / z_0)**2)

    psi = np.arctan(z0 / z_R)

    x = np.linspace(-1E-3, 1E-3, 500)
    y = np.linspace(-1E-3, 1E-3, 500)
    X, Y = np.meshgrid(x, y)

    k = 2*np.pi / lambda0

    E = (w0 / w_z)\
        * np.exp(-(X**2 + Y**2) / (w_z**2))\
        * np.exp(-1j * k * z0)\
        * np.exp(+1j * psi)\
        * np.exp(-1j * k * (X**2 + Y**2) / (2 * R_z))

    I = np.abs(E)**2

    im1 = axes[1][0].contourf(X*1E3, Y*1E3, I, levels=30, cmap='jet')
    axes[1][0].set_title(f"Intensity at z = {z_0} mm", fontsize=35, color='white')
    axes[1][0].set_xlabel("x (μm)", fontsize=30, color='white')
    axes[1][0].set_ylabel("y (μm)", fontsize=30, color='white')
    axes[1][0].tick_params(axis="both", labelsize=22, colors="white")

    cbar1 = fig.colorbar(im1, ax=axes[1][0])
    cbar1.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar1.ax.axes, 'yticklabels'), color='white')

    z_0 = z0

    w_z  = w0 * np.sqrt(1 + (z_0 / z_R)**2)

    if z_0 == 0:
        R_z = np.inf
    else:
        R_z = z_0 * (1 + (z_R / z_0)**2)

    psi = np.arctan(z_0 / z_R)


    x = np.linspace(-1E-3, 1E-3, 500)
    y = np.linspace(-1E-3, 1E-3, 500)
    X, Y = np.meshgrid(x, y)

    k = 2*np.pi / lambda0

    E = (w0 / w_z)\
        * np.exp(-(X**2 + Y**2) / (w_z**2))\
        * np.exp(-1j * k * z_0)\
        * np.exp(+1j * psi)\
        * np.exp(-1j * k * (X**2 + Y**2) / (2 * R_z))

    I = np.abs(E)**2 / (2 * 377)

    im2 = axes[1][1].contourf(X*1E3, Y*1E3, I, levels=30, cmap='jet')
    axes[1][1].set_title(f"Intensity at z = {z_0} mm", fontsize=35, color='white')
    axes[1][1].set_xlabel("x (μm))", fontsize=30, color='white')
    axes[1][1].set_ylabel("y (μm)", fontsize=30, color='white')
    axes[1][1].tick_params(axis="both", labelsize=22, colors="white")

    cbar1 = fig.colorbar(im2, ax=axes[1][1])
    cbar1.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar1.ax.axes, 'yticklabels'), color='white')


    #######
    r_max = w_max
    x = np.linspace(-1E-2, 1E-2, 500)
    y = np.linspace(-1E-2, 1E-2, 500)
    X, Y = np.meshgrid(x, y)
    k = 2*np.pi / lambda0
    r = np.linspace(-r_max, r_max, 500)

    R = np.sqrt(X**2 + Y**2)
    z0 = -d/2
    z1 = d/2
    z = np.linspace(z0, z1, 500)

    R, Z = np.meshgrid(r, z)

    w_z = w0 * np.sqrt(1 + (Z / z_R)**2)
    R_z = Z * (1 + (z_R / Z)**2)

    psiZ = np.arctan(Z / z_R)

    # field
    A_r_z = (w0 / w_z) \
        * np.exp(-(R**2) / (w_z**2)) \
        * np.exp(-1j * k * Z) \
        * np.exp(+1j * psiZ) \
        * np.exp(-1j * k * (R**2) / (2 * R_z))

    c = constants.speed_of_light

    omega = 2 * np.pi * c/ lambda0

    t = 0

    E_total = A_r_z * np.exp(1j * omega * t)

    E_real = np.absolute(E_total.real)

    E = np.absolute(E_real)**2

    I = (w0 / w_z)**2 * np.exp(-2*R**2 / w_z**2)


    im3 = axes[2][0].contourf(Z, R, I, levels=30, cmap='jet')
    axes[2][0].set_title(f"$I/I_0$", fontsize=35, color='white')
    axes[2][0].set_xlabel("z (mm)", fontsize=30, color='white')
    axes[2][0].set_ylabel("r (mm)", fontsize=30, color='white')
    axes[2][0].tick_params(axis="both", labelsize=22, colors="white")

    cbar1 = fig.colorbar(im3, ax=axes[2][0])
    cbar1.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar1.ax.axes, 'yticklabels'), color='white')

    im4 = axes[2][1].contourf(Z, R, E, levels=30, cmap='jet')
    axes[2][1].set_title(f"$|Re(E(r, z, t={t}s))|^2$", fontsize=35, color='white')
    axes[2][1].set_xlabel("z (mm)", fontsize=30, color='white')
    axes[2][1].set_ylabel("r (mm)", fontsize=30, color='white')
    axes[2][1].tick_params(axis="both", labelsize=22, colors="white")

    cbar1 = fig.colorbar(im4, ax=axes[2][1])
    cbar1.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar1.ax.axes, 'yticklabels'), color='white')


    fig.set_facecolor(color='black')
    fig.savefig("figures.png", dpi=300, bbox_inches='tight')

    return fig


if __name__ == "__main__":
    main()
