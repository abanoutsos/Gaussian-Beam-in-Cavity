# Gaussian Beams in a Cavity
#### Video Demo:  <URL HERE>
#### Description:
The main objective of this project is to plot the real part of the electric field, as well as the beam waist $w(z)$, the wavefront curvature $R(z)$ and the intensity at different positions inside an optical cavity. In order to achieve this, the program requires the user to input the following parameters:

- The wavelength $λ$ of the laser used in the experiment.
- The Radius of Curvature (or their focal lengths) for both mirrors.
- Their respective reflectivities.
- The length of the cavity.

All of these values are asked to be inserted in millimetres, with the exception of the reflectivities, as these are dimensionless. If the user inputs a string instead of a number, the program will handle this error and prompt them to type the value again. Moreover, the cavity has to be stable, something that is handled in the main function. If the cavity is not stable, then the user is prompted to give a different value for the cavity length.

As soon as the values given by the user satisfy the stability condition, the user can see the values of $w_0$ and $z_R$ in the terminal window. They will also be able to know if the two different methods used to calculate $w_0$ produce the same value (explained in $\textit{'The Program'}$ section).

When the program ends running, the user will find two new files, a .txt and a .png file. In the former, the user can extract important information on the cavity, such as the aforementioned $w_0$ and $z_R$, but also on the finesse of the cavity, the photon lifetime, etc. In the latter file, the user will find a picture of six graphs, each one of these depicting different physical parameters.

In order to make this project as concise as possible, the program includes a class for all the necessary parameters of the cavity and six custom functions, which handles operations such as matrix multiplication and graphs.

The libraries that I have used for this project are numpy, sympy, matplotlib and scipy.


## Theory of Gaussian Beams
#### Preliminaries
The electric field of the $\text{TEM}_{00}$ mode is given by this equation

$$
\begin{equation}
    E(r, z) = \frac{w_0}{w(z)}e^{-r^2/w^2(z)}e^{-ikr^2/2R(z)}e^{-i(kz - \text{tan}^{-1}(\frac{z}{z_R}))}
\end{equation}
$$

where

- $w_0$ is the size of the beam waist at $z_0$ (minimum beam waist)
- $w(z) = w_0 (1 + (\frac{z}{z_R})^2)^{1/2}$, the beam waist everywhere
- $R(z) = z (1 + (\frac{z_R}{z})^2)^{1/2}$ is the wavefront curvature

The intensity is given by this formula:

$$
\begin{equation}
    I = \frac{\vert{E(r, z)}\vert^2}{2η}
\end{equation}
$$

where $η = 377 Ω$ is the impedence of free space.

We introduce the parameter $q$ as

$$
\begin{equation}
    q = z + iz_R
\end{equation}
$$

The parameter $z$ has information on how far you are from the waist, while the parameter $z_R$ is the Rayleigh distance and is equal to

$$
\begin{equation}
    z_R = \frac{π w_0^{2}}{λ/n}
\end{equation}
$$

where $n$ is the refractive index of the medium. In this project, it is assumed that the cavity is in free space and there is no medium between the two mirror, so $n=1$, i.e. the refractive index of air.

Although eq.(3) is useful, in literature it is preferable to use

$$
\begin{equation}
    \frac{1}{q} = \frac{1}{R(z)} - i \frac{λ}{nπw^2(z)}
\end{equation}
$$

From ABCD Matrix theory, we know that

$$
\begin{equation}
    \frac{1}{q_2} = \frac{C + D/q_1}{A + B/q_1}
\end{equation}
$$

or, equivalently

$$
\begin{equation}
    q_2 = \frac{A/q_1 + B}{C/q_1 + D}
\end{equation}
$$

The values A, B, C, D are given by the total transfer matrix of the optical system

$$
\begin{bmatrix}
    A & B\\
    C & D
\end{bmatrix}
$$

Knowing $1/q_1$ and using eq. (5) and (6) will help us find $R(z)$ (the real part of the fraction) and $w(z)$ (the imaginary part of the fraction).


#### Mode Locking in a Cavity

Inside the cavity, we want to have the same wave "pattern", that is a stable, reproducible and defined structure. Suppose that the length of the cavity is L, and we measure the properties of the beam at a distance $y$ inside the cavity. When the cavity is mode-locked, after a round trip, the beam will look exactly the same as when it was measured previously. Using the ABCD Law, eq.(4) and (5) are now equal, thus we get

$$
\begin{equation}
    \frac{A}{q} + \frac{B}{q^2} = C + \frac{D}{q} \Longrightarrow
\end{equation}
$$

$$
\begin{equation}
    \frac{B}{q^2} + \frac{A-D}{q^2} - C = 0
\end{equation}
$$

Using the self-consisting condition

$$
\begin{equation}
    AD-BC = 1
\end{equation}
$$

and the formulas for the quadratic equation, we get the following solutions

$$
\begin{equation}
    \frac{1}{q} = -\frac{A-D}{2B} \pm i\frac{\sqrt{1-(\frac{A+D}{2})^2}}{B}
\end{equation}
$$

and using eq.(5)

$$
\begin{equation}
    -\frac{A-D}{2B} \pm i\frac{\sqrt{1-(\frac{A+D}{2})^2}}{B} = \frac{1}{R(z)} - i \frac{λ}{nπw^2(z)}
\end{equation}
$$

In order to have a correct physical derivation, the imaginary part from the left hand size has to be positive. As a result,

$$
\begin{equation}
    w(z) = \sqrt{\frac{λ}{π}} \sqrt{\frac{\vert{B}\vert}{\sqrt{1 - (\frac{A+D}{2})^2}}}
\end{equation}
$$

$$
\begin{equation}
    \frac{1}{R(z)} = \frac{D-A}{2B}
\end{equation}
$$


## Matrix Theory
The program handles matrix multiplication between different optical elements. The most important are

1. $$\begin{bmatrix}
    1 & d\\
    0 & 1
\end{bmatrix} $$ is the ray transfer matrix for free space propagation between two points of distance d.

2. $$\begin{bmatrix}
    1 & 0\\
    -1/f & 1
\end{bmatrix} =\begin{bmatrix} 1 & 0\\ -2/R & 1 \end{bmatrix} $$
is the ray matrix for a thin lens of focal length $f$ (or radius of curvature $R$, $f = R/2$).

In order to calculate the total transfer matrix, we multiply the matrices as follows: The first element is the last object that the beam reaches ($M_{n}$), and the last one is the starting point of the beam propagation ($M_{1}$). Thus, $M_{total} = M_{n} \cdot M_{n-1} \cdot M_{n-2} \cdot ... \cdot M_{2} \cdot M_{1}$.


## The Program

The program features a class function containing the most important parameters of the cavity, a main function and six custom functions.

The class is written as follows:

```python
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
```

The class' arguments are the wavelength (lambda0), the radii of curvature (R1, R2), the cavity length (d), the reflectivities of the mirrors (ref1, ref2) and also g1 and g2, which will be used when evaluating cavity stability. A cavity is stable when

$$
\begin{equation}
    0 \leq g_1 \cdot g_2 \leq 1 \ , \ g_i = 1 - \frac{d}{R_i}
\end{equation}
$$

The class is helpful when we want to call multiple parameters in our custom functions.

The first custom function handles matrix multiplication for a photon round-trip inside the cavity

```python
def matrices(cavity, y):  #Define ABCD Matrices for the Round-Trip
    M1 = sp.Matrix([[1, y], [0, 1]])
    M2 = sp.Matrix([[1, 0], [-(2/cavity.R1), 1]])
    M3 = sp.Matrix([[1, cavity.d], [0, 1]])
    M4 = sp.Matrix([[1, 0], [-(2/cavity.R2), 1]])
    M5 = sp.Matrix([[1, (cavity.d)-(y)], [0, 1]])

    M_rt = sp.Matrix(M1 * M2 * M3 * M4 * M5) #Round trip matrix
    return M_rt
```

which is shown also in this schematic

![cavity](https://cdn.shopify.com/s/files/1/1026/4509/files/Screenshot_2023-12-20_at_6.12.01_PM.png?v=1703124732)

The distance between the first mirror (left) and the vertical dashed line is of distance $y$. Because the starting point of the round trip is the dashed line, $M_1$ (as seen in the function) is the first element of multiplication, respecting optical elements multiplication. The output of this function is the total round-trip matrix $M_{rt}$.

Having found the total matrix for the round trip, we can calculate the value of $y$ in order to find the position $z_0$ of the minimum beam waist, $w_0$. At $z_0$, the wavefront is flat. As a result, $R(z_0) = \infty$, and using eq. (13), $D=A$. This is implemented in this function

```python
def minimum_waist_position(cavity):  #Find z0 where w(z) is w0 (minimum waist)
    M_rt = matrices(cavity, y)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    z0 = solve(D - A, y)
    z0 = float(z0[0])
    return z0
```

Having found $z_0$ from the above function, we can find $w_0$ with two methods. The first one is shown in this function

```python
def minimum_waist_1(cavity, lambda0):      #First Method to find w0, this uses the solution to
    z0 = minimum_waist_position(cavity)    #1/q=(C + D(1/q))/(A + B(1/q))
    M_rt = matrices(cavity, z0)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    w0 = np.sqrt(lambda0/ np.pi) * \
         np.sqrt(abs(float(B)) / (np.sqrt(1 - ((float(A) + float(D))/2)**2)))
    return w0
```

This function uses the theory I have described above and utilizes eq.(12).

Moreover, we can test if the result is correct with another method, shown in the next function. The $q$ factor, introduced at eq.(3), can be used in the following pattern

```python
def minimum_waist_2(cavity):      #q=z+izR, z=0 at the start
    z0 = minimum_waist_position(cavity)
    M_rt = matrices(cavity, z0)
    A, B, C, D = M_rt[0, 0], M_rt[0, 1], M_rt[1, 0], M_rt[1, 1]
    q1_2 = solve(y - (float(A) * y + float(B))/(float(C) * y + float(D)), y)
    z_R = q1_2[0] * 1j
    z_R = float(z_R)
    w_0 = np.sqrt(cavity.lambda0 * z_R / np.pi)
    return w_0, z_R
```

With this function, we can also find the Rayleigh range, which is very useful as at this distance, the beam waist is $w(z) = \sqrt{2}w_0$ and we can also extract information for the total angular spread $\Theta$, but this lies outside the scope of the project.

The second to last function handles calculations for:

- $\text{FSR} = \frac{c}{2d}$, free spectral range

- $\mathcal{F} = \frac{π(R_1 R_2)^{1/4}}{1 - \sqrt{R_1 R_2}} $, the finesse of the cavity, with $R_i$ being the reflectivities

- $Δv = \frac{FSR}{\mathcal{F}}$, the Full Width at Half Maximum

- $v_0 = \frac{c}{λ}$, the frequency of the laser

- $Q = \frac{v_0}{Δv}$, the Quality factor, and

- $τ_γ = \frac{1}{2π Δv}$, the photon lifetime inside the cavity.

These calculations are handled in the next functions and they are saved in a .txt file in the main() function.


```python
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
```

These functions provide the values required by the final function, as this one handles the graphs of $w(z), R(z)$, the intensities of the beam inside the cavity at $z=0$ (at the first mirror) and at $z=z_0$ (where $w(z) = w_0$) (eq.(2)), the physical electric field, which is given by the real part of the phasor field amplitude (eq.(1)) times a time factor $e^{iωt}$ and the intensity given by

$$
\begin{equation}
    I = I_0 (\frac{w_0}{w(z)})^2 \text{exp}(-\frac{r^2}{w(z)^2})
\end{equation}
$$

```python
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
```

The output of this function is a .png file, containing all six plots.
