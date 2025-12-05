import pytest
import sympy as sp
from scipy import constants
from project import matrices, minimum_waist_position, calculations, minimum_waist_1, minimum_waist_2

class Dummy():
        lambda0 = 1064e-6 #Wavelength in mm
        R1 = 100          #ROC of Mirror 1 in mm
        R2 = 100          #ROC of Mirror 2 in mm
        ref1 = 0.99       #Reflectivity of Mirror 1
        ref2 = 0.999      #Reflectivity of Mirror 2
        d = 150           #Cavity Length in mm

def test_matrices():
    cavity = Dummy()
    y = 10
    M_rt = matrices(cavity, y)
    M_expected = sp.Matrix([[1, y],
                [0, 1]])* sp.Matrix([[1, 0],
                [-(2/cavity.R1), 1]])* sp.Matrix([[1, cavity.d],
                [0, 1]])* sp.Matrix([[1, 0],
                [-(2/cavity.R2), 1]])*sp.Matrix([[1, (cavity.d)-(y)], [0, 1]])
    assert M_expected == M_rt

def test_minimum_waist_position():
    cavity = Dummy()
    z0 = minimum_waist_position(cavity)
    assert z0 == float(cavity.d/2) #R1=R2, so z0 is in the middle of the cavity

def test_calculations():
    cavity = Dummy()
    values = calculations(cavity)
    assert "FWHM" in values
    assert "Photon lifetime" in values
    assert "Q" in values
    assert "Finesse" in values
    assert "FSR" in values
    c = constants.speed_of_light
    fsr_expected = c / (2.0 * (cavity.d * 1E-3))
    assert f"{fsr_expected:.3e}" in values

def test_minimum_waist_1():
    cavity = Dummy()
    w0 = minimum_waist_1(cavity, cavity.lambda0)
    w0_test = 0.1211
    assert str(w0_test) == f"{w0:.04f}"

def test_minimum_waist_2():
    cavity = Dummy()
    w0, z_R = minimum_waist_2(cavity, cavity.lambda0)
    w0_test = 0.1211
    z_R_test = 43.3013
    assert str(w0_test) == f"{w0:.04f}"
    assert str(z_R_test) == f"{z_R:.04f}"
