import numpy as np
import sys
import cmath as m
from matplotlib import pyplot as plt

def square(x):
    return x*x

''' 
get the complex speed of sound given atten. coefficient alpha
'''
def complexC(cReal, Alpha):
    delta = Alpha / 54.58
    cImag = delta * cReal
    return complex(cReal, cImag)


''' impedance of a material layer with attenuation coefficient alpha and density rho and sound speed c and angle theta
'''

def Z(alpha, rho, c, theta):
    ret = rho * c / m.sin(theta)
    return ret


''' snell's law:
cos(theta2) = c2/ c2 cos(theta1)
'''
def Snell(c1, c2, theta1):
    arg = c2 / c1 * m.cos(theta1)
    # branch cut
    if (arg.real > 1) and (arg.imag == 0):
        print("on the branch cut")
        print(theta1)
#        raise ValueError("arg is on branch cut")
        return m.acos(arg)
    return m.acos(arg)


def Ztot(Zp, Zs, thetaS):
    return Zp * square(m.cos(2*thetaS)) + Zs * square(m.sin(2*thetaS))


class Medium:
    # use -1 default thickness to allow for checking
    def __init__(self, rho, c, alpha=0, thickness=-1):
        self.rho = rho
        self.c = complexC(c, alpha)
        self.alpha=alpha
        self.thickness = thickness
        self.delta = alpha / 54.58

    # get impedance for an incident ray with refracted angle theta
    def getZ(self, theta):
        Z  = self.rho * self.c / m.sin(theta)
        return Z


class SolidMedium(Medium):
    def __init__(self, rho, c, alpha, thickness, c_shear, alpha_shear):
        Medium.__init__(self, rho, c, alpha, thickness)
        self.c_shear = c_shear
        self.alpha_shear = alpha_shear

    def getZ(self, theta_shear, theta_p):
        Zs = self.rho * self.c_shear / m.sin(theta_shear)
        Zp = self.rho * self.c / m.sin(theta_p)
        Z = Zp * square(m.cos(2*theta_shear)) + Zs * square(m.sin(2*theta_shear))
        
        

class Ray:
    def __init__(self, f, medium, A, theta):
        self.f = f
        self.medium = medium
        self.A = A
        self.theta = theta
        self.snell_constant = m.cos(self.theta) / medium.c

    def getRefractedAngle(self, new_medium):
        snell = self.snell_constant
        if type(new_medium) is SolidMedium:
            c = new_medium.c_shear
        else:
            c = new_medium.c
        angle = m.acos(c * snell)
        return angle

def RfromZ(Z1, Z2):
    return (Z2 - Z1)/(Z2+Z1)


class LayeredBottomInteraction:
    def __init__(self, incidentRay, medium2, medium3):
        self.incidentRay = incidentRay
        self.medium2 = medium2
        self.medium3 = medium3

    def getRCoeff(self):
        primaryRay = self.incidentRay
        med1 = primaryRay.medium
        med2 = self.medium2
        med3 = self.medium3
        theta1 = primaryRay.theta
        theta2 = primaryRay.getRefractedAngle(med2)
        R12 = RfromZ(med1.getZ(theta1), med2.getZ(theta2))
        secondaryA = primaryRay.A * R12
        secondaryRay = Ray(primaryRay.f, med2, secondaryA, theta2)
        
        theta3 = secondaryRay.getRefractedAngle(med3)
        R23 = RfromZ(med2.getZ(theta2), med3.getZ(theta3))
        tertiaryA = R23 * secondaryRay.A
        tertiaryRay = Ray(primaryRay.f, med3, tertiaryA, theta3)

        Z1 = med1.getZ(theta1)
        Z2 = med2.getZ(theta2)
        Z3 = med3.getZ(theta3)

        print("Z1, Z2, Z3: ", ",".join([str(x) for x in [Z1, Z2, Z3]]))

        thickness = med2.thickness
        c = med2.c
        wavelength = c / primaryRay.f
        phi = 2*m.pi*(thickness / wavelength) * m.sin(theta2)
        RCoeff_num = complex(Z2*(Z3 - Z1), -(square(Z2) - Z1*Z3)*m.tan(phi))
        RCoeff_den = complex(Z2*(Z3 + Z1), -(square(Z2) + Z1*Z3)*m.tan(phi))
        RCoeff = RCoeff_num / RCoeff_den
        return RCoeff

''' function to compute the reflectivity coefficient for a given solid layer in terms of
the solid layer shear wave speed, compression wvae speed, p-attenuation, s-attenuation, and medium density
and the grazing angle theta
we assume liquid layer has speed of sound 1500m/s and rho = 1000 kg /m^3
'''
def R(Cp, Cs, AlphaP, AlphaS, rho, theta):
    Cp = complexC(Cp, AlphaP)
    Cs = complexC(Cs, AlphaS)
    if (theta == 0):
        return 1
    try:
        thetaP = Snell(1500, Cp, theta)
    except ValueError:
        print("hi")
    try:
        thetaS = Snell(1500, Cs, theta)
    except ValueError:
        print("hi")
    Zp = Z(AlphaP, rho, Cp, thetaP)
    Zs = Z(AlphaS, rho, Cs, thetaS)
    Zwater = Z(0, 1000, 1500, theta)
    Ztotal = Ztot(Zp, Zs, thetaS)
    numerator = Ztotal - Zwater
    denom = Ztotal + Zwater
    ret = (Ztotal - Zwater) / (Ztotal + Zwater)
   # return (Ztotal - Zwater) / (Ztotal + Zwater)
    return ret
   

def BottomLoss(RCoeff):
    ret =  (-10 * m.log10(square(abs(RCoeff))))
#    if ret.real < 0:
 #       print("negative ret")
  #      print("Rcoeff: ", abs(RCoeff))
    return ret


if __name__ == '__main__':
    plt.figure()
    plt.subplot(4, 2, 1)
    showFiga(100)
    plt.subplot(4, 2, 2)
    showFigb(100)
    plt.subplot(4, 2, 3)
    showFigc(100)
    plt.subplot(4, 2, 4)
    showFigd(100)
