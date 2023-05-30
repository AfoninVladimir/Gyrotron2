import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import sys
import time
from scipy.special import jn as Bessel_JN
from scipy.special import factorial

from band.zbndset import zbndset
from band.zbndcpy import zbndcpy
from band.zbndtrn_v1_1 import zbndtrn
from band.zbndcrm_v1_1 import zbndcrm
from band.zbndzdrm_v1_1 import zbndzdrm
from band.zbndapcol_v1_1 import zbndapcol
from scipy import special
import numpy as np

from LinearGyrotronFDMF_0_11 import LinearGyrotronFDMFv0_11

Ri = []
c0 = 2.99792458 ** 8  # скорость света
mu0 = 4 * np.pi * 10**7
epsilon0 = 1.602176620898 * 10 ** (-12)
e0 = 1.602176620898 * 10 ** (-19)  # заряд электрона
m0 = 9.10938356 * 10 ** (-31)  # масса электрона
mm = 10 ** (-3)
ms = 10 ** (-3)
GHz = 10 ** 9
kV = 10 ** 3
nS = 10 ** (-9)
MHz = 10 ** 6
KHz = 10 ** 3
VBeam = 100000.0  # напряжение пучка
IBeam = 5.2  # ток пучка
RBeam = 8.271 * mm  # радиус пучка
g = 1.2  # питч-фактор
"""релятивистские факторы"""
gamma0 = ((1 + e0) * VBeam) / (m0 * c0 ** 2)  # гамма-фактор
beta0 = np.sqrt(1 - (1 / gamma0 ** 2))  # beta-фактор
betalong = beta0 / np.sqrt(1 + (g ** 2))  # продольный бета -фактор
betatrans = g * betalong  # поперечный бета-фактор

m = 28  # азимутальный индекс
n = 12  # радиальный индекс

R0 = 20.77 * mm  # радиус однородной части резонатора
B0 = 7.2
# numn = special.jn_zeros(n, m)[0]
numn = 73.9521
# omega0 = numn * c0 / R0
omega0 = 0.0+0.0j
f0 = omega0 / (2 / np.pi)
gyroratio = e0 / (m0 * gamma0)  # гироотношение
omegaH = gyroratio * B0
fH = omegaH/(2 * np.pi)
s = 1  # номер циклотронного резонанса
L = 12. * mm  # длина однородной части резонатора
zleft = -35.0 * mm  # координата левой границы расчетной области
zright = 85.0 * mm  # координата правой границы расчетной области
zstr = -25.0 * mm  # координата начала пространства взаимодействия
zfin = 5.0 * mm  # координата конца пространства взаимодействия
dz = 0.2 * mm  # шаг по координате
NModes = 4  # число мод в расчетах квазисобственных колебаний в резонаторе
sigma = 10.0 ** 12  # проводимость
OmegaL2divR2 = 0.0
BStart = 6.9
BFin = 7.7
DeltaB = 0.025
fB = 1.0

zi = []
num = zleft
while num < zright:
    # print(round(num, 10))
    zi.append(round(num, 10))
    num += dz
zi.append(zright)

zi = np.array(zi)
NN = len(zi)
lenB = NN
Bi = np.array([1.0] * lenB)

RZ = ["z <= (-31.5 * mm)", "11.0 * mm",
      "(-31.5 * mm) < z < (-24.0 * mm)", "(11.0 + (19.736 - 11) * (z/mm - (-31.5))/7.5) * mm",
      "(-24.0 * mm) <= z <= (-23.0 * mm)", "19.736 * mm",
      "(-23 * mm) < z < (-12 * mm)", "(-38.23014642003139 + np.sqrt(3481.0 - (12.0 + z/mm) ** 2)) * mm",
      "(-12.0 * mm) <= z <= (0 * mm)", "20.77 * mm",
      "(0.0 * mm) < z <= (17.7873 * mm)", "(295.7702855739705 - np.sqrt(75625.0 - (0.0 + z/mm)**2)) * mm",
      "(17.7873 * mm) < z < (80.0 * mm)", "(-938.4552895706714 + np.sqrt(925090.40258244 - (-80.0 + z/mm)**2)) * mm",
      "(80 * mm) <= z", "23.36 * mm"]


# Построение профиля резонатора
def Rz(zi, RZ):
    global Ri
    Rzt = np.array([0.0] * len(zi))
    n = -1
    for i in zi:
        n += 1
        z = i
        for j in range(0, len(RZ), 2):
            if eval(RZ[j]):
                Rzt[n] = eval(RZ[j + 1])
    Ri = Rzt

    # plt.figure(figsize = (12, 6.75))
    # plt.xlabel('z, мм')
    # plt.ylabel('R(z), мм')
    # plt.grid(True)
    # plt.plot(zi * 1000, Rzt * 1000, "b")
    # plt.plot(zi * 1000, -Rzt * 1000, "b")
    # plt.subplots_adjust(top = 0.98, bottom = 0.08, left = 0.08, right = 0.98, hspace = 0.2, wspace = 0.2)
    # plt.show()


Rz(zi, RZ)

# plt.plot(zi, Ri)
# plt.show()


istart = np.argsort(abs(zi-zstr))[0] + 1
ifin = np.argsort(abs(zi-zfin))[0] + 1
# print(istart, ifin)

Eigenfunctions = 1
PrintInfo = 1
PrintDebug = 1
PrintArraysDebug = 1
PrintCalcInfo = 1
ConditionNumberFound = 1
Criteria = 1

tol = 0

Shift0 = 0 + 0j

ConductanceThreshold = 1.0 * 10**11
iparam = [m, n, s, NN, NModes, Eigenfunctions, PrintInfo, PrintDebug, istart, ifin, PrintArraysDebug, PrintCalcInfo,
          ConditionNumberFound, Criteria]
rparam = [R0, B0, L, sigma, zleft, zright, zstr, zfin, dz, tol, omega0.real, omega0.imag, Shift0.real, Shift0.imag,
          VBeam, IBeam, RBeam, g, numn, ConductanceThreshold]

LinearGyrotronFDMFv0_11(iparam, rparam, zi, Ri, Bi)
