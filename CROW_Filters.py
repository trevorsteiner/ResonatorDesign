import numpy as np
from scipy import interpolate
from scipy.constants import pi, c
import math
from scipy.io import loadmat
import pandas as pd
import matplotlib.pyplot as plt

lambda0 = np.linspace(1200,1800,200)*1e-9
height_wg=0.4e-6
w_rr=0.69e-6
n_wg= 3.2658
n_clad=1.444
radius = 60*1e-6
alpha_dB_cm = 0.40
kappaList=np.array([0.5,0.14,0.14,0.5])

## CONVERSION ##############################################################
alpha = alpha_dB_cm*10*math.log(10)                  # per m


## EFFECTIVE AND GROUP INDEX CALCULATIONS ##################################
#Eventually plan to replace this section with some Lumerical scripting
xtest=loadmat('480nmbend_neff.mat')
data=np.asarray(xtest['ans'])

f=np.sort([data[i][0] for i in range(len(data))])
neff=np.sort([data[i][1] for i in range(len(data))])

ang_neff=interpolate.interp1d(f,neff,kind='cubic')
fnew=np.linspace(f[0],f[-1],len(lambda0))
L_rt=2*pi*radius
beta=ang_neff * 2*pi/lambda0
gamma=0 #lossless coupling
C=[]
S=[]
for i in range(0,4):
     C.append(math.sqrt(((1-kappaList[i])*(1-gamma))))
     S.append(1/1j*math.sqrt(((1-gamma*1j)*kappaList[i])))
xi=math.exp(alpha*L_rt/2)*exp(-1j*beta*L_rt)
T_drop=S(1)*S(2)*S(3)*S(4)*xi./(1-C(1)*C(2)*xi.^2+C(1)*C(3)*xi.^2+C(2)*C(4)*xi.^3+C(1)*C(2)*C(3)*C(4)*xi.^2)
T_drop_dB=10*math.log10(T_drop)

plt.plot(lambda0,T_drop_dB)
plt.show()