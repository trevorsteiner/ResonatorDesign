import numpy as np
from scipy import interpolate
from scipy.constants import pi, c
import math
import cmath
from scipy.io import loadmat
import pandas as pd
import matplotlib.pyplot as plt

lambda0 = np.linspace(1200,1800,200)*1e-9
height_wg=0.4e-6
w_wg=0.69e-6
n_clad=1.444
radius = 600*1e-6
alpha_dB_cm = 0.40
kappaList=np.array([0.5,0.14,0.14,0.5])

## CONVERSION ##############################################################
alpha = alpha_dB_cm*10*math.log(10)                  # per m


## LUMERICAL Automation for effective index
import imp
lumapi = imp.load_source("lumapi.py", "C:\\Program Files\\Lumerical\\v211\\api\python\\lumapi.py") # Windows
mode=lumapi.MODE(hide=False)

## Add AlGaAs material with Al concentration x to Lumerical
from FindIndex import FindIndex
x=0.2
wgMat='AlGaAs'
material=mode.addmaterial('Sampled 3D Data')
mode.setmaterial(material,"Name",wgMat)
indexlist=FindIndex(x,lambda0)
indexdata=[[c/lambda0[i],np.power(indexlist[i],2)] for i in range(len(lambda0))]
mode.setmaterial(wgMat,"sampled data",np.asarray(indexdata))

mode.switchtolayout()
mode.deleteall()
mode.addrect()
mode.set('name','Waveguide')
mode.set('x span',w_wg)
mode.set('y span',height_wg)
mode.set('x',0)
mode.set('y',0)
mode.set('material',wgMat)
    
mode.addrect()
mode.set('name','Oxide')        
mode.set('x span',w_wg+40e-6)
mode.set('x',0)
mode.set('y',0)
mode.set('y span',height_wg+6e-6)
mode.set('material','SiO2 (Glass) - Palik')
mode.set('override mesh order from material database',1)
mode.set('mesh order',3)
mode.set('alpha',0.2)

mode.addmesh()  # mesh override, higher resolution in the waveguide.
mode.set('x span',w_wg+0.1e-6)
mode.set('y span',height_wg+0.1e-6)
mode.set('x',0)
mode.set('y',0)
mode.set('dx',10e-9)
mode.set('dy',10e-9)
    
mode.addfde() #  create simulate mesh
mode.set('x span',w_wg+10e-6)
mode.set('y span',height_wg+6e-6)
mode.set('x',0)
mode.set('y',0)
mode.set('mesh cells x',100) 
mode.set('mesh cells y',200)
mode.set('x min bc','PML')
mode.set('y min bc','PML')
mode.set('x max bc','PML')
mode.set('y max bc','PML')

mode.run()



# ## EFFECTIVE INDEX CALCULATIONS ##################################
# #Eventually plan to replace this section with some Lumerical scripting
# xtest=loadmat('480nmbend_neff.mat')
# data=np.asarray(xtest['ans'])

# f=np.sort([data[i][0] for i in range(len(data))])
# neff=np.sort([data[i][1] for i in range(len(data))])

## Increased sampling of neff
ang_neff=interpolate.interp1d(f,neff,kind='cubic')
fnew=np.linspace(f[0],f[-1],len(lambda0))

## Round trip distance, propagation constant
L_rt=2*pi*radius
beta=ang_neff(fnew) * 2*pi/lambda0
gamma=0 #lossless coupling

## Transfer matrix elements
C=[]
S=[]
for i in range(0,4):
     C.append(cmath.sqrt(((1-kappaList[i])*(1-gamma))))
     S.append(1/1j*cmath.sqrt(((1-gamma*1j)*kappaList[i])))
xi=np.multiply(np.exp(alpha*L_rt/2),np.exp(-1j*beta*L_rt))
## Compute transmission
T_drop=S[0]*S[1]*S[2]*S[3]*xi/(1-C[0]*C[1]*np.power(xi,2)+C[0]*C[2]*np.power(xi,2)+\
    C[1]*C[3]*np.power(xi,3)+C[0]*C[1]*C[2]*C[3]*np.power(xi,2))

## Data in dB
T_drop_dB=10*np.log10(T_drop)

#Plot
plt.plot(lambda0,T_drop_dB)
plt.show()