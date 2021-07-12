import numpy as np
from scipy import interpolate
from scipy.constants import pi, c
import math
import cmath
from scipy.io import loadmat, savemat
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lambda0 = np.linspace(1545,1555,2000)*1e-9
height_wg=0.4e-6
w_wg=0.69e-6
radius = 20*1e-6
alpha_dB_cm = 0.40
loaddata=1 #boolean specifying if a data file can be loaded

## CONVERSION ##############################################################
alpha = alpha_dB_cm*10*math.log(10)                  # per m

if loaddata:
    from tkinter import Tk     # from tkinter import Tk for Python 3.x
    from tkinter.filedialog import askopenfilename

    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

    ## EFFECTIVE INDEX CALCULATIONS ##################################
    #Eventually plan to replace this section with some Lumerical scripting
    data=loadmat(filename)
    f=data['f']
    neff=data['neff']
    ng=data['ng']
    f=[f[i] for i in range(len(f))][0]
    neff=[neff[i] for i in range(len(neff))][0]
    ng=[ng[i] for i in range(len(ng))][0]
    #f=np.sort([data[i][0] for i in range(len(data))])
    #neff=np.sort([data[i][1] for i in range(len(data))])
else:
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
    mode.setanalysis('wavelength',lambda0[0])
    mode.findmodes()
    mode.selectmode(1)
    mode.setanalysis("track selected mode",1)
    mode.setanalysis('stop wavelength',lambda0[-1])
    mode.setanalysis('number of points',10)
    mode.setanalysis('number of test modes',4)
    mode.frequencysweep()
    neff=mode.getdata('frequencysweep','neff')
    f=mode.getdata('frequencysweep','f')
    vg=mode.getdata('frequencysweep','vg')
    f=[f[i][0] for i in range(len(f))]
    neff=[neff[i][0] for i in range(len(neff))]
    ng=[c/vg[i][0] for i in range(len(vg))]
    mdic={"f":f,"neff":neff,"ng":ng}
    savemat('neff_f_{}_{}.mat'.format(height_wg*1e9,w_wg*1e9),mdic)


## Increased sampling of neff
interp_neff=interpolate.interp1d(f,neff,kind='cubic')



## Round trip distance, propagation constant
L_rt=2*pi*radius*3
beta=interp_neff(c/lambda0) * 2*pi/lambda0
gamma=0 #lossless coupling
FSR=c/(ng[0]*L_rt)
print('FSR={} GHz'.format(np.round(FSR*1e-9)))

##Loop space for varying kappa
kappa1=0.6
kappa2=0.1
#for kappa2 in np.linspace(0.1,0.5,4):
kappaList=np.array([kappa1,kappa2,kappa2,kappa1])
## Transfer matrix elements
C=[]
S=[]
for i in range(0,4):
    C.append(cmath.sqrt(((1-kappaList[i])*(1-gamma))))
    S.append(1/1j*cmath.sqrt(((1-gamma*1j)*kappaList[i])))
xi=np.multiply(np.exp(alpha*L_rt/2),np.exp(-1j*beta*L_rt))
## Compute transmission
E_drop=S[0]*S[1]*S[2]*S[3]*np.power(xi,1.5)/(1-C[0]*C[1]*xi-C[1]*C[2]*xi-C[2]*C[3]*xi+C[0]*C[2]*np.power(xi,2)+C[1]*C[3]*np.power(xi,2)-C[0]*C[3]*np.power(xi,3)+\
    C[0]*C[1]*C[2]*C[3]*np.power(xi,2))
T_drop=np.power(abs(E_drop),2)

E_thru=(C[0]-C[1]*xi-C[0]*C[1]*C[2]*xi+C[2]*np.power(xi,2)-C[0]*C[2]*C[3]*xi+C[1]*C[2]*C[3]*np.power(xi,2)+C[0]*C[1]*C[3]*np.power(xi,2)-C[3]*np.power(xi,3))\
    /(1-C[0]*C[1]*xi-C[1]*C[2]*xi-C[2]*C[3]*xi+C[0]*C[2]*np.power(xi,2)+C[1]*C[3]*np.power(xi,2)-C[0]*C[3]*np.power(xi,3)+C[0]*C[1]*C[2]*C[3]*np.power(xi,2))
T_thru=np.power(abs(E_thru),2)
## Data in dB
T_drop_dB=10*np.log10(T_drop)
T_thru_dB=10*np.log10(T_thru)

#Plot
plt.rcParams.update({'font.size': 22})
plt.plot(lambda0*1e9,T_thru_dB)

plt.plot(lambda0*1e9,T_drop_dB)
plt.legend(['T_thru','T_drop'])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (dB)')
#plt.legend(['\u03BA_2=0.1','\u03BA_2=0.2','\u03BA_2=0.3','\u03BA_2=0.4'],loc="upper right")
plt.show()

"""
### Space to create a heat map of important parameters: ER, IL, shape factor vs kappa1,kappa2
numsplit=18
ER=np.zeros((numsplit,numsplit))
IL=np.zeros((numsplit,numsplit))
SF=np.zeros((numsplit,numsplit))
kappa1List=np.linspace(0.1,.9,numsplit)
kappa2List=np.linspace(0.1,.9,numsplit)

for j in range(len(kappa1List)):
    for k in range(len(kappa2List)):
        kappa1=kappa1List[j]
        kappa2=kappa2List[k]
        kappaList=[kappa1,kappa2,kappa2,kappa1]
        ## Transfer matrix elements
        C=[]
        S=[]
        for i in range(0,4):
            C.append(cmath.sqrt(((1-kappaList[i])*(1-gamma))))
            S.append(1/1j*cmath.sqrt(((1-gamma*1j)*kappaList[i])))
        xi=np.multiply(np.exp(alpha*L_rt/2),np.exp(-1j*beta*L_rt))
        ## Compute transmission
        E_drop=S[0]*S[1]*S[2]*S[3]*np.power(xi,1.5)/(1-C[0]*C[1]*xi-C[1]*C[2]*xi-C[2]*C[3]*xi+C[0]*C[2]*np.power(xi,2)+C[1]*C[3]*np.power(xi,2)-C[0]*C[3]*np.power(xi,3)+\
            C[0]*C[1]*C[2]*C[3]*np.power(xi,2))
        T_drop=np.power(abs(E_drop),2)

        E_thru=(C[0]-C[1]*xi-C[0]*C[1]*C[2]*xi+C[2]*np.power(xi,2)-C[0]*C[2]*C[3]*xi+C[1]*C[2]*C[3]*np.power(xi,2)+C[0]*C[1]*C[3]*np.power(xi,2)-C[3]*np.power(xi,3))\
            /(1-C[0]*C[1]*xi-C[1]*C[2]*xi-C[2]*C[3]*xi+C[0]*C[2]*np.power(xi,2)+C[1]*C[3]*np.power(xi,2)-C[0]*C[3]*np.power(xi,3)+C[0]*C[1]*C[2]*C[3]*np.power(xi,2))
        T_thru=np.power(abs(E_thru),2)
        ## Data in dB
        T_drop_dB=10*np.log10(T_drop)
        T_thru_dB=10*np.log10(T_thru)
        ER[j][k]=np.max(T_drop_dB)-np.min(T_drop_dB)
        IL[j][k]=np.max(T_drop_dB)
        try:
            peakidx=np.argwhere(T_drop_dB==np.max(T_drop_dB))
            temp=np.where(T_drop_dB[peakidx]-1>T_drop_dB)
            idx_1db=temp[1][np.argmax(temp[1]>peakidx)]
            temp=np.where(T_drop_dB[peakidx]-10>T_drop_dB)
            idx_10db=temp[1][np.argmax(temp[1]>peakidx)]
            SF[j][k]=(np.abs(lambda0[idx_1db]-lambda0[peakidx])/(np.abs(lambda0[idx_10db]-lambda0[peakidx])))
        except:
            SF[j][k]=0



fig,(ax1,ax2,ax3)=plt.subplots(1,3)

x_axis_labels=np.around(kappa2List,2)
y_axis_labels=np.around(kappa1List,2)
sns.heatmap(ER,cmap='rocket',ax=ax1,cbar=True,cbar_kws={"orientation": "horizontal"},xticklabels=x_axis_labels,yticklabels=y_axis_labels)
sns.heatmap(IL,cmap='rocket',ax=ax2,cbar=True,cbar_kws={"orientation": "horizontal"},xticklabels=x_axis_labels,yticklabels=y_axis_labels)
sns.heatmap(SF,cmap='rocket',ax=ax3,cbar=True,cbar_kws={"orientation": "horizontal"},xticklabels=x_axis_labels,yticklabels=y_axis_labels)
ax1.invert_yaxis()
ax2.invert_yaxis()
ax3.invert_yaxis()
ax1.set_ylabel('Kappa1')
ax1.set_xlabel('Kappa2')
ax2.set_ylabel('Kappa1')
ax2.set_xlabel('Kappa2')
ax3.set_ylabel('Kappa1')
ax3.set_xlabel('Kappa2')
ax1.collections[0].colorbar.set_label("Extinction ratio (dB)")
ax2.collections[0].colorbar.set_label("Insertion Loss (dB)")
ax3.collections[0].colorbar.set_label("Shape Factor")
plt.show()"""