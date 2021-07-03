import numpy as np
def FindIndex(x,lambda0):
    #x: fraction of Al content
    #lambda0: wavelength (m)
    h_plank = 4.135667662e-15
    c_speed = 2.998e8   #m/s

    # The following are the coefficients for refractive index computation
    def A0(y):
        return 6.3+19.0*y
    def B0(y):
        return 9.4-10.2*y
    def E0(y):
        return 1.425 + 1.155*y + 0.37*np.square(y)
    def E0_plus_D0(y):
        return 1.765 + 1.115*y + 0.37*np.square(y)

    chi_s0    = h_plank*c_speed/(lambda0*(E0_plus_D0(x)))
    chi       = h_plank*c_speed/(lambda0*E0(x))
    f_chi     = (2 - np.sqrt(1 + chi) - np.sqrt(1 - chi))/np.square(chi)
    f_chi_s0  = (2 - np.sqrt(1 + chi_s0) - np.sqrt(1 - chi_s0))/np.square(chi_s0)
    # n is the refractive index
    n =  np.sqrt(A0(x)*(f_chi + f_chi_s0/2*np.power(E0(x)/E0_plus_D0(x),1.5)) + B0(x))
    #print('index of refraction for this material is ' + str(n))
    return n
