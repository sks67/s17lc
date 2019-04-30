#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Light curve from paper. This is to make sure all the equations
# as they appear in the paper are consistent with my program.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import matplotlib.pyplot as plt

#<><><><><><><><><> F U N C T I O N S <><><><><><><><><><><><>#

def radius_velocity(t, n0, mej, e51, sntype='cc'):
    """
    Returns the shock radius and velocity for a given density, 
    energy and ejecta mass.

    Parameters:
    ----------

    t : age of SNR (in years)
    n0 : density (atoms/cc)
    e51 : energy (in units of 10^51 ergs)
    mej : ejecta mass (in units of Msun)
    sntype : ['ia' or 'cc']
    """

#Characteristic variables - TM99                                                                                                                      
    tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years                                                                                      
    rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs                                                                                                     
    vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s                                                                                                         
    
    #Check whether to apply Ia or CC physics of TM99                                                                                                         

    if (sntype=='ia'):
        tstar_st = 0.481
        t_st0 = tstar_st*tch
        tstar = t/tch
        vstar = 0.805*tstar**(-(3.0/10.0)) if t<t_st0 else 0.569*((1.42*tstar - 0.2935)**(-3.0/5.0))
        rstar = 1.15*tstar**(7.0/10.0) if t<t_st0 else (1.42*tstar - 0.2935)**(2.0/5.0)
        vs = vstar*vch  #in km/s                                                                                                                 
        rs = rstar*rch  #in pc                                                                                                                      
    
    elif(sntype=='cc'):
        tstar_st = 0.424
        t_st0 = tstar_st*tch
        tstar = t/tch
        vstar = 0.906*tstar**(-(1.0/4.0)) if t<t_st0 else 0.569*((1.42*tstar - 0.28)**(-3.0/5.0))
        rstar = 1.20906*tstar**(3.0/4.0) if t<t_st0 else (1.42*tstar - 0.28)**(2.0/5.0)
        vs = vstar*vch  #in km/s                                                                                                                
        rs = rstar*rch  #in pc

    return (rs, vs)

    
def luminosity(rs, vs, n0, epse, p, nu):
    """
    Luminosity of SNR, given radius, velocity and other parameters

    (NOTE: the electron acceleration and magnetic field terms have a separate dependence
           on ISM density, n0. The radius and velocity also depend on n0. So be careful about
           pluggin in an arbitrary value of 'rs' and 'vs' and getting luminosity, because the
           'rs' and 'vs' has to be consistent with 'n0'.)

    Parameters:
    ----------

    rs: shock radius (in pc)
    vs: shock velocity (in km/s)
    n0: ISM density (in atoms/cc)
    epse: electron acceleration efficiency 
    p: electron spectral index (stay within p=2.2-2.5)
    nu: observing frequency (in GHz)

    """
    #~~~~~~ Constants ~~~~~~~~~~~~~~~#
    ef = 0.38  #Emission filling factor
    c1 = 6.27e18 #Pacholyzyk (1980) constants. Note that these depend on p.
    c5 = 9.68e-24 #Note for SKS: Correction in paper
    c6 = 8.1e-41
    mu = 1.4 
    mp = 1.67e-24 #g
    rho0 = mu*mp*n0 #density in g/cm^3
    me =  9.11e-28  #in grams                                                                                                                   
    c = 3.0e10 #cm/s
    eta = 4. #compression factor
    Em = me*c**2 #Rest energy of electrons
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~ CGS units for r,v  ~~~~~~~~~~~~~~~#
    rs_cgs = rs*3.086e18 #cm
    vs_cgs = vs*1.0e5 #cm/s

    #~~~~~~ Electron Acceleration Stuff (Eq A2) ~~~~~~~~~~~~~~~#
    N_0 = (p-2.)*epse*rho0*(vs_cgs**2)*(Em**(p-2))

    #~~~~~~ Magnetic Field Stuff~~~~~~~~~~~~~~~#
    B0 = 9.0e-6*n0**0.47 #Eq A4 #Note for SKS: Correct this in the paper
    vA = B0/np.sqrt(4.*3.14*rho0) #Alfven Velocity term (Sec A2, 2nd para)
    MA = vs_cgs/vA #Alfven Mach Number
    
    epsCR = 0.01*(0.15*MA + 6.)  #Cosmic ray acceleration efficiency
    epsCR[epsCR>0.1] = 0.1  #See Section A2
    epsub = (epsCR/2.)*((vs_cgs/c) + (1./MA)) #Eq A7
    Bu = np.sqrt(8.*3.14*epsub*rho0*vs_cgs**2) #Eq A3
    B = Bu*np.sqrt((1. + 2.*(eta**2))/3.) #Eq A8
    B[np.where(B<4.*B0)] = 4.*B0 #We assert that SNR magnetic field should atleast be a simple compression of the ISM field (B0)

    #~~~~~~ Luminosity ~~~~~~~~~~~~~~~~~~~#
    s = ef*(4./3.)*(rs_cgs) #Eq A9
    nu1 = 2.*c1*((s*c6*N_0)**(2./(p+4)))*(B**((p+2)/(p+4))) #Eq A11
    
    L1 = 4.0*(3.14**2)*(rs_cgs**2)*(B**(-0.5))*(c5/c6) #Breaking Eq A10 into 3 parts - L1, L2, L3 - for simplicity
    L2 = 1 - np.exp(-(nu/nu1)**(-(p+4)/2.))

    #This is a short step I am adding because the above equation for L2 runs into numerical errors
    #for nu>>nu1 (i.e. much later in the SNR age). This is particularly apparent for very low densities (e.g. n0<0.01)
    #To avoid this for now, I enforce the approximation: e^-x = 1 - x (for x<<1). May revisit this at some point
    L2[nu>20.*nu1] = ((nu/nu1)**(-1.0*(p+4)/2.))[nu>20.*nu1]
    L3 = (nu/(2.*c1))**(5./2.)

    return L1*L2*L3 #Eq A10
    
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#


