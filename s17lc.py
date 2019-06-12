#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Supplementary Python Program to calculate radio light curves
# of supernova remnants from ejecta-dominated through Sedov-Taylor
# stages.
#
# Paper: Sarbadhicary et al (2017), MNRAS, 464, 2326
# Erratum to Paper: https://doi.org/10.1093/mnras/stz1490
#
# Creators: S. Sarbadhicary, L. Chomiuk (MSU)
# Date: June 12, 2019
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import matplotlib.pyplot as plt

#<><><><><><><><><> F U N C T I O N S <><><><><><><><><><><><>#

def radius_velocity(t, n0, mej, e51, sntype='cc'):
    """
    Returns the shock radius and velocity for a given density,
    energy and ejecta mass.  The parameters can be any combination
    of matching-length arrays and scalars.

    Parameters:
    ----------

    t : [float, list, ndarray] age of SNR (in years)
    n0 : [float, list, ndarray] density (atoms/cc)
    e51 : [float, list, ndarray] energy (in units of 10^51 ergs)
    mej : [float, list, ndarray] ejecta mass (in units of Msun)
    sntype : ['ia' or 'cc']
    """

    t = np.asarray(t)
    n0 = np.asarray(n0)
    mej = np.asarray(mej)
    e51 = np.asarray(e51)

#Characteristic variables - TM99
    tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
    rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs
    vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s
    
    #Check whether to apply Ia or CC physics of TM99

    if (sntype=='ia'):
        tstar_st = 0.481
        t_st0 = tstar_st*tch
        tstar = t/tch
        if np.isscalar(tstar):
            vstar = 0.805*tstar**(-(3.0/10.0)) if t<t_st0 else 0.569*((1.42*tstar - 0.2935)**(-3.0/5.0))
            rstar = 1.15*tstar**(7.0/10.0) if t<t_st0 else (1.42*tstar - 0.2935)**(2.0/5.0)
        else:
            w = np.where(t>=t_st0)
            vstar = 0.805*tstar**(-(3.0/10.0))
            vstar[w] = 0.569*((1.42*tstar[w] - 0.2935)**(-3.0/5.0))
            rstar = 1.15*tstar**(7.0/10.0)
            rstar[w] = (1.42*tstar[w] - 0.2935)**(2.0/5.0)
        vs = vstar*vch  #in km/s
        rs = rstar*rch  #in pc
    
    elif(sntype=='cc'):
        tstar_st = 0.424
        t_st0 = tstar_st*tch
        tstar = t/tch
        if np.isscalar(tstar):
            vstar = 0.906*tstar**(-(1.0/4.0)) if t<t_st0 else 0.569*((1.42*tstar - 0.28)**(-3.0/5.0))
            rstar = 1.20906*tstar**(3.0/4.0) if t<t_st0 else (1.42*tstar - 0.28)**(2.0/5.0)
        else:
            w = np.where(t>=t_st0)
            vstar = 0.906*tstar**(-(1.0/4.0))
            vstar[w] = 0.569*((1.42*tstar[w] - 0.28)**(-3.0/5.0))
            rstar = 1.20906*tstar**(3.0/4.0)
            rstar[w] = (1.42*tstar[w] - 0.28)**(2.0/5.0)
        vs = vstar*vch  #in km/s
        rs = rstar*rch  #in pc

    return (rs, vs)


def luminosity(rs, vs, n0, epse, p, nu, b0ref=9.0):
    """
    Luminosity of SNR, given radius, velocity and other parameters.  The parameters
    can be any combination of matching-length arrays and scalars.

    (NOTE: the electron acceleration and magnetic field terms have a separate dependence
           on ISM density, n0. The radius and velocity also depend on n0. So be careful about
           pluggin in an arbitrary value of 'rs' and 'vs' and getting luminosity, because the
           'rs' and 'vs' has to be consistent with 'n0'.)

    Parameters:
    ----------

    rs: [float, list, ndarray] shock radius (in pc)
    vs: [float, list, ndarray] shock velocity (in km/s)
    n0: [float, list, ndarray] ISM density (in atoms/cc)
    epse: [float, list, ndarray] electron acceleration efficiency 
    p: [float, list, ndarray] electron spectral index (stay within p=2.2-2.5)
    nu: [float, list, ndarray] observing frequency (in Hz)
    b0ref: [float, list, ndarray] magnetic field at density n0=1 (in uG)

    """

    rs = np.asarray(rs)
    vs = np.asarray(vs)
    n0 = np.asarray(n0)
    epse = np.asarray(epse)
    p = np.asarray(p)
    nu = np.asarray(nu)
    b0ref = np.asarray(b0ref)

    #~~~~~~ Constants ~~~~~~~~~~~~~~~#
    ef = 0.38  #Emission filling factor
    c1 = 6.27e18 #Pacholyzyk (1980) constants. Note that these depend on p.
    c5 = 9.68e-24 
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
    B0 = (b0ref*1.0e-6)*n0**0.47 #Eq A4 #Note for SKS: Correct this in the paper
    vA = B0/np.sqrt(4.*3.14*rho0) #Alfven Velocity term (Sec A2, 2nd para)
    MA = vs_cgs/vA #Alfven Mach Number
    
    epsCR = 0.01*(0.15*MA + 6.)  #Cosmic ray acceleration efficiency
    epsCR = np.clip(epsCR,None,0.1) #See Section A2
    epsub = (epsCR/2.)*((vs_cgs/c) + (1./MA)) #Eq A7
    Bu = np.sqrt(8.*3.14*epsub*rho0*vs_cgs**2) #Eq A3
    B = Bu*np.sqrt((1. + 2.*(eta**2))/3.) #Eq A8
    B = np.clip(B,4.*B0,None) #We assert that SNR magnetic field should at least be a simple compression of the ISM field (B0)

    #~~~~~~ Luminosity ~~~~~~~~~~~~~~~~~~~#
    s = ef*(4./3.)*(rs_cgs) #Eq A9
    nu1 = 2.*c1*((s*c6*N_0)**(2./(p+4)))*(B**((p+2)/(p+4))) #Eq A11
    
    L1 = 4.0*(3.14**2)*(rs_cgs**2)*(B**(-0.5))*(c5/c6) #Breaking Eq A10 into 3 parts - L1, L2, L3 - for simplicity
    L2 = 1 - np.exp(-(nu/nu1)**(-(p+4)/2.))

    #This is a short step I am adding because the above equation for L2 runs into numerical errors
    #for nu>>nu1 (i.e. much later in the SNR age). This is particularly apparent for very low densities (e.g. n0<0.01)
    #To avoid this for now, I enforce the approximation: e^-x = 1 - x (for x<<1). May revisit this at some point
    if np.isscalar(L2):
        if nu>20.*nu1:
            L2 = ((nu/nu1)**(-1.0*(p+4)/2.))
    else:
        L2[nu>20.*nu1] = ((nu/nu1)**(-1.0*(p+4)/2.))[nu>20.*nu1]
    L3 = (nu/(2.*c1))**(5./2.)

    return L1*L2*L3 #Eq A10

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#


