#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# This contains radius/velocity/luminosity in the ejecta
# and Sedov phases as calculated 

import numpy as np
import matplotlib.pyplot as plt

#<><><><><><><><><> F U N C T I O N S <><><><><><><><><><><><>#

def r_ej_Ia(t2, n0, e51, mej):
    """
    Ejecta dominated radius of Ia in parsec
    """

    return 1.29*(t2**0.7)*(e51**0.35)*(n0**(-0.1))*(mej**(-0.25))

def r_ej_cc(t2, n0, e51, mej):
    """
    Ejecta dominated radius of CC in parsec
    """

    return 1.26*(t2**0.75)*(e51**0.38)*(n0**(-0.08))*(mej**(-0.29))

def r_st(t4, e51, n0):
    """
    Sedov-Taylor radius in parsec
    """

    return 12.5*(t4**0.4)*(e51**(0.2))*(n0**(-0.2))

def v_ej_Ia(t2, n0, e51, mej):
    """
    Ejecta dominated velocity of Ia in km/s
    """
    
    return 8797.*(t2**(-0.3))*(e51**0.35)*(n0**(-0.1))*(mej**(-0.25))

def v_ej_cc(t2, n0, e51, mej):
    """
    Ejecta dominated velocity of CC in km/s
    """

    return 9213.*(t2**(-0.25))*(e51**0.38)*(n0**(-0.08))*(mej**(-0.29))

def v_st(t4, e51, n0):
    """
    Sedov-Taylor velocity in km/s
    """
    return 490.*(t4**(-0.6))*(e51**0.2)*(n0**(-0.2))

def t_ejtosedov(n0, e51, mej, sntype='cc'):
    """
    Age [in years] at which SN transitions from Ejecta dominated to Sedov phase
    """
    ts = 0.481 if sntype=='ia' else 0.424 
    return (423.*ts)*(e51**(-0.5))*(mej**(5./6.))*(n0**(-1./3.))

def luminosity(rs, vs, n0, e51, mej, epse, p, nu):
    """
    Since this is just for test purposes, I am only calculating luminosities
    for core-collapse supernovae, not Type Ia. Change it accordingly. 
    """
    #~~~~~~ Constants ~~~~~~~~~~~~~~~#
    ef = 0.38
    c1 = 6.27e18
    c5 = 9.68e-24 #Correction in paper
    c6 = 8.1e-41
    mu = 1.4
    mp = 1.67e-24 #g
    rho0 = mu*mp*n0
    me =  9.11e-28  #in grams                                                                                                                   
    c = 3.0e10 #cm/s
    eta = 4.
    Em = me*c**2
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~ CGS units for r,v  ~~~~~~~~~~~~~~~#
    rs_cgs = rs*3.086e18 #cm
    vs_cgs = vs*1.0e5 #cm/s

    #~~~~~~ Electron Acceleration Stuff (Eq A2) ~~~~~~~~~~~~~~~#
    N_0 = (p-2.)*epse*rho0*(vs_cgs**2)*(Em**(p-2))

    #~~~~~~ Magnetic Field Stuff~~~~~~~~~~~~~~~#
    B0 = 9.0e-6*n0**0.47 #Eq A4 #Note for SKS: Correct this in the paper
    vA = B0/np.sqrt(4.*3.14*rho0) #Alfven Velocity term (Sec A2, 2nd para)
    MA = vs_cgs/vA #Mach Number
    
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
    L3 = (nu/(2.*c1))**(5./2.)

    return L1*L2*L3 #Eq A10
    
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
