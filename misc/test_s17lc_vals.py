#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Iterating over all values of density, energy, ejecta
# mass to make sure s17lc.py is consistent with the 
# light curve models I have used before
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# In[3]:
import numpy as np
import matplotlib.pyplot as plt
import Sarba17_SNRLC_oneisentwhite as snrlc
import s17lc
import s17paperlc as s17pap

#Output path
outp = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/newRadSNRs/diagnostics/erratum_checks/luminosity/'

# In[15]:
#Case and Bhattacharya SNRs
path = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/RADSNRS/Inputs/'
table = np.genfromtxt(path+'KnownSNRS_CaseBhat98.dat',skip_header=19,usecols=(2,3,4))  
sb_knownSNRs = table[:,0]
diam_knownSNRs = table[:,2]
#Column 2: SB, Col3: Distance, Col4: Diameter
pc = 3.086e16 #m
joule = 1.0e7 #ergs
lum_knownSNRs = sb_knownSNRs*(3.14/4.0)*((diam_knownSNRs*pc)**2.0)*joule


# In[20]:

n0arr = [0.01, 1., 10.]
e51arr = [0.1, 1., 10.]
mejarr = [0.6, 1., 1.4]
p = 2.2
nu = 1.4e9
epse = 0.0042

for n0 in n0arr:
    for e51 in e51arr:
        for mej in mejarr:
#From my original code
            print 'n0={0}, e51={1}, mej={2}'.format(n0, e51, mej)
            #FROM Light Curve I sent Rick
            #tim, rad, lum, vel = snrlc.lightcurve_Full(n0=n0, mej=mej, e51=e51, epse=epse, pp=p, sntype='cc', nu=nu)
            
            #FROM Light curve built from all the equations as they appear in Sarbadhicary et al (2017)
            trans = s17pap.t_ejtosedov(n0, e51, mej, sntype='ia')
            tej = np.linspace(1, trans, 400)
            tsedov = np.linspace(trans, 1.0e5, 400)
            rej = s17pap.r_ej_Ia(tej*1.0e-2, n0, e51, mej)
            rst = s17pap.r_st(tsedov*1.0e-4, e51, n0)
            vej = s17pap.v_ej_Ia(tej*1.0e-2, n0, e51, mej)
            vst = s17pap.v_st(tsedov*1.0e-4, e51, n0)
            lej = s17pap.luminosity(rej, vej, n0, e51, mej, epse, p, nu)
            lst = s17pap.luminosity(rst, vst, n0, e51, mej, epse, p, nu)

            #FROM s17lc.py
            tim = np.logspace(0, 5, 800)
            rs = np.zeros(tim.size)
            vs = np.zeros(tim.size)
            for i in range(tim.size):
                rs[i], vs[i] = s17lc.radius_velocity(tim[i], n0, mej, e51, sntype='ia')
            ls = s17lc.luminosity(rs, vs, n0, epse, p, nu)
            ls[vs<200.]=0 #this is from the assumption in the paper that after shock decelerates to 200 km/s, 
              #electrons aren't accelerated efficiently, so no radio emission. 

            plt.figure(figsize=(14,6))
            plt.subplot(1,2,1)
            plt.title(r'n0={0}, e51={1}, mej={2}'.format(n0, e51, mej), fontsize=15)
            #plt.plot(tim[lum>0], lum[lum>0], 'r-', lw=5., alpha=0.4, label='Code I sent')
            plt.plot(tej, lej, 'r-', lw=5., alpha=0.5, label='Eqs from paper')
            plt.plot(tsedov, lst, 'r-', lw=5., alpha=0.5)
            plt.plot(tim[ls>0], ls[ls>0], 'k-', label='New Code')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Age [years]', fontsize=15)
            plt.ylabel('Luminosity [ergs/s]', fontsize=15)
            plt.xlim(1., 1.0e6)
            plt.ylim(1.0e21, 1.0e27)
            plt.legend()

            plt.subplot(1,2,2)
            plt.title(r'n0={0}, e51={1}, mej={2}'.format(n0, e51, mej), fontsize=15)
#            plt.plot(2.*rad[lum>0], lum[lum>0], 'r-', lw=5., alpha=0.4, label='Code I sent')
            plt.plot(2.*rej, lej, 'r-', lw=5., alpha=0.5, label='Eqs from paper')
            plt.plot(2.*rst, lst, 'r-', lw=5., alpha=0.5)
            plt.plot(2.*rs[ls>0], ls[ls>0], 'k-', label='New Code')
            plt.plot(diam_knownSNRs, lum_knownSNRs, 'ro', label='CB98')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Diameter [pc]', fontsize=15)
            plt.ylabel('Luminosity [ergs/s]', fontsize=15)
            plt.xlim(1., 200.)
            plt.ylim(1.0e21, 1.0e27)
            plt.legend()
            plt.savefig(outp + 'Ia_paper_lum_age_diam_n0_{0}_e51_{1}_mej_{2}.png'.format(n0, e51, mej), bbox_inches='tight')




