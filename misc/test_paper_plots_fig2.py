#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Checking Paper Figure #2 with the new light curve code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# In[3]:
import numpy as np
import matplotlib.pyplot as plt
import Sarba17_SNRLC_oneisentwhite as snrlc
import s17lc
reload(s17lc)

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


print 'CHECKING DENSITIES..'

n0arr = [0.001, 0.1, 10.]
lsarr = ['--', '-', '-.']
e51 = 0.5
mej = 1.4
p = 2.2
epse = 5.0e-3
nu = 1.4e9 


plt.figure(figsize=(7,6))
plt.title('Mej = {0}, e51 = {1}, p = {2}'.format(mej, e51, p), fontsize=12)
plt.plot(diam_knownSNRs, lum_knownSNRs, 'ro')
for j, n0 in enumerate(n0arr):
    tim, rad, lum, vel = snrlc.lightcurve_Full(n0=n0, mej=mej, e51=e51, epse=epse, pp=p, sntype='ia', nu=nu)
    rs = np.zeros(tim.size)
    vs = np.zeros(tim.size)

    for i in range(tim.size):
        rs[i], vs[i] = s17lc.radius_velocity(tim[i], n0, mej, e51, sntype='ia')

    l = s17lc.luminosity(rs, vs, n0, epse, p, nu)
    l[vs<200.]=0 #this is from the assumption in the paper that after shock decelerates to 200 km/s, 
              #electrons aren't accelerated efficiently, so no radio emission. 

    plt.plot(2.*rad[lum>0], lum[lum>0], color='k', ls=lsarr[j], label='n0 = {}'.format(n0))
    plt.plot(2.*rs[l>0], l[l>0], color='k', ls='-', lw=5., alpha=0.4, label='new Code n0 = {}'.format(n0))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Diameter [pc]', fontsize=15)
plt.ylabel('Luminosity [ergs/s]', fontsize=15)
plt.xlim(0.5, 500.)
plt.ylim(2.0e19, 4.0e26)
plt.legend()
plt.savefig(outp + 'check_paper_plot_Figure3_density.png', bbox_inches='tight')

print 'CHECKING ENERGIES..'

n0 = 0.5
lsarr = ['--', '-', '-.']
e51arr = [2., 0.5, 0.08]
mej = 5.
p = 2.2
epse = 5.0e-3
nu = 1.4e9 


plt.figure(figsize=(7,6))
plt.title('Mej = {0}, n0 = {1}, p = {2}'.format(mej, n0, p), fontsize=12)
plt.plot(diam_knownSNRs, lum_knownSNRs, 'ro')
for j, e51 in enumerate(e51arr):
    tim, rad, lum, vel = snrlc.lightcurve_Full(n0=n0, mej=mej, e51=e51, epse=epse, pp=p, sntype='cc', nu=nu)
    rs = np.zeros(tim.size)
    vs = np.zeros(tim.size)

    for i in range(tim.size):
        rs[i], vs[i] = s17lc.radius_velocity(tim[i], n0, mej, e51, sntype='cc')

    l = s17lc.luminosity(rs, vs, n0, epse, p, nu)
    l[vs<200.]=0 #this is from the assumption in the paper that after shock decelerates to 200 km/s, 
              #electrons aren't accelerated efficiently, so no radio emission. 

    plt.plot(2.*rad[lum>0], lum[lum>0], color='k', ls=lsarr[j], label='E = {}'.format(e51))
    plt.plot(2.*rs[l>0], l[l>0], color='k', ls='-', lw=5., alpha=0.4, label='new Code E = {}'.format(e51))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Diameter [pc]', fontsize=15)
plt.ylabel('Luminosity [ergs/s]', fontsize=15)
plt.xlim(0.5, 500.)
plt.ylim(2.0e19, 4.0e26)
plt.legend()
plt.savefig(outp + 'check_paper_plot_Figure3_energy.png', bbox_inches='tight')

print 'CHECKING MASSES...'

n0 = 0.1
lsarr = ['--', '-', '-.']
e51 = 1.
mejarr= [2., 10., 20.]
p = 2.2
epse = 5.0e-3
nu = 1.4e9 


plt.figure(figsize=(7,6))
plt.title('e51 = {0}, n0 = {1}, p = {2}'.format(e51, n0, p), fontsize=12)
plt.plot(diam_knownSNRs, lum_knownSNRs, 'ro')
for j, mej in enumerate(mejarr):
    tim, rad, lum, vel = snrlc.lightcurve_Full(n0=n0, mej=mej, e51=e51, epse=epse, pp=p, sntype='cc', nu=nu)
    rs = np.zeros(tim.size)
    vs = np.zeros(tim.size)

    for i in range(tim.size):
        rs[i], vs[i] = s17lc.radius_velocity(tim[i], n0, mej, e51, sntype='cc')

    l = s17lc.luminosity(rs, vs, n0, epse, p, nu)
    l[vs<200.]=0 #this is from the assumption in the paper that after shock decelerates to 200 km/s, 
              #electrons aren't accelerated efficiently, so no radio emission. 

    plt.plot(2.*rad[lum>0], lum[lum>0], color='k', ls=lsarr[j], label='M = {}'.format(mej))
    plt.plot(2.*rs[l>0], l[l>0], color='k', ls='-', lw=5., alpha=0.4, label='new Code M = {}'.format(mej))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Diameter [pc]', fontsize=15)
plt.ylabel('Luminosity [ergs/s]', fontsize=15)
plt.xlim(0.5, 500.)
plt.ylim(2.0e19, 4.0e26)
plt.legend()
plt.savefig(outp + 'check_paper_plot_Figure3_Mass.png', bbox_inches='tight')


