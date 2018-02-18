# -----------------------------------------------------------------------------
# Problem Description:
#   - 1 cm thick slab
#   - pure absorber (sigma = total cross section of U-235)
#   - isotropic fission source on the left boundary, which turns off after tau seconds
#     (phi(x=0,E,t) = chi(E) for t=0 to t=tau, where chi is the fission spectrum for U-235)
#   
# Solution:
#   - time-of-flight signal (n/cm^2) at the right boundary (x = 1cm)  
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from numpy.polynomial.legendre import leggauss

def get_index(array, value):
    index = len(array) - 1
    for i in range(len(array)):
        if array[i] >= value:
            index = i
            break
    return index

def time_of_flight(M, N, q, sigt, v, tau, option):
    # -----------------------------------------------------------------------------
    # Inputs:
    #   M ......... number of angles for Gauss-Legendre quadrature
    #   N ......... number of time steps
    #   q[g] ...... source term (n/cm^2-s) for group g @ left boundary (x=0cm)  
    #   sigt[g] ... pure absorber cross section (1/cm) for group g
    #   v[g] ...... speed (cm/s) for group g
    #   tau ....... lifespan of pulse (s), pulse lasts from t=0 to t=tau
    #   option .... 'lin' (linear) or 'log' (logarithmic) time bins for time-of-flight calculation 
    #
    # Outputs:
    #   t[n] .............. midpoint of time bin n
    #   num_neutrons[n] ... number of neutrons @ right boundary (x=1cm) for time bin n
    # -----------------------------------------------------------------------------

    print "Starting time of flight calculation..."

    if option == 'log':
        detector_start_time = 1e-10
        detector_end_time = 1e-6
        times = np.logspace(np.log10(detector_start_time), np.log10(detector_end_time), num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    elif option == 'lin':
        detector_start_time = 0
        detector_end_time = 1e-6
        times = np.linspace(detector_start_time, detector_end_time, num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    else:
        print " ** Fatal Error: option must be 'lin' or 'log'" 

    mu, w = leggauss(M)
    G = len(sigt)
    num_neutrons = np.zeros(N)
    
    for m in range(M/2,M):
        for g in range(G): 
            if g % (G/10) == 0:
                print "  "+str(int(float(g+1)/float(G)*100.))+"%"

            # time at which pulse first reaches right boundary (x=1cm)
            pulse_start_time = 1/(mu[m]*v[g]) 

            # time at which neutrons are no longer at right boundary (x=1cm)
            pulse_end_time = 1/(mu[m]*v[g]) + tau

            for n in range(get_index(t,pulse_start_time)-1, get_index(t,pulse_end_time)+1):
                #delta_t = dt[n] 
                delta_t = min(pulse_end_time, t[n]+0.5*dt[n]) - max(pulse_start_time, t[n]-0.5*dt[n]) 
                num_neutrons[n] += delta_t*de[g]*w[m]*q[g]*np.exp(-sigt[g]/mu[m])

    print "Time of flight calculation completed! \n"
                    
    return t, num_neutrons

###############################################################################
# fission spectrum for U-235
chi = lambda E: 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

# obtain velocity for a particular energy in (cm/s)
vel = lambda E: np.sqrt(2.*E/938.280)*3e10 

# energy at midpoint of energy group (MeV)
emid = np.loadtxt('xs/emid', skiprows=1) 
# width of energy group (MeV)
de = np.loadtxt('xs/de', skiprows=1) 
# speed of energy group (cm/s)
spgrp = 1e6 * np.loadtxt('xs/spgrp', skiprows=1) 
# cross section (1/cm)
sigt = np.loadtxt('xs/92235/sigt', skiprows=1) 

M = 2
N = 1000000
q = chi(emid)/2.
v = spgrp
tau = 1e-10
option = 'log'

(t,mg) = time_of_flight(M, N, q, sigt, v, tau, option)
plt.loglog(t,mg)
###############################################################################
# Open cross-section files
sigma_235_t = np.genfromtxt('xs/raw_endf_92235/u235_total.csv', delimiter=",")

# Make interpolation functions
sig_235_t_interp = interpolate.interp1d(sigma_235_t[:,0], sigma_235_t[:,1], bounds_error=False, fill_value=sigma_235_t[-1,1])

# Create energy grid
e = np.logspace(2, -11, 100000)
de = e[0:99999] - e[1:100000]
energies = 0.5*(e[0:99999] + e[1:100000])


M = 2
N = 1000
q = chi(energies)/2.
sigt = sig_235_t_interp(energies)
v = vel(energies)
tau = 1e-10
option = 'log'

(t,exact) = time_of_flight(M, N, q, sigt, v, tau, option)
plt.loglog(t,exact)
###############################################################################

plt.show()
plt.close()

'''
#plt.loglog(t,feds2,'-r')
plt.ylabel("$\phi$")
plt.xlabel("$t$ (sec)")
#plt.axis([1e-10, 1e-6, 1e-18, 1e-6])
#plt.savefig("Exact_vs_MG.pdf")
#plt.axis([1e-10, 1e-8, 1e-14, 1e-6])
#plt.savefig("Exact_vs_MG_zoom1.pdf")
#plt.axis([3e-10, 3e-9, 1e-8, 1e-6])
plt.savefig("Exact_vs_MG_zoom2.pdf")
plt.show()

plt.close() 


fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(energies, chi(energies))
plt.loglog(emid, chi(emid))
plt.ylabel("$\chi$")
plt.xlabel("E (MeV)")
plt.savefig("U-235_chi.pdf")
plt.close() 

fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(energies, sig_235_t_interp(energies))
plt.loglog(emid, sigt)
plt.ylabel("$\sigma$ (barns)")
plt.xlabel("E (MeV)")
plt.savefig("U-235_xsect.pdf")
plt.close() 
'''
