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

def time_of_flight(M, N, q, sigt, v, de, tau, option):
    # -----------------------------------------------------------------------------
    # Inputs:
    #   M ......... number of angles for Gauss-Legendre quadrature
    #   N ......... number of time steps
    #   q[g] ...... source term (n/cm^2-s) for group g @ left boundary (x=0cm)  
    #   sigt[g] ... pure absorber cross section (1/cm) for group g
    #   v[g] ...... speed (cm/s) for group g
    #   de[g] ..... energy width for group g
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
        print "\n\n  m =", m, "mu =", mu[m] 
        for g in range(G): 
            if g % (G/20) == 0:
                print "  "+str(int(float(g)/float(G)*100.))+"%"

            # time at which pulse first reaches right boundary (x=1cm)
            pulse_start_time = 1./(mu[m]*v[g]) 

            # time at which neutrons are no longer at right boundary (x=1cm)
            pulse_end_time = 1./(mu[m]*v[g]) + tau

            #num_neutrons[get_index(t,pulse_start_time)-1: get_index(t,pulse_end_time)] += de[g]*w[m]*q[g]*np.exp(-sigt[g]/mu[m])
            #num_neutrons[get_index(t,pulse_start_time)-1: get_index(t,pulse_end_time)] += w[m]*q[g]*np.exp(-sigt[g]/mu[m])

            #for n in range(get_index(t,pulse_start_time)-1, get_index(t,pulse_end_time)+1):
            for n in range(N):
                delta_t = max(0, min(pulse_end_time, t[n]+0.5*dt[n]) - max(pulse_start_time, t[n]-0.5*dt[n]) )
                num_neutrons[n] += delta_t*de[g]*w[m]*q[g]*np.exp(-sigt[g]/mu[m]) # (n/cm^2)

                
    print "Time of flight calculation completed! \n"

    print v
                    
    return t, num_neutrons
###############################################################################
def time_of_flight_subelements(M, N, q, sigt, v, de_sub, tau, option):
    # -----------------------------------------------------------------------------
    # Inputs:
    #   M ......... number of angles for Gauss-Legendre quadrature
    #   N ......... number of time steps
    #   q[g] ...... source term (n/cm^2-s) for group g @ left boundary (x=0cm)  
    #   sigt[g] ... pure absorber cross section (1/cm) for group g
    #   v[g] ...... speed (cm/s) for group g
    #   de[g] ..... energy width for group g
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
    S = len(v)
    num_neutrons = np.zeros(N)
    
    for m in range(M/2,M):
        print "\n\n  m =", m, "mu =", mu[m] 
        for s in range(S): 
            if s % (S/20) == 0:
                print "  "+str(int(float(s)/float(S)*100.))+"%"
            
            g = int(de_sub[s,0])
            
            # time at which pulse first reaches right boundary (x=1cm)
            pulse_start_time = 1./(mu[m]*v[s]) 

            # time at which neutrons are no longer at right boundary (x=1cm)
            pulse_end_time = 1./(mu[m]*v[s]) + tau

            #for n in range(get_index(t,pulse_start_time)-1, get_index(t,pulse_end_time)+1):
            for n in range(N):
                delta_t = max(0, min(pulse_end_time, t[n]+0.5*dt[n]) - max(pulse_start_time, t[n]-0.5*dt[n]) )
                num_neutrons[n] += delta_t*de_sub[s,1]*w[m]*q[s]*np.exp(-sigt[g]/mu[m]) # (n/cm^2)
                

    print "Time of flight calculation completed! \n"

    print v
                    
    return t, num_neutrons
'''###############################################################################
# fission spectrum for U-235
chi = lambda E: 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

# obtain velocity for a particular energy in (cm/s)
vel = lambda E: np.sqrt(2.*E/938.280)*3e10 

# energy at midpoint of energy group (MeV)
emid = np.loadtxt('xs/xs_old/emid', skiprows=1) 
# width of energy group (MeV)
de = np.loadtxt('xs/xs_old/de', skiprows=1) 
# speed of energy group (cm/s)
spgrp = 1e6 * np.loadtxt('xs/xs_old/spgrp', skiprows=1) 
# cross section (1/cm)
sigt = np.loadtxt('xs/xs_old/92235/sigt', skiprows=1) 

M = 32
N = 250
q = chi(emid)/2.
v = spgrp
tau = 1e-10
option = 'log'

(t,mg) = time_of_flight(M, N, q, sigt, v, de, tau, option)
plt.loglog(t,mg)#,'.')
'''###############################################################################
# fission spectrum for U-235
chi = lambda E: 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

# obtain velocity for a particular energy in (cm/s)
vel = lambda E: np.sqrt(2.*E/938.280)*3e10 

# energy at midpoint of energy group (MeV)
emid = np.loadtxt('xs/xs_FEDS_100g_200e/emid', skiprows=1) 
# width of energy group (MeV)
de = np.loadtxt('xs/xs_FEDS_100g_200e/de', skiprows=1) 
# speed of energy group (cm/s)
spgrp = 1e6 * np.loadtxt('xs/xs_FEDS_100g_200e/spgrp', skiprows=1) 
# cross section (1/cm)
sigt = np.loadtxt('xs/xs_FEDS_100g_200e/92235/sigt', skiprows=1) 

M = 32
N = 250
q = chi(emid)/2.
v = spgrp
tau = 1e-10
option = 'log'

(t,feds) = time_of_flight(M, N, q, sigt, v, de, tau, option)
plt.loglog(t,feds)#,'.')\
#plt.show()

###############################################################################
clust = np.loadtxt('xs/xs_FEDS_100g_200e/clust', skiprows=2, usecols=[0,1])
clust[:,1] *= 1e-6 #converting values from clust file from eV to MeV
num_subelements = len(clust) - 1

de_sub = np.zeros((len(clust)-1,2)) 
emid_sub = np.zeros(len(clust)-1) 
de_elem = np.zeros((max(clust[:,0])+1,2)) 
for i in range(len(de_sub)):
    de_sub[i,0] = clust[i,0]
    de_sub[i,1] = clust[i,1] - clust[i+1,1] 
    emid_sub[i] = 0.5*(clust[i,1] + clust[i+1,1]) 
    de_elem[int(clust[i,0]), 0] = clust[i,0]
    de_elem[int(clust[i,0]), 1] += de_sub[i,1]
    print de_sub[i,1], de_elem[int(de_sub[i,0]), 1]

# obtaining the subelement mapping from the custom group structure file
# (the 1st column is the subelement, the 2nd column is the corresponding element)
subelement_mapping = np.zeros((2,num_subelements))
subelement_mapping[0] = [i for i in range(num_subelements)]
subelement_mapping[1] = [clust[i,0] for i in range(num_subelements)] 
num_elements = int(max(subelement_mapping[1])) + 1

map = np.zeros(( num_elements, num_subelements ))
for s in range(num_subelements):
    i = int(subelement_mapping[1,s])
    j = int(s)  
    # map[i,j] = de_subelement/de_element from clust-*.txt file 
    #map[i,j] = de_sub[j,1]/de_elem[int(de_sub[s,0]),1]   
    map[i,j] = 1.
#print map
#np.savetxt('map.txt', map, fmt='%.1f', delimiter=' ')

M = 32
N = 250
q = chi(emid_sub)/2.
v = vel(emid_sub)
tau = 1e-10
option = 'log'

(t,feds) = time_of_flight_subelements(M, N, q, sigt, v, de_sub, tau, option)     
plt.loglog(t,feds)#,'.')\
plt.show()

#print np.size(feds*map)
#print np.shape(feds*map)

'''###############################################################################
# Open cross-section files
sigma_235_t = np.genfromtxt('xs/xs_old/raw_endf_92235/u235_total.csv', delimiter=",")

# Make interpolation functions
sig_235_t_interp = interpolate.interp1d(sigma_235_t[:,0], sigma_235_t[:,1], bounds_error=False, fill_value=1.)

# Create energy grid
e = np.logspace(np.log10(20), -11, 100000)
de = e[0:99999] - e[1:100000]
energies = 0.5*(e[0:99999] + e[1:100000])

M = 32
N = 1000
q = chi(energies)/2.
sigt = sig_235_t_interp(energies)
v = vel(energies)
tau = 1e-10
option = 'log'

(t,exact) = time_of_flight(M, N, q, sigt, v, de, tau, option)
plt.loglog(t,exact)#,'.')
exact_matrix = np.zeros((2,len(t)))
exact_matrix[0] = t
exact_matrix[1] = exact
np.savetxt('exact.txt', exact_matrix.transpose(), fmt='%.18e', delimiter=' ')
###############################################################################
'''

#exact_matrix = np.loadtxt('exact_M32_N1000.txt', delimiter=' ').transpose()
#t = exact_matrix[0]
#exact = exact_matrix[1]
#plt.loglog(t,exact)#,'.')
#plt.show()

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

#fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(energies, chi(energies))
plt.loglog(emid, chi(emid))
plt.ylabel("$\chi$")
plt.xlabel("E (MeV)")
plt.savefig("U-235_chi.pdf")
#plt.show()
plt.close() 

#fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(energies, sig_235_t_interp(energies))
plt.loglog(emid, sig_235_t_interp(emid))
plt.ylabel("$\sigma$ (barns)")
plt.xlabel("E (MeV)")
plt.savefig("U-235_xsect.pdf")
#plt.show()
plt.close() 
'''