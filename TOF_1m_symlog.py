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

def time_of_flight(mu, I, N, q, sigt_Fe, sigt_U, v, de, tau, option):
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
        detector_end_time = 4e-7
        times = np.logspace(np.log10(detector_start_time), np.log10(detector_end_time), num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    elif option == 'lin':
        detector_start_time = 0
        detector_end_time = 4e-7
        times = np.linspace(detector_start_time, detector_end_time, num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    else:
        print " ** Fatal Error: option must be 'lin' or 'log'" 

    if len(sigt_U) != len(sigt_Fe):
        print "Fatal Error: len(sigt_U) != len(sigt_Fe)"
        return 
    
    G = len(sigt_U)
    num_neutrons = np.zeros(N)
    x_cells = np.linspace(0, 0.05, num=I+1)
    dx = x_cells[1:I+1] - x_cells[0:I]
    x_source =  0.5*(x_cells[1:I+1] + x_cells[0:I])
    
    #GL_points, GL_weights = leggauss(I)
    #x_source = (GL_points + 1)*0.05/2.
    #dx = GL_weights*0.05/2.
    
    print x_source, dx
    
    for i in range(I):
        print "\n\n  i =", i, "x_source =", x_source[i]
        
        for g in range(G): 
            if g % (G/20) == 0:
                print "  "+str(int(float(g)/float(G)*100.))+"%"
    
            # time at which pulse first reaches right boundary (x=1cm)
            pulse_start_time = (100.-x_source[i])/(mu*v[g]) 
    
            # time at which neutrons are no longer at right boundary (x=1cm)
            pulse_end_time = (100.-x_source[i])/(mu*v[g]) + tau
    
            N_start = max(0, -1+np.searchsorted(t, pulse_start_time, side='left'))
            N_end = min(N, 1+np.searchsorted(t, pulse_end_time, side='right'))
            for n in range(N_start, N_end):
                #if t[n]-0.5*dt[n] < pulse_end_time and t[n]+0.5*dt[n] > pulse_start_time:  
                delta_t = max(0, min(pulse_end_time, t[n]+0.5*dt[n]) - max(pulse_start_time, t[n]-0.5*dt[n]) )
                num_neutrons[n] += delta_t*dx[i]*de[g]*q[g]*np.exp(-sigt_U[g]*(0.05-x_source[i]))*np.exp(-2.06*sigt_Fe[g]/mu) # (n/cm^2)
                
    print "Time of flight calculation completed! \n"

    print v
                    
    return t, num_neutrons
###############################################################################
def time_of_flight_subelements(mu, I, N, q, sigt_Fe, sigt_U, v, de_sub, tau, option):
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
        detector_end_time = 4e-7
        times = np.logspace(np.log10(detector_start_time), np.log10(detector_end_time), num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    elif option == 'lin':
        detector_start_time = 0
        detector_end_time = 4e-7
        times = np.linspace(detector_start_time, detector_end_time, num=N+1)
        dt = times[1:N+1] - times[0:N]
        t = 0.5*(times[1:N+1] + times[0:N])
    else:
        print " ** Fatal Error: option must be 'lin' or 'log'" 

    S = len(v)
    num_neutrons = np.zeros(N)
    
    x_cells = np.linspace(0, 0.05, num=I+1)
    dx = x_cells[1:I+1] - x_cells[0:I]
    x_source =  0.5*(x_cells[1:I+1] + x_cells[0:I])
    
    for i in range(I):
        print "\n\n  i =", i, "x_source =", x_source[i]
    
        for s in range(S): 
            if s % (S/20) == 0:
                print "  "+str(int(float(s)/float(S)*100.))+"%"
            
            g = int(de_sub[s,0])

            # time at which pulse first reaches right boundary (x=1cm)
            pulse_start_time = (100.-x_source[i])/(mu*v[s]) 
    
            # time at which neutrons are no longer at right boundary (x=1cm)
            pulse_end_time = (100.-x_source[i])/(mu*v[s]) + tau

            N_start = max(0, -1+np.searchsorted(t, pulse_start_time, side='left'))
            N_end = min(N, 1+np.searchsorted(t, pulse_end_time, side='right'))
            for n in range(N_start, N_end):
                delta_t = max(0, min(pulse_end_time, t[n]+0.5*dt[n]) - max(pulse_start_time, t[n]-0.5*dt[n]) )
                num_neutrons[n] += delta_t*dx[i]*de_sub[s,1]*q[s]*np.exp(-sigt_U[g]*(0.05-x_source[i]))*np.exp(-2.06*sigt_Fe[g]/mu) # (n/cm^2)
                
    print "Time of flight calculation completed! \n"

    print v
                    
    return t, num_neutrons
###############################################################################
# fission spectrum for U-235
chi = lambda E: 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

# obtain velocity for a particular energy in (cm/s)
vel = lambda E: np.sqrt(2.*E/938.280)*3e10 

# energy at midpoint of energy group (MeV)
emid = np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/emid', skiprows=1) 
# width of energy group (MeV)
de = np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/de', skiprows=1) 
# speed of energy group (cm/s)
spgrp = 1e6 * np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/spgrp', skiprows=1) 
# cross section (1/cm)
sigt_Fe = np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/26000/sigt', skiprows=1) 
sigt_U = np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/92235/sigt', skiprows=1) 

#sigt_Fe *= 7.874*0.6022/55.847
#sigt_U  *= 18.39374*0.6022/235.

mu = 1.
I = 20
N = 250
q = chi(emid)/2.
v = spgrp
#v = vel(emid)
tau = 1e-10
option = 'lin'

(t,feds) = time_of_flight(mu, I, N, q, sigt_Fe, sigt_U, v, de, tau, option)
plt.semilogy(t,feds+1e-21*np.ones(len(feds)),'b')
###############################################################################
clust = np.loadtxt('xs/xs_FEDS_1t_200r_600b_1000u/clust', skiprows=2, usecols=[0,1])
clust[:,1] *= 1e-6 #converting values from clust file from eV to MeV
num_subelements = len(clust) - 1

de_sub = np.zeros((len(clust)-1,2)) 
emid_sub = np.zeros(len(clust)-1) 
for i in range(len(de_sub)):
    de_sub[i,0] = clust[i,0]
    de_sub[i,1] = clust[i,1] - clust[i+1,1] 
    emid_sub[i] = 0.5*(clust[i,1] + clust[i+1,1]) 

mu = 1.
I = 20
N = 250
q = chi(emid_sub)/2.
v = vel(emid_sub)
tau = 1e-10
option = 'lin'

(t,subelements) = time_of_flight_subelements(mu, I, N, q, sigt_Fe, sigt_U, v, de_sub, tau, option)     
plt.semilogy(t,subelements+1e-21*np.ones(len(subelements)),'r')
plt.show()
##############################################################################
# Open cross-section files
sigma_Fe_t = np.loadtxt('xs/xs_NJOY_pointwise/26000', skiprows=1, usecols=[0,1])
sigma_U_t = np.loadtxt('xs/xs_NJOY_pointwise/92235', skiprows=1, usecols=[0,1])

sigma_Fe_t[:,0] *= 1e-6
sigma_U_t[:,0]  *= 1e-6
#sigma_Fe_t[:,1] *= 7.874*0.6022/55.847
#sigma_U_t[:,1]  *= 18.39374*0.6022/235.

# Make interpolation functions
sigma_Fe_t_interp = interpolate.interp1d(sigma_Fe_t[:,0], sigma_Fe_t[:,1], bounds_error=False, fill_value=1.)
sigma_U_t_interp = interpolate.interp1d(sigma_U_t[:,0], sigma_U_t[:,1], bounds_error=False, fill_value=1.)

# Create energy grid
G_njoy = 10000
#G_njoy = 100000
#e = np.logspace(np.log10(5e-1), np.log10(1e-2), G_njoy)
e = np.logspace(np.log10(14.2), np.log10(1e-2), G_njoy)
de = e[0:G_njoy-1] - e[1:G_njoy]
emid = 0.5*(e[0:G_njoy-1] + e[1:G_njoy])

mu = 1.
I = 20
N = 250
#I = 20
#N = 10000
sigt_Fe = sigma_Fe_t_interp(emid)
sigt_U = sigma_U_t_interp(emid)
q = chi(emid)/2.
v = vel(emid)
tau = 1e-10
#option = 'log'
option = 'lin'

(t,exact) = time_of_flight(mu, I, N, q, sigt_Fe, sigt_U, v, de, tau, option)
L2_error_feds = np.linalg.norm(exact-feds)/np.linalg.norm(exact)
L2_error_subelements = np.linalg.norm(exact-subelements)/np.linalg.norm(exact)
print L2_error_feds
print L2_error_subelements
plt.semilogy(t,exact+1e-21*np.ones(len(exact)),'g')#,'.')
plt.show()
plt.close()
exact_matrix = np.zeros((2,len(t)))
exact_matrix[0] = t
exact_matrix[1] = exact
np.savetxt('exact.txt', exact_matrix.transpose(), fmt='%.18e', delimiter=' ')
'''###############################################################################
exact_matrix = np.loadtxt('exact.txt', delimiter=' ').transpose()
t = exact_matrix[0]
exact = exact_matrix[1]
plt.plot(t,exact,'purple')
plt.yscale('log')
#plt.semilogy(t,exact,'purple')
plt.xlim([0, 4e-7])
plt.ylim([1e-23, 1e-11])
plt.show()
'''
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
'''
'''
#fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(emid, chi(emid))
plt.ylabel("$\chi$")
plt.xlabel("E (MeV)")
plt.xlim([1e-6, 10])
plt.savefig("U-235_chi.pdf")
plt.show()
plt.close() 

#fig = plt.figure(figsize=(8,6), dpi=1600)
plt.loglog(emid, sigma_Fe_t_interp(emid))
plt.loglog(emid, sigma_U_t_interp(emid))
plt.ylabel("$\sigma (\frac{1}{\text{cm}})$")
plt.xlabel("E (MeV)")
plt.xlim([1e-6, 10])
plt.savefig("U-235_xsect.pdf")
plt.show()
plt.close() 
'''
