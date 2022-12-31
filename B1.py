# -*- coding: utf-8 -*-
"""
@author: HuidobroMG

Description:
    
    In this code I am solving the equation of motion for the B = 1 skyrmion
    in spherical coordinates using a shooting method. The hedgehog ansatz
    leads to a 1-dimensional second order differential equation which can
    be solved using a 4th order Runge-Kutta method. You may consider the 
    different terms in the generalized Skyrme lagrangian.
    The system is solved in Skyrme units:
        E2MeV = 3*np.pi**2*fpi/e
        x2fm = hbarc/(fpi*e)
        
"""

#------------------------------------------------------------------------------

# Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scin
import scipy.optimize as scop

#------------------------------------------------------------------------------

# The code does not work for c2 = c4 = 0 because in that case the soliton
# is a compacton, then f(r > R) must be exactly 0.

# Parameters
c2 = 1.0
c4 = 1.0
c6 = 0.0 # 2*lambda2*fpi**2*e**4/hbarc**3
c0 = 0.0 # 2*mpi**2/(fpi*e)**2

coefs = [c2, c4, c6, c0]

#------------------------------------------------------------------------------

# Baryon density
def B0(r, f, df):
    return -np.sin(f)**2*df/(2*np.pi**2*r**2)

# Energy density
def ED(r, f, df, coefs):
    ED_2 = coefs[0]*(df**2 + 2*np.sin(f)**2/r**2)
    ED_4 = coefs[1]*(8*(np.sin(f)*df/r)**2 + 4*(np.sin(f)/r)**4)
    ED_6 = coefs[2]*(np.sin(f)**2*df/r**2)**2
    ED_0 = coefs[3]*(1 - np.cos(f))
    return (ED_2 + ED_4 + ED_6 + ED_0)/(24*np.pi**2)

# Equation of motion
def EL(r, inits):
    f, fp = inits
    
    df = fp
    ddf = (c2*(-2*r**3*df + r**2*np.sin(2*f)) + 
           c4*(-4*r**2*np.sin(2*f)*df**2 + 4*np.sin(f)**2*np.sin(2*f)) + 
           c6*(-df**2*np.sin(f)**2*np.sin(2*f) + 2*np.sin(f)**4/r*df) + 
           c0/2*r**4*np.sin(f))/(c2*r**4 + 8*c4*r**2*np.sin(f)**2 + 
                                 c6*np.sin(f)**4)

    # Cut the second derivative since it could diverge
    if ddf > 5:
        ddf = 5
    if ddf < -5:
        ddf = -5
    return np.array([df, ddf])

#------------------------------------------------------------------------------

# Grid parameters
r_init = 0.01
r_end = 15.0
dr = 1e-3
r = np.arange(r_init, r_end, dr)
N = len(r)

# Initial conditions
f0 = np.pi
f_tol = 1e-10

# Shooting method
param_max = 10.0
param_min = -8.0
iterations = 200
for i in range(iterations):
    param = (param_max + param_min)/2
    
    vinic = [f0 + param*r_init, param]
    sol = scin.solve_ivp(EL, (r_init, r_end), vinic, method = 'RK45', 
                         t_eval = r, rtol = 1e-12)
    f, df = sol.y[0], sol.y[1]
    
    if abs(f[-1]) < f_tol:
        print('Finished at i = ', i)
        break
    
    if f[-1] > 0:
        param_max = param
    else:
        param_min = param

# Compute the second derivative
ddf = np.zeros(N)
for i in range(N):
    ddf[i] = EL(r[i], [f[i], df[i]])[1]

#------------------------------------------------------------------------------

# Energy and Baryon densities
Bd = B0(r, f, df)
Ed = ED(r, f, df, coefs)
E2 = ED(r, f, df, [coefs[0], 0, 0, 0])
E4 = ED(r, f, df, [0, coefs[1], 0, 0])
E6 = ED(r, f, df, [0, 0, coefs[2], 0])
E0 = ED(r, f, df, [0, 0, 0, coefs[3]])

# Integral with the corresponding volume form
B = scin.simpson(4*np.pi*r**2*Bd, r)

# Energy
E = scin.simpson(4*np.pi*r**2*Ed, r)

# Calculate the mean square radius of each term
rRMS = np.sqrt(scin.simpson(4*np.pi*r**4*Bd, r)/B)

print('B = ', np.round(B, 5))
print('E = ', np.round(E, 5))
print('r_RMS = ', np.round(rRMS, 5))

#------------------------------------------------------------------------------

# Plots of energy density
fig = plt.figure(figsize = (12, 8))
fig.subplots_adjust(top = 0.93, bottom = 0.11, left = 0.07, right = 0.96, 
                    hspace = 0, wspace = 0.14)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(r, f, 'r-')
ax1.plot(r, df, 'b-')
ax1.plot(r, ddf, 'g-')

ax2.plot(r, 4*np.pi*r**2*Ed, 'k-')
ax2.plot(r, 4*np.pi*r**2*E2, 'b-', label = 'E2')
ax2.plot(r, 4*np.pi*r**2*E4, 'g-', label = 'E4')
ax2.plot(r, 4*np.pi*r**2*E6, 'r-', label = 'E6')
ax2.plot(r, 4*np.pi*r**2*E0, 'c-', label = 'E0')

ax2.legend(fontsize = 13)

