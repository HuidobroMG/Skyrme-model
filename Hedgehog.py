"""
@author: HuidobroMG

We find the B = 1 skyrmion solution under spherical symmetry using a shooting method for the value of
the first derivative of the profile function.
Additionally, we compute further observables like the root mean square radius, the pion-nucleon-nucleon coupling
constant and the isospin moment of inertia under the standard collective coordinates quantization method.
Finally, we construct the O(4) representation of the Skyrme field for this solution.
"""

# Import the modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scin
import scipy.optimize as scop
import scipy.interpolate as scinter

# Parameters
hbarc = 197.3269804 # [MeV*fm]
fpi = 108 # Pion decay constant [MeV]
e = 4.84 # Skyrme parameters
lambda2 = 0.0 # Omega-Pion coupling constant [MeV*fm**3]
mpi = 138.0 # Pion mass [MeV]

# Adimensional coupling constants
c2 = 1.0
c4 = 1.0
c6 = 2*lambda2*fpi**2*e**4/hbarc**3
c0 = 2*mpi**2/(fpi*e)**2

coefs = [c2, c4, c6, c0]

# Energy and length scales
E2MeV = 3*np.pi**2*fpi/e # Energy scale
x2fm = hbarc/(fpi*e) # Length scale
hbarc_adim = e**2/(3*np.pi**2) # Quantum scale

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

# Equations of motion
def EL(r, inits):
    f, df = inits
    
    ddf = (c2*(-2*r**3*df + r**2*np.sin(2*f)) + 
           c4*(-4*r**2*np.sin(2*f)*df**2 + 4*np.sin(f)**2*np.sin(2*f)) + 
           c6*(-df**2*np.sin(f)**2*np.sin(2*f) + 2*np.sin(f)**4/r*df) + 
           c0/2*r**4*np.sin(f))/(c2*r**4 + 8*c4*r**2*np.sin(f)**2 + 
                                 c6*np.sin(f)**4)
    
    # Avoid too large derivatives               
    if ddf > 5:
        ddf = 5
    if ddf < -5:
        ddf = -5
    return np.array([df, ddf])

# Profile function decay at large distances
def Decay(x, a):
    return a*(np.sqrt(2.0/c0)/x**2 + 1.0/x)*np.exp(-np.sqrt(c0/2.0)*x)

# Initial conditions
f0 = np.pi

# Space grid
r_init = 1e-4
r_end = 35.0
dr = 1e-3
r = np.arange(r_init, r_end, dr)
N = len(r)

# Tolerance estimated from previous solutions
f_tol = Decay(r_end, 1.9)

# Shooting method
param_max = 10.0
param_min = -8.0
iterations = 100
for i in range(iterations):
    if i%20 == 0:
        print('i_shoot = ', i)
    param = (param_max + param_min)/2
    
    vinic = [f0 + param*r_init, param]
    sol = scin.solve_ivp(EL, (r_init, r_end), vinic, method = 'RK45', 
                         t_eval = r, rtol = 2.5e-14, atol = 1e-10)
    f, df = sol.y[0], sol.y[1]
    
    if abs(f[-1]) < f_tol:
        print('Finished at i = ', i)
        break

    if f[-1] > 0:
        param_max = param
    else:
        param_min = param

ddf = np.zeros(N)
for i in range(N):
    ddf[i] = EL(r[i], [f[i], df[i]])[1]

# Energy and Baryon densities
Bd = B0(r, f, df)
Ed = ED(r, f, df, coefs)
Ed2 = ED(r, f, df, [coefs[0], 0, 0, 0])
Ed4 = ED(r, f, df, [0, coefs[1], 0, 0])
Ed6 = ED(r, f, df, [0, 0, coefs[2], 0])
Ed0 = ED(r, f, df, [0, 0, 0, coefs[3]])

# Integral with the corresponding volume form
B = scin.simpson(4*np.pi*r**2*Bd, x = r)

# Energies
E2 = scin.simpson(4*np.pi*r**2*Ed2, x = r)
E4 = scin.simpson(4*np.pi*r**2*Ed4, x = r)
E6 = scin.simpson(4*np.pi*r**2*Ed6, x = r)
E0 = scin.simpson(4*np.pi*r**2*Ed0, x = r)
E = scin.simpson(4*np.pi*r**2*Ed, x = r)

# Calculate the root mean square radius of each term
rRMS = np.sqrt(scin.simpson(4*np.pi*r**4*Bd, x = r)/B)

# Compute the Pion-Nucleon-Nucleon coupling constant
cutoff = np.where(f < 0.1)[0][0]
a, err = scop.curve_fit(Decay, r[cutoff:], f[cutoff:])
g_piNN = 4*np.pi*E2MeV*E*a[0]/(3.0*e*mpi)

# Print the results
print('B = ', np.round(B, 4))
print('E = ', np.round(E, 4))
print('r_RMS = ', np.round(rRMS, 4))
print('g_piNN = ', np.round(g_piNN, 3))
print('Virial = ', abs(E4/(E2 + 3*(E0-E6)) - 1))

# Plots of energy density
fig = plt.figure(figsize = (12, 8))
fig.subplots_adjust(top = 0.93, bottom = 0.11, left = 0.07, right = 0.96, 
                    hspace = 0, wspace = 0.14)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(r, f, 'r-', label = 'f')
ax1.plot(r, df, 'b-', label = 'df')
ax1.plot(r, ddf, 'g-', label = 'ddf')

ax2.plot(r, 4*np.pi*r**2*Ed, 'k-', label = 'Ed')
ax2.plot(r, 4*np.pi*r**2*Ed2, 'b-', label = 'Ed2')
ax2.plot(r, 4*np.pi*r**2*Ed4, 'g-', label = 'Ed4')
ax2.plot(r, 4*np.pi*r**2*Ed6, 'r-', label = 'Ed6')
ax2.plot(r, 4*np.pi*r**2*Ed0, 'c-', label = 'Ed0')

ax1.legend(fontsize = 13)
ax2.legend(fontsize = 13)

# Isospin quantum contributions
def Isospin(r, f, df, coefs):
    Lambda_2 = 16*np.pi/3*np.sin(f)**2
    Lambda_4 = 64*np.pi*np.sin(f)**2/(3*r**2)*(r**2*df**2 + np.sin(f)**2)
    Lambda_6 = 16*np.pi*np.sin(f)**4*df**2/(3*r**2)
    
    Lambda = (coefs[0]*Lambda_2 + coefs[1]*Lambda_4 + coefs[2]*Lambda_6)/(24*np.pi**2)
    return scin.simpson(r**2*Lambda, x = r)

Lambda = Isospin(r, f, df, coefs)
print('Lambda = ', np.round(Lambda, 5))

# Extract the SO(4) representation of the Skyrme field
f = np.insert(f, 0, np.pi)
r = np.insert(r, 0, 0)
f_interp = scinter.interp1d(r, f)

# Create the cartesian grid
dx = 0.2
xmax = 15.0
x = np.arange(-xmax, xmax+dx, dx)
Nx = len(x)

X, Y, Z = np.meshgrid(x, x, x, indexing = 'ij')

# Calculate the spherical coordinates
RHO = np.sqrt(X**2 + Y**2)
R = np.sqrt(X**2 + Y**2 + Z**2)
THETA = np.arctan2(RHO, Z)
PHI = np.arctan2(Y, X)

# Calculate sigma, pi1, pi2, and pi3
sigma = np.cos(f_interp(R))
pi1 = np.sin(f_interp(R))*np.sin(THETA)*np.cos(PHI)
pi2 = np.sin(f_interp(R))*np.sin(THETA)*np.sin(PHI)
pi3 = np.sin(f_interp(R))*np.cos(THETA)

# Save the data
#np.savetxt('sigma.dat', sigma.flatten())
#np.savetxt('pi1.dat', pi1.flatten())
#np.savetxt('pi2.dat', pi2.flatten())
#np.savetxt('pi3.dat', pi3.flatten())

plt.show()
