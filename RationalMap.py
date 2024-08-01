"""
@author: HuidobroMG

We construct and compute the baryon and energy densities of a desired Rational Map configuration.
Additionally, a contour plot is performed to visualize the energy configuration at the height of maximal value.
"""

# Import the modules
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Nfields = 4
Ndims = 3

dx = 0.2
N = 101
xmax = (N-1)*dx/2
x = np.arange(-xmax, xmax+dx, dx)

c0 = 0 # c0 = 2*mpi**2/(fpi*e)**2
c6 = 0 # c6 = 2*lambda**2*fpi**2*e**4

# Levi-Civita symbol
def LeviCivita(a,b,c,d):
    return (a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d)/12.

# Derivatives
def D(phi):
    # dphi[a][i] = d_i n_a
    dphi = np.zeros((Nfields, Ndims, N-4, N-4, N-4))
    for a in range(Nfields):
        dphi[a][0] = (phi[a][:-4, 2:-2, 2:-2]-8*phi[a][1:-3, 2:-2, 2:-2]+8*phi[a][3:-1, 2:-2, 2:-2]-phi[a][4:, 2:-2, 2:-2])/(12*dx)
        dphi[a][1] = (phi[a][2:-2, :-4, 2:-2]-8*phi[a][2:-2, 1:-3, 2:-2]+8*phi[a][2:-2, 3:-1, 2:-2]-phi[a][2:-2, 4:, 2:-2])/(12*dx)
        dphi[a][2] = (phi[a][2:-2, 2:-2, :-4]-8*phi[a][2:-2, 2:-2, 1:-3]+8*phi[a][2:-2, 2:-2, 3:-1]-phi[a][2:-2, 2:-2, 4:])/(12*dx)
    
    return dphi

# Densities
def BaryonDensity(phi):
    dphi = D(phi) # Size N - 4
    B0 = 0
    for a in range(Nfields):
        for b in range(Nfields):
            for c in range(Nfields):
                for d in range(Nfields):
                    B0 += LeviCivita(a,b,c,d)*phi[a][2:-2, 2:-2, 2:-2]*dphi[b][0]*dphi[c][1]*dphi[d][2]
    
    return B0/(2*np.pi**2)

def EnergyDensity(phi):
    dphi = D(phi) # Size N - 4
    E2 = np.zeros((N-4, N-4, N-4))
    E4 = np.zeros((N-4, N-4, N-4))
    E6 = (2*np.pi**2*BaryonDensity(phi))**2
    for a in range(Nfields):
        aux = dphi[a][0]**2 + dphi[a][1]**2 + dphi[a][2]**2
        E2 += aux
        for b in range(Nfields):
            E4 -= 2*(dphi[a][0]*dphi[b][0] + dphi[a][1]*dphi[b][1] + dphi[a][2]*dphi[b][2])**2
    E4 += 2*E2**2
    E0 = 1-phi[0][2:-2, 2:-2, 2:-2]
    
    return (E2 + E4 + c6*E6 + c0*E0)/(24*np.pi**2)

# Rational Map initial configuration
def RationalMap(x): 
    X, Y, Z = np.meshgrid(x, x, x)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(rho, Z)
    pangle = np.arctan2(Y, X)
    if c0 == 0: # Instanton configuration
        f = np.pi*(1 - (1 + 2.11/(10*r/xmax)**2)**(-1/2))
    else: # Exponential decay
        f = np.pi*np.exp(-10*r/xmax)
    
    Re_z = np.tan(theta/2)*np.cos(pangle)
    Im_z = np.tan(theta/2)*np.sin(pangle)
    
    Denom = 1.0
    R_Re = Re_z**2 - Im_z**2
    R_Im = 2*Re_z*Im_z
    R_Re /= Denom
    R_Im /= Denom
    R_mod2 = R_Re**2 + R_Im**2
    
    s = np.cos(f)
    p1 = np.sin(f)*2*R_Re/(1 + R_mod2)
    p2 = np.sin(f)*2*R_Im/(1 + R_mod2)
    p3 = np.sin(f)*(1 - R_mod2)/(1 + R_mod2)
    return np.array([s, p1, p2, p3])

# Test the code
phi0 = RationalMap(x)

# Compute the densities
B0 = BaryonDensity(phi0)
ED = EnergyDensity(phi0)

# Integrate the densities
B = np.sum(B0)*dx**3
E = np.sum(ED)*dx**3
print('B = ', np.round(B, 6), 'E = ', np.round(E, 6))

# Plot the densities
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.contourf(x[2:-2], x[2:-2], B0[:, :, N//2], levels = 20)
ax2.contourf(x[2:-2], x[2:-2], ED[:, :, N//2], levels = 20)

plt.show()