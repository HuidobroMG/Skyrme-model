"""
@author: HuidobroMG

We solve the field equations of the Skyrme model for a generic value of the topological charge using a
gradient descent minimization method. For the initial configuration we use the rational map approximation,
which fix the symmetry of the final solution.
The system is solved using the Skyrme units:
    E2MeV = 3*np.pi**2*fpi/e
    x2fm = hbarc/(fpi*e)
and the O(4)-representation field solutions are saved in .dat files.
"""

# Import the modules
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

# Grid parameters
Nfields = 4
Ndims = 3

dx = 0.2
N = 51
xmax = (N-1)*dx/2
x = np.arange(-xmax, xmax+dx, dx)

# Adimensional coupling constants
c0 = 0 # c0 = 2*mpi**2/(fpi*e)**2
c6 = 0 # c6 = 2*lambda**2*fpi**2*e**4

# Levi-Civita symbol
def LeviCivita(a,b,c,d):
    return (a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d)/12.

# Finite differences first derivative of the fields
def D(phi):
    # dphi[a][i] = d_i n_a
    dphi = np.zeros((Nfields, Ndims, N-4, N-4, N-4))
    for a in range(Nfields):
        dphi[a][0] = (phi[a][:-4, 2:-2, 2:-2]-8*phi[a][1:-3, 2:-2, 2:-2]+8*phi[a][3:-1, 2:-2, 2:-2]-phi[a][4:, 2:-2, 2:-2])/(12*dx)
        dphi[a][1] = (phi[a][2:-2, :-4, 2:-2]-8*phi[a][2:-2, 1:-3, 2:-2]+8*phi[a][2:-2, 3:-1, 2:-2]-phi[a][2:-2, 4:, 2:-2])/(12*dx)
        dphi[a][2] = (phi[a][2:-2, 2:-2, :-4]-8*phi[a][2:-2, 2:-2, 1:-3]+8*phi[a][2:-2, 2:-2, 3:-1]-phi[a][2:-2, 2:-2, 4:])/(12*dx)
    return dphi

# Finite differences second derivative of the fields
def DD(phi):
    # ddphi[a][i][j] = d_i d_j n_a
    ddphi = np.zeros((Nfields, Ndims, Ndims, N-4, N-4, N-4))
    for a in range(Nfields):
        ddphi[a][0][0] = (-phi[a][4:, 2:-2, 2:-2] + 16.0*phi[a][3:-1, 2:-2, 2:-2] - 30.0*phi[a][2:-2, 2:-2, 2:-2] + 16.0*phi[a][1:-3, 2:-2, 2:-2] - phi[a][:-4, 2:-2, 2:-2])/(12.0*dx*dx)
        ddphi[a][0][1] = (phi[a][4:, 2:-2, 2:-2] + phi[a][:-4, 2:-2, 2:-2] + phi[a][2:-2, :-4, 2:-2] + phi[a][2:-2, 4:, 2:-2] - phi[a][4:, 4:, 2:-2] - phi[a][:-4, :-4, 2:-2] - 16.0*phi[a][3:-1, 2:-2, 2:-2] + 30.0*phi[a][2:-2, 2:-2, 2:-2] - 16.0*phi[a][1:-3, 2:-2, 2:-2] - 16.0*phi[a][2:-2, 1:-3, 2:-2] - 16.0*phi[a][2:-2, 3:-1, 2:-2] + 16.0*phi[a][3:-1, 3:-1, 2:-2] + 16.0*phi[a][1:-3, 1:-3, 2:-2])/(24.0*dx*dx)
        ddphi[a][0][2] = (phi[a][4:, 2:-2, 2:-2] + phi[a][:-4, 2:-2, 2:-2] + phi[a][2:-2, 2:-2, :-4] + phi[a][2:-2, 2:-2, 4:] - phi[a][4:, 2:-2, 4:] - phi[a][:-4, 2:-2, :-4] - 16.0*phi[a][3:-1, 2:-2, 2:-2] + 30.0*phi[a][2:-2, 2:-2, 2:-2] - 16.0*phi[a][1:-3, 2:-2, 2:-2] - 16.0*phi[a][2:-2, 2:-2, 1:-3] - 16.0*phi[a][2:-2, 2:-2, 3:-1] + 16.0*phi[a][3:-1, 2:-2, 3:-1] + 16.0*phi[a][1:-3, 2:-2, 1:-3])/(24.0*dx*dx)
        ddphi[a][1][1] = (-phi[a][2:-2, 4:, 2:-2] + 16.0*phi[a][2:-2, 3:-1, 2:-2] - 30.0*phi[a][2:-2, 2:-2, 2:-2] + 16.0*phi[a][2:-2, 1:-3, 2:-2] - phi[a][2:-2, :-4, 2:-2])/(12.0*dx*dx)
        ddphi[a][1][2] = (phi[a][2:-2, 4:, 2:-2] + phi[a][2:-2, :-4, 2:-2] + phi[a][2:-2, 2:-2, :-4] + phi[a][2:-2, 2:-2, 4:] - phi[a][2:-2, 4:, 4:] - phi[a][2:-2, :-4, :-4] - 16.0*phi[a][2:-2, 3:-1, 2:-2] + 30.0*phi[a][2:-2, 2:-2, 2:-2] - 16.0*phi[a][2:-2, 1:-3, 2:-2] - 16.0*phi[a][2:-2, 2:-2, 1:-3] - 16.0*phi[a][2:-2, 2:-2, 3:-1] + 16.0*phi[a][2:-2, 3:-1, 3:-1] + 16.0*phi[a][2:-2, 1:-3, 1:-3])/(24.0*dx*dx)
        ddphi[a][2][2] = (-phi[a][2:-2, 2:-2, 4:] + 16.0*phi[a][2:-2, 2:-2, 3:-1] - 30.0*phi[a][2:-2, 2:-2, 2:-2] + 16.0*phi[a][2:-2, 2:-2, 1:-3] - phi[a][2:-2, 2:-2, :-4])/(12.0*dx*dx)
        ddphi[a][1][0] = 1.0*ddphi[a][0][1]
        ddphi[a][2][0] = 1.0*ddphi[a][0][2]
        ddphi[a][2][1] = 1.0*ddphi[a][1][2]
    return ddphi

# Baryon density
def BaryonDensity(phi):
    dphi = D(phi) # Size N - 4
    
    B0 = 0
    for a in range(Nfields):
        for b in range(Nfields):
            for c in range(Nfields):
                for d in range(Nfields):
                    B0 += LeviCivita(a,b,c,d)*phi[a][2:-2, 2:-2, 2:-2]*dphi[b][0]*dphi[c][1]*dphi[d][2]
    return B0/(2*np.pi**2)

# Energy density
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

# Field equations
def EulerLagrange(phi):
    dphi = D(phi) # Size N - 4
    ddphi = DD(phi) # Size N - 4
    
    EL = np.zeros((Nfields, N-4, N-4, N-4))
    EL_2 = np.zeros((N-4, N-4, N-4))
    EL_4 = np.zeros((N-4, N-4, N-4))
    EL_6 = np.zeros((N-4, N-4, N-4))
    for a in range(Nfields):
        EL_4[:] = 0
        EL_2 = 2*(ddphi[a][0][0] + ddphi[a][1][1] + ddphi[a][2][2])
        for b in range(Nfields):
            # There is a factor 2 from EL_2
            EL_4 += 4*EL_2*(dphi[b][0]**2 + dphi[b][1]**2 + dphi[b][2]**2)
            for i in range(Ndims):
                EL_4 += 16*dphi[a][i]*(ddphi[b][i][0]*dphi[b][0] + ddphi[b][i][1]*dphi[b][1] + ddphi[b][i][2]*dphi[b][2])
                EL_4 -= 8*dphi[b][i]*(dphi[b][0]*ddphi[a][0][i] + dphi[b][1]*ddphi[a][1][i] + dphi[b][2]*ddphi[a][2][i])
                EL_4 -= 8*ddphi[b][i][i]*(dphi[a][0]*dphi[b][0] + dphi[a][1]*dphi[b][1] + dphi[a][2]*dphi[b][2])
                EL_4 -= 8*dphi[b][i]*(ddphi[b][i][0]*dphi[a][0] + ddphi[b][i][1]*dphi[a][1] + ddphi[b][i][2]*dphi[a][2])
    
        EL[a] = EL_2 + EL_4 + c6*EL_6 + c0*(a == 0)
    return EL

# Initial configurations
def Hedgehog(x):   
    X, Y, Z = np.meshgrid(x, x, x)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(rho, Z)
    pangle = np.arctan2(Y, X)
    f = np.pi/(1.0 + r**2)
    
    # B = 1 hedgehog ansatz
    s = np.cos(f)
    p1 = np.sin(f)*np.sin(theta)*np.cos(pangle)
    p2 = np.sin(f)*np.sin(theta)*np.sin(pangle)
    p3 = np.sin(f)*np.cos(theta)
    
    return np.array([s, p1, p2, p3])

def RationalMap(x):   
    X, Y, Z = np.meshgrid(x, x, x)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(rho, Z)
    pangle = np.arctan2(Y, X)
    f = np.pi/(1.0 + r**2)
    
    # This is the B = 1 case
    Re_z = np.tan(theta/2)*np.cos(pangle)
    Im_z = np.tan(theta/2)*np.sin(pangle)
    
    R_Re = Re_z
    R_Im = Im_z
    R_mod2 = R_Re**2 + R_Im**2
    
    s = np.cos(f)
    p1 = np.sin(f)*2*R_Re/(1 + R_mod2)
    p2 = np.sin(f)*2*R_Im/(1 + R_mod2)
    p3 = np.sin(f)*(1 - R_mod2)/(1 + R_mod2)
    return np.array([s, p1, p2, p3])

def Init(x):
    fields = Hedgehog(x)
    #fields = RationalMap(x)
    return fields

# Construct the field
phi0 = Init(x)
B0 = BaryonDensity(phi0)
ED = EnergyDensity(phi0)

B = np.sum(B0)*dx**3
E = np.sum(ED)*dx**3
print('B = ', np.round(B, 6), 'E = ', np.round(E, 6))

# Gradient flow minimization
iterations = 1000
delta = 2e-5
delta_phi = np.zeros((Nfields, N, N, N))
for i in range(iterations):
    print('iteration = ', i)
    delta_phi[:, 2:-2, 2:-2, 2:-2] = 1.0*EulerLagrange(phi0)
    phi_new = phi0 + delta*delta_phi
    
    # Renormalize the field
    norma = np.sqrt(phi_new[0]**2 + phi_new[1]**2 + phi_new[2]**2 + phi_new[3]**2)
    phi_new /= norma
    
    B0 = BaryonDensity(phi_new)
    ED = EnergyDensity(phi_new)
    
    B = np.sum(B0)*dx**3
    E = np.sum(ED)*dx**3
    print('B = ', np.round(B, 6), 'E = ', np.round(E, 6))
    print('-----------------------------------------------')
    phi0 = 1.0*phi_new

# Save the field solution
#np.savetxt('sigma.dat', phi0[0].flatten())
#np.savetxt('pi1.dat', phi0[1].flatten())
#np.savetxt('pi2.dat', phi0[2].flatten())
#np.savetxt('pi3.dat', phi0[3].flatten())