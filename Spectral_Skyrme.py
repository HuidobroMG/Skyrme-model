"""
@author: HuidobroMG

We find the B = 1 skyrmion profile function using the spectral methods approach.
This problem is extremely useful since it involves a two-domain and non-linear resolution,
so it represents a sufficiently general an interesting application of spectral methods.
As a comment, the system is solved in Skyrme units:
    E2MeV = 3*np.pi**2*fpi/e
    x2fm = hbarc/(fpi*e) 
"""

# Import the modules
import numpy as np
import Polynomials as Pol
import matplotlib.pyplot as plt
from numpy.linalg import inv

# Boundary between domains
r_cut = 4.0

# First domain, r1 in [0, r_cut], x1 in [-1, 1]
# r1 = r_cut/2*(1 + x1)
N1 = 100
Ncoefs1 = 20

# Grid 1
x1 = -np.cos(np.pi*np.arange(0, N1+1, 1)/N1) # (CGL Nodes)
w1 = np.pi*np.ones(N1+1)/N1 # (CGL Weights)
w1[0] /= 2
w1[-1] /= 2
N1 = len(x1)

# Chebyshev polynomials and the coefficients of the expansions
T1 = Pol.Ti(Ncoefs1, x1)
dT1 = Pol.dTi(Ncoefs1)
xT1 = Pol.xTi(Ncoefs1)
T_x1 = Pol.Ti_x(Ncoefs1)

# Construct the linear operators
Ld1 = Pol.Derivative(Ncoefs1, dT1)
Lx1 = Pol.Product(Ncoefs1, xT1)
L_x1 = Pol.Cocient(Ncoefs1, T_x1)
Ldd1 = np.dot(Ld1, Ld1)
Lx1_I = np.eye(Ncoefs1, Ncoefs1) + Lx1

# Linear operator of the ODE
L1 = Lx1_I.dot(Lx1_I).dot(Ldd1) + 2*Lx1_I.dot(Ld1)

# Initial configuration
r1 = r_cut/2*(1 + x1)
dr_dx = r_cut/2

f1 = np.pi/(1+r1**2) # Profile function initialization

# Find the coefficients of the interpolation and construct the derivatives
ci_1 = Pol.interpolate(f1, Ncoefs1, T1, w1, N1)
dci_1 = np.dot(ci_1, dT1)
ddci_1 = np.dot(dci_1, dT1)
f1 = np.dot(ci_1, T1)
df1 = np.dot(dci_1, T1) # This is df/dx, df/dr = df/dx * dx/dr
ddf1 = np.dot(ddci_1, T1)

# Define the source of the ODE
def source1(ci):
    # Construct the field and its derivatives
    dci = np.dot(ci, dT1)
    ddci = np.dot(dci, dT1)
    f = np.dot(ci, T1)
    df = np.dot(dci, T1)
    ddf = np.dot(ddci, T1)
    # Change from df/dx to df/dr
    df *= 1/dr_dx
    ddf *= (1/dr_dx)**2
    # Now write the source term
    s = np.sin(2*f) - 8*np.sin(f)**2*ddf - 4*np.sin(2*f)*df**2 + 4*np.sin(f)**2*np.sin(2*f)/r1**2
    s[0] = 0 # Singularity in r = 0
    return s

s1 = source1(ci_1)
si_1 = Pol.interpolate(s1, Ncoefs1, T1, w1, N1)


# Second domain, r2 in [r_cut, inf], u in [1/r_cut, 0], x2 in [-1, 1]
# u = 1/r2 = (1 - x2)/(2*r_cut)
N2 = 100
Ncoefs2 = 20

# Grid 2
x2 = -np.cos(np.pi*np.arange(0, N2+1, 1)/N2) # (CGL Nodes)
w2 = np.pi*np.ones(N2+1)/N2 # (CGL Weights)
w2[0] /= 2
w2[-1] /= 2
N2 = len(x2)

# Chebyshev polynomials
T2 = Pol.Ti(Ncoefs2, x2)
dT2 = Pol.dTi(Ncoefs2)
xT2 = Pol.xTi(Ncoefs2)
T_x2 = Pol.Ti_x(Ncoefs2)

# Construct the linear operators
Ld2 = Pol.Derivative(Ncoefs2, dT2)
Lx2 = Pol.Product(Ncoefs2, xT2)
L_x2 = Pol.Cocient(Ncoefs2, T_x2)
Ldd2 = np.dot(Ld2, Ld2)
I_Lx2 = np.eye(Ncoefs2, Ncoefs2) - Lx2

# Linear operator in the second domain
L2 = I_Lx2.dot(I_Lx2).dot(Ldd2)

# Initial configuration
r2 = 2*r_cut/(1 - x2)
u2 = 1/r2
du_dx = -1/(2*r_cut)

f2 = np.pi/(1 + r2**2) # Profile function initialization

# Find the coefficients of the interpolation and construct the derivatives
ci_2 = Pol.interpolate(f2, Ncoefs2, T2, w2, N2)
dci_2 = np.dot(ci_2, dT2)
ddci_2 = np.dot(dci_2, dT2)
f2 = np.dot(ci_2, T2)
df2 = np.dot(dci_2, T2)
ddf2 = np.dot(ddci_2, T2)

# Define the source of the ODE
def source2(ci):
    # Construct the fields and its derivatives
    dci = np.dot(ci, dT2)
    ddci = np.dot(dci, dT2)
    f = np.dot(ci, T2)
    df = np.dot(dci, T2)
    ddf = np.dot(ddci, T2)
    # Change from df/dx to df/du
    df *= 1/du_dx
    ddf *= (1/du_dx)**2
    # Now write the source term
    s = np.sin(2*f) - 8*u2**3*np.sin(f)**2*(2*df + u2*ddf) - 4*np.sin(2*f)*u2**4*df**2 + 4*np.sin(2*f)*u2**2*np.sin(f)**2
    return s

s2 = source2(ci_2)
si_2 = Pol.interpolate(s2, Ncoefs2, T2, w2, N2)


# Impose the boundary conditions
# f(r = 0) = np.pi, f(r = inf) = 0
A = np.pi
B = 0

# The first B.C. corresponds to the first domain
# f1(x = -1) = np.pi
for i in range(Ncoefs1):
    L1[-2,i] = (-1)**i
si_1[-2] = A

# The second B.C. corresponds to the second domain
# f2(x = 1) = 0
L2[-2,:] = 1
si_2[-2] = B

# Now concatenate both domains
Nt = Ncoefs1 + Ncoefs2
L = np.zeros((Nt, Nt))
L[:Ncoefs1, :Ncoefs1] = L1
L[Ncoefs1:, Ncoefs1:] = L2

si = np.concatenate((si_1, si_2))
ci = np.concatenate((ci_1, ci_2))

# Impose the continuity conditions
for i in range(Ncoefs1):
    # f1(x = 1) = f2(x = -1)
    L[Ncoefs1-1, i] = -1
    
    # df1/dx(x = 1) = df2/dx(x = -1)
    L[-1, i] = sum(Ld1[:, i])

for i in range(Ncoefs2):
    # f1(x = 1) = f2(x = -1)
    L[Ncoefs1-1, i+Ncoefs1] = (-1)**i
    
    # df1/dx(x = 1) = df2/dx(x = -1)
    L[-1, i+Ncoefs1] = (-1)**i*sum(Ld2[:, i])

si[Ncoefs1-1] = 0
si[-1] = 0


# Residuals function
def Res(ci):
    ci_1 = 1.0*ci[:Ncoefs1]
    s1 = source1(ci_1)
    si_1 = Pol.interpolate(s1, Ncoefs1, T1, w1, N1)
    
    ci_2 = 1.0*ci[Ncoefs1:]
    s2 = source2(ci_2)
    si_2 = Pol.interpolate(s2, Ncoefs2, T2, w2, N2)
    
    si_1[-2] = A
    si_2[-2] = B
    
    si = np.concatenate((si_1, si_2))
    si[Ncoefs1-1] = 0
    si[-1] = 0
    return np.dot(L, ci) - si

# Obtain the Jacobian matrix numerically, J_ij = dRes_i/dc_j
def Jij(ci):
    J = np.zeros((Nt, Nt))
    
    ci_new = 1.0*ci
    for j in range(Nt):
        ci_new[j] = 1.001*ci[j]
        R1 = Res(ci_new)
        ci_new[j] = 0.999*ci[j]
        R2 = Res(ci_new)
        ci_new[j] = 1.0*ci[j]
        eps = 0.001*ci[j]
        J[:, j] = (R1 - R2)/(2*eps)
    
    J[Ncoefs1-2, :] = 1.0*L[Ncoefs1-2, :]
    J[Ncoefs1-1, :] = 1.0*L[Ncoefs1-1, :]
    J[-2, :] = 1.0*L[-2, :]
    J[-1, :] = 1.0*L[-1, :]
            
    return J


# Start the iteration
iterations = 500
for i in range(iterations):   
    R = Res(ci)
    J = Jij(ci)
    
    err = np.sqrt(np.sum(R**2))
    print('|Res| = ', err)
    if err < 1e-7:
        break
    Jinv = inv(J)    
    X = np.dot(Jinv, R)
    ci -= 0.1*X
    
    ci_1 = 1.0*ci[:Ncoefs1]
    f1 = np.dot(ci_1, T1)
    ci_2 = 1.0*ci[Ncoefs1:]
    f2 = np.dot(ci_2, T2)
    
    s1 = source1(ci_1)
    si_1 = Pol.interpolate(s1, Ncoefs1, T1, w1, N1)
    s2 = source2(ci_2)
    si_2 = Pol.interpolate(s2, Ncoefs2, T2, w2, N2)
    
    # Refresh the boundary and continuity conditions
    si_1[-2] = A
    si_2[-2] = B
    
    si = np.concatenate((si_1, si_2))
    si[Ncoefs1-1] = 0
    si[-1] = 0