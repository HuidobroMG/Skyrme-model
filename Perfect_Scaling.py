"""
@author: HuidobroMG

We propose a systematic fitting of the parameters in the generalized Skyrme model for the minimal energy
configuration: the FCC half Skyrme crystal, based on an approximation.
This approach allows to fit the energy, baryon density and the symmetry energy at the saturation point,
extremely fast and with sufficient accuracy.
The approximation assumes that around the minimum of energy, each term in the energy functional
scales exactly with the number of derivatives of the Skyrme fields.
"""

# Import the modules
import numpy as np
import scipy.optimize as scop
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

Fit_2params = False
Fit_3params = False
Minimizer = False

# Fitting values
E0_val = 923 # Energy scale of INM, [MeV]. The value is 923 = 939 + aV(= -16 MeV).
n0_val = 0.16 # Nuclear saturation density, [1/fm**3]
S0_val = 31.7 # Symmetry energy of INM, [MeV]
J0_val = 57.7 # First multipole of the symmetry energy curve, [MeV]

# Parameters of the model
hbarc = 197.3269804 # [MeV*fm]
fpi = 129 # [MeV]
e = 5.45 # Skyrme parameter
lambda2 = 2 # Omega-Pion coupling constant, [MeV*fm**3]
mpi = 138 # Pion mass, [MeV]

# Adimensional cpupling constants
hbar_adim = e**2/(3*np.pi**2)
c0 = 2*(mpi/(fpi*e))**2
c6 = 2*lambda2*fpi**2*e**4/hbarc**3

# Perfect Scaling Universal constants
model = '2460'

if model == '24':
    #E_24
    A = 0.11077375
    B = 2.437328
    C = 0
    D = 0
    #Lambda_24
    Ap = 0.03927
    Bp = 1.276225
    Cp = 0
    
elif model == '240':
    # E_240
    A = 0.11363425
    B = 2.408945
    C = 0
    D = 0.00849656
    #Lambda_240
    Ap = 0.0386186
    Bp = 1.33636
    Cp = 0
    
elif model == '246':
    # E_246
    A = 0.1111275
    B = 2.43187
    C = 1.2435328
    D = 0
    #Lambda_246
    Ap = 0.0391755
    Bp = 1.2841125
    Cp = 0.829236

elif model == '2460':
    # E_2460
    A = 0.116453
    B = 2.40419
    C = 1.082434
    D = 0.008449063
    #Lambda_2460
    Ap = 0.038043875
    Bp = 1.3931
    Cp = 0.883148
    
# Energy
def Energy(L):
    return A*L + B/L + c6*C/L**3 + c0*D*L**3

# dE/dL
def dEnergy(L):
    return A - B/L**2 - 3*c6*C/L**4 + 3*c0*D*L**2

# d2E/dL2
def ddEnergy(L):
    return 2*B/L**3 + 12*c6*C/L**5 + 6*c0*D*L

# Isospin moment of inertia
def Isospin(L):
    return Ap*L**3 + Bp*L + c6*Cp/L

# Physical scales
E2MeV = 3*np.pi**2*fpi/e
x2fm = hbarc/(fpi*e)
L0 = scop.fsolve(dEnergy, 1.0)[0]
n0 = 1/(2*(x2fm*L0)**3)
E0 = E2MeV*Energy(L0)
S0 = E2MeV*hbar_adim**2/(2*Isospin(L0))
print('The minimum has energy E = ', np.round(E0/E2MeV, 3))
print('The minimum is in L = ', np.round(L0, 3))
print('The minimum has energy E (MeV)', np.round(E0, 3))
print('The minimum is in n (1/fm**3) = ', np.round(n0, 3))
print('The minimum has symmetry energy S (MeV)', np.round(S0, 3))

# Construct the curves
nmax = 7*n0
nmin = 0.2*n0
Lmin = 1/(2*nmax)**(1/3)/x2fm
Lmax = 1/(2*nmin)**(1/3)/x2fm
L = np.linspace(Lmin, Lmax, 20)
E = Energy(L)
S = E2MeV*hbar_adim**2/(2*Isospin(L))

if Fit_2params == True:
    # Energy
    def Energy(L, params):
        c6 = 2*lambda2*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A*L + B/L + c6*C/L**3 + c0*D*L**3
    
    # dE/dL
    def DerEnergy(L, params):    
        c6 = 2*lambda2*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A - B/L**2 - 3*c6*C/L**4 + 3*c0*D*L**2
    
    # Lambda
    def Lambda(L, params):
        c6 = 2*lambda2*params[0]**2*params[1]**4/hbarc**3
        return Ap*L**3 + Bp*L + c6*Cp/L
    
    # Fit fpi and e
    def Fit_fpi_e(pars):        
        # Change to physical units
        E2MeV = 3*np.pi**2*pars[0]/pars[1]
        x2fm = hbarc/(pars[0]*pars[1])
        hbarc_adim = pars[1]**2/(3*np.pi**2)
        
        # Find the minimum
        Ltilde = scop.fsolve(DerEnergy, 1.0, args=pars)[0]
        print('Lmin = ', Ltilde)
        
        L0 = Ltilde*x2fm
        S0 = E2MeV*hbarc_adim**2/(2*Lambda(Ltilde, pars))
        E0 = Energy(Ltilde, pars)*E2MeV
        n0 = 1/(2*L0**3)
        
        print(np.round(E0, 4), np.round(n0, 4), np.round(S0, 4))
        print('\n')
        
        return np.array([E0-E0_val, n0-n0_val])
    
    roots = scop.fsolve(Fit_fpi_e, [130, 6])
    print('fpi, e = ', roots)

if Fit_3params == True:    
    # Energy
    def Energy(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A*L + B/L + c6*C/L**3 + c0*D*L**3
    
    # dE/dL
    def DerEnergy(L, params):    
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A - B/L**2 - 3*c6*C/L**4 + 3*c0*D*L**2
    
    # Lambda
    def Lambda(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        return Ap*L**3 + Bp*L + c6*Cp/L
    
    def DerLambda(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        return 3*Ap*L**2 + Bp - c6*Cp/L**2

    # Fit fpi, e and lambda2
    def Fit_fpi_e_lambda2(pars):
        # Change to physical units
        E2MeV = 3*np.pi**2*pars[0]/pars[1]
        x2fm = hbarc/(pars[0]*pars[1])
        hbarc_adim = pars[1]**2/(3*np.pi**2)
        
        # Find the minimum
        Ltilde = scop.fsolve(DerEnergy, 1.0, args=pars)[0]
        print('Lmin = ', Ltilde)
        
        L0 = Ltilde*x2fm
        S0 = E2MeV*hbarc_adim**2/(2*Lambda(Ltilde, pars))
        J0 = E2MeV*Ltilde*hbarc_adim**2/(2*Lambda(Ltilde, pars)**2)*DerLambda(Ltilde, pars)
        E0 = Energy(Ltilde, pars)*E2MeV
        n0 = 1/(2*L0**3)
        
        print(np.round(E0, 3), np.round(S0, 3), np.round(n0, 3), np.round(J0, 3))
        print('\n')
        
        return np.array([E0-E0_val, n0-n0_val, S0-S0_val])
    
    roots = scop.fsolve(Fit_fpi_e_lambda2, [130, 6, 3])
    print('fpi, e, lambda2 = ', roots)

if Minimizer == True:    
    # Energy
    def Energy(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A*L + B/L + c6*C/L**3 + c0*D*L**3
    
    # dE/dL
    def DerEnergy(L, params):    
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        c0 = 2*mpi**2/(params[0]*params[1])**2
        return A - B/L**2 - 3*c6*C/L**4 + 3*c0*D*L**2
    
    # Lambda
    def Lambda(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        return Ap*L**3 + Bp*L + c6*Cp/L
    
    def DerLambda(L, params):
        c6 = 2*params[2]*params[0]**2*params[1]**4/hbarc**3
        return 3*Ap*L**2 + Bp - c6*Cp/L**2
    
    # Fit fpi, e and lambda2
    def Min_fpi_e_lambda2(pars):
        print('VALUES of fpi, e, lambda2 = '+str(pars))
        
        # Change to physical units
        E2MeV = 3*np.pi**2*pars[0]/pars[1]
        x2fm = hbarc/(pars[0]*pars[1])
        hbarc_adim = pars[1]**2/(3*np.pi**2)
        
        # Find the minimum
        Ltilde = abs(scop.fsolve(DerEnergy, 1.0, args=pars)[0])
        print('Lmin = ', Ltilde)
        
        L0 = Ltilde*x2fm
        S0 = E2MeV*hbarc_adim**2/(2*Lambda(Ltilde, pars))
        J0 = E2MeV*Ltilde*hbarc_adim**2/(2*Lambda(Ltilde, pars)**2)*DerLambda(Ltilde, pars)
        E0 = Energy(Ltilde, pars)*E2MeV
        n0 = 1/(2*L0**3)
        
        print(np.round(E0, 3), np.round(S0, 3), np.round(n0, 3), np.round(J0, 3))
        print('\n')
        
        return np.sqrt((1-abs(E0/E0_val))**2 + (1-abs(n0/n0_val))**2 + (1-abs(S0/S0_val))**2 + (1-abs(J0/J0_val))**2)
    
    roots = scop.minimize(Min_fpi_e_lambda2, [130, 6, 3], bounds = ((0, np.inf),
                                                                    (0, np.inf),
                                                                    (0, np.inf))).x
    print('fpi, e, lambda2 = ', roots)