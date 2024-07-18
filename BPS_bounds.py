"""
@author: HuidobroMG

We compute the different topological BPS energy bounds from the 4 possible modifications of the Skyrme model.
"""

# Import the modules
import numpy as np
import scipy.optimize as scop

# Parameters
hbarc = 197.3269804 # [MeV*fm]
fpi = 117.23 # Pion decay constant [MeV]
e = 4.81 # Skyrme parameter
lambda2 = 0.0 # Omega-Pion coupling constant [MeV*fm**3]
mpi = 138 # Pion mass [MeV]

# Energy scale
E2MeV = 3*np.pi**2*fpi/e

# Redefined constants
mu2 = (2*np.pi**2)**(2/3)*fpi**2/(8*hbarc)
mu4 = (2*np.pi**2)**(4/3)*hbarc/(2*e**2)
mu6 = lambda2*np.pi**4
mu0 = (mpi*fpi)**2/(8*hbarc**3)

# Second power of the integral of the square root of the pion-mass potential term on the target space
I2 = 2.1333

# Fourth power of the integral of the quartic root of the pion-mass potential term on the target space
I4 = 1.80739

# L24 Bound
L24 = 1 # By definition of the energy and length scale

# L240 Bound
a = np.sqrt(mu2**2/(mu0*mu4*(I4*2/np.pi)**4))
Bound_240 = np.sqrt(mu4)*(4*np.sqrt(mu2/a)*(1 - a**2/2*(np.sqrt(1+4/a**2)-1))**(3/4) + 3*np.sqrt(2*mu2)*a*(np.sqrt(1+4/a**2)-1)**0.5)

Bound_240 /= E2MeV
print('Bound_240 = ', np.round(Bound_240, 3))

# L246 Bound        
c = np.sqrt(mu4**2/(mu2*mu6))
Bound_246 = np.sqrt(mu2)*(3*np.sqrt(2*mu4)*c*np.sqrt(np.sqrt(1+4/c**2)-1)+4*np.sqrt(mu4/c)*(1-c**2/2*(np.sqrt(1+4/c**2)-1))**(3/4))

Bound_246 /= E2MeV
print('Bound_246 = ', np.round(Bound_246, 3))

# L2460 Bound
def fun(pars):
    alpha0, alpha2, alpha4, alpha6 = pars
    
    Bound = 2*(I2*2/np.pi)*(mu0*alpha0*mu6*alpha6)**0.5 + 4*(I4*2/np.pi)*(mu4*alpha4)**(3/4)*(mu0*(1-alpha0))**(1/4) + 6*(mu4*(1-alpha4)*mu2*alpha2)**0.5 + 4*(mu6*(1-alpha6))**(1/4)*(mu2*(1-alpha2))**(3/4)
    return 1/Bound

inits = np.ones(4)*0.5

results = scop.minimize(fun, inits, 
                        bounds = ((0,1), (0,1), (0,1), (0,1)),
                        method = 'Nelder-Mead')

Bound_2460 = 1/fun(results.x)

Bound_2460 /= E2MeV
print('Bound_2460 = ', np.round(Bound_2460, 3))