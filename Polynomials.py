# -*- coding: utf-8 -*-
"""
@author: HuidobroMG

Description:
    
    In this code we define everything we may need to solve a problem using
    spectral methods. Concretely, we define some operations on the Chebyshev
    and Legendre polynomials.
        
"""

#------------------------------------------------------------------------------

# Packages
import numpy as np
import scipy.special as scsp

#-----------------------------------------------------------------------------

# Polynomials stuff

# Chebyshev-Gauss-Lobatto quadrature integration
def CGL(f, w, N):
    return sum(f[i]*w[i] for i in range(N))

# Integration of Chebyshev polynomials <Ti|Tj>
def TiTj(i, j):
    if i != j:
        return 0
    else:
        if i == 0:
            return np.pi
        else:
            return np.pi/2

# Generator of Chebyshec polynomials
def Ti(Ncoefs, x):
    N = len(x)
    
    P = np.zeros((Ncoefs, N))
    P[0] = 1
    P[1] = x
    for i in range(2, Ncoefs):
        P[i] = 2*x*P[i-1]-P[i-2] # Chebyshev polynomials
        
    return P

# Derivative matrix of Chebyshev polynomials
def dTi(Ncoefs):
    dP = np.zeros((Ncoefs, Ncoefs))
    dP[1,0] = 1
    dP[2,1] = 4
    for i in range(3, Ncoefs):
        dP[i][i-1] = 2*i
        for j in range(0, i-2):
            dP[i][j] = i/(i-2)*dP[i-2][j]

    return dP

# Matrix of the product x * P(x)
def xTi(Ncoefs):
    xP = np.zeros((Ncoefs, Ncoefs))
    xP[0,1] = 1
    xP[1,2] = xP[1,0] = 0.5
    for i in range(2, Ncoefs-1):
        xP[i][i+1] = 0.5
        xP[i][i-1] = 0.5

    xP[-1][-2] = 0.5
    return xP

# Matrix of the product 1/x * P(x)
def Ti_x(Ncoefs):
    P_x = np.zeros((Ncoefs, Ncoefs))
    P_x[1,0] = 1
    for i in range(2, Ncoefs):
        P_x[i][i-1] = 2
        for j in range(0, i-2):
            P_x[i][j] = -P_x[i-2][j]

    return P_x

# Legendre polynomials generator
def Legendre(Ncoefs, x):
    N = len(x)
    
    Lp = np.zeros((Ncoefs, N))
    for i in range(Ncoefs):
        Lp[i] = scsp.eval_legendre(i, x)
    
    return Lp

# Change of basis Chebyshev to Legendre
def Lp_T(Ncoefs, Lp, T, w, N):
    M = np.zeros((Ncoefs, Ncoefs))
    for i in range(Ncoefs):
        for j in range(Ncoefs):
            M[i,j] = CGL(Lp[i]*T[j], w, N)/CGL(T[j]*T[j], w, N)
    
    return M

# Interpolation of a function
def interpolate(f, Ncoefs, T, w, N):
    # Coefficients of the interpolation
    fi = np.zeros(Ncoefs)
    
    for i in range(Ncoefs):
        fi[i] = CGL(f*T[i], w, N)/CGL(T[i]*T[i], w, N)
    
    return fi


# Matrices operations
def Derivative(Ncoefs, dT):
    Ld = np.zeros((Ncoefs, Ncoefs)) # Derivative
    for i in range(Ncoefs):
        for j in range(Ncoefs):
            val = 0
            for k in range(Ncoefs):
                val += (2-np.eye(Ncoefs)[i,0])*dT[j,k]*TiTj(i,k)
                
            Ld[i,j] = val/np.pi
    
    return Ld

def Product(Ncoefs, xT):
    Lx = np.zeros((Ncoefs, Ncoefs)) # Derivative
    for i in range(Ncoefs):
        for j in range(Ncoefs):
            val = 0
            for k in range(Ncoefs):
                val += (2-np.eye(Ncoefs)[i,0])*xT[j,k]*TiTj(i,k)
                
            Lx[i,j] = val/np.pi
    
    return Lx


def Cocient(Ncoefs, T_x):
    L_x = np.zeros((Ncoefs, Ncoefs)) # Derivative
    for i in range(Ncoefs):
        for j in range(Ncoefs):
            val = 0
            for k in range(Ncoefs):
                val += (2-np.eye(Ncoefs)[i,0])*T_x[j,k]*TiTj(i,k)
                
            L_x[i,j] = val/np.pi
    
    return L_x