# -*- coding: utf-8 -*-
"""
@author: HuidobroMG

Description:
    
    This code is just to plot a given solution of a skyrmion. It takes the
    baryon and energy densities as well as the field solutions to represent
    the skyrmion.

"""

#------------------------------------------------------------------------------

# Packages
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------

# Read the data. Choose the skyrmion you want to plot, n = 1, 2, 3, 4
n = '3'

dataED = np.loadtxt('ED_'+n+'.dat').T
dataBD = np.loadtxt('BD_'+n+'.dat').T
datasigma = np.loadtxt('sigma_'+n+'.dat').T
datapi1 = np.loadtxt('pi1_'+n+'.dat').T
datapi2 = np.loadtxt('pi2_'+n+'.dat').T
datapi3 = np.loadtxt('pi3_'+n+'.dat').T

#------------------------------------------------------------------------------

# Set up the parameters
N = int(np.round(len(dataED)**(1/3), 2))
dx = 0.2
xmax = (N-1)*dx/2

x = np.arange(-xmax, xmax+dx, dx)
X, Y, Z = np.meshgrid(x, x, x)

ED = np.reshape(dataED, (N, N, N))

E = np.sum(ED)*dx**3

# Plot the contours
idx = np.where(ED == np.max(ED))[0][0]
plt.contourf(x, x, ED[:, :, idx], levels = 15, vmin = np.min(ED), vmax = np.max(ED))

#------------------------------------------------------------------------------

from mayavi import mlab

colour = []    
for i in range(len(datapi1)):
    if abs(datapi1[i])>abs(datapi2[i]) and abs(datapi1[i])>abs(datapi3[i]):
        if datapi1[i] > 0:
            c = 0
        else:
            c = 1
    elif abs(datapi2[i])>abs(datapi1[i]) and abs(datapi2[i])>abs(datapi3[i]):
        if datapi2[i] > 0:
            c = 2
        else:
            c = 3
    else:
        if datapi3[i] > 0:
            c = 4
        else:
            c = 5
            
    colour.append(c)

newcol = np.reshape(np.array(colour),(N, N, N))

#------------------------------------------------------------------------------

# 3D plot
x, y, z = np.mgrid[-xmax:xmax:N*1j,-xmax:xmax:N*1j,-xmax:xmax:N*1j]

mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
src = mlab.pipeline.scalar_field(x, y, z, ED)
src.image_data.point_data.add_array(newcol.T.ravel())
src.image_data.point_data.get_array(1).name = 'angle'
src.update()
src2 = mlab.pipeline.set_active_attribute(src, point_scalars='scalar')
contour = mlab.pipeline.contour(src2)
contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='angle')
surf = mlab.pipeline.surface(contour2, transparent = False)

# Our colour map following Leed's colouring scheme
lut = np.array([[255, 255, 0, 255], [0, 255, 255, 255], [255, 0, 0, 255], [0, 255, 0, 255], [255, 255, 255, 255], [0, 0, 0, 255]])
surf.module_manager.scalar_lut_manager.lut.table = lut   
mlab.colorbar(title = r'$\pi_i$', orientation = 'vertical', nb_labels = len(lut))
#mlab.outline(surf, color = (0,0,0), extent = (-xmax,xmax,-xmax,xmax,-xmax,xmax))
