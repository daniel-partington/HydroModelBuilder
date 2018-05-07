# -*- coding: utf-8 -*-
import numpy as np
from pyevtk.hl import gridToVTK 

# Dimensions 

nz, ny, nx = 2, 3, 4

ncells = nx * ny * nz 

npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Coordinates 

#x1 = np.arange(x0, ncol*delc, delc) 

#y1 = np.arange(y0-nrow*delr, nrow*delr, delr)  

#z1 = np.arange(0, (layers+1)*1, 1)


# Coordinates 

delc = 1.0
delr = 1.0
x0 = 0.0
y0 = 4.0

X = np.arange(x0, x0+(nx+1)*delc, delc, dtype='float64') 

Y = np.arange(y0-(ny+1)*delr, y0, delr, dtype='float64') 

Z = np.arange(0, (nz+1)*1, 1, dtype='float64') 

x = np.zeros((nx + 1, ny + 1, nz + 1)) 

y = np.zeros((nx + 1, ny + 1, nz + 1)) 

z = np.zeros((nx + 1, ny + 1, nz + 1))


for k in range(nz + 1): 
    for j in range(ny + 1):
        for i in range(nx + 1): 
 
            x[i,j,k] = X[i] 

            y[i,j,k] = Y[j]
            
            (i_ind, j_ind, k_ind) = (i,j,k)
            if i==nx:
                i_ind = nx-1
            if j==ny:
                j_ind = ny-1
            if k==nz:
                z_ind = nz-1
                
            z[i,j,k] = Z[k] # mesh_pad[i_ind,j_ind,k_ind] #Z[k]
                

gridToVTK("./mesh_test", x, y, z, pointData = {"z_elev" : z}) # cellData = {"zone" : zone}, 



