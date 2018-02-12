import os

import numpy as np
from pyevtk.hl import gridToVTK


def build_vtk_from_array(grid_info, mesh_z, point_array_names, point_arrays, cell_array_names, cell_arrays,
                         out_path, vtk_out):
    """

    :param grid_info:

    :param mesh_z:

    :param point_array_names:

    :param point_arrays:

    :param cell_array_names:

    :param cell_arrays:

    :param out_path:

    :param vtk_out:
    """

    # Dimensions
    ncol = grid_info[0]
    nrow = grid_info[1]
    delc = grid_info[2]
    delr = grid_info[3]
    x0 = grid_info[4]
    y0 = grid_info[5]

    (nz, ny, nx) = mesh_z.shape

    # Coordinates

    X = np.arange(x0, x0 + (ncol + 1) * delc, delc, dtype='float64')
    Y = np.arange(y0 - (nrow + 1) * delr, y0, delr, dtype='float64')

    x = np.zeros((nx, ny, nz))
    y = np.zeros((nx, ny, nz))
    z = np.zeros((nx, ny, nz))

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x[i, j, k] = X[i]

                y[i, j, k] = Y[j]

                (i_ind, j_ind, k_ind) = (i, j, k)
                if i == nx:
                    i_ind = nx - 1
                if j == ny:
                    j_ind = ny - 1
                if k == nz:
                    k_ind = nz - 1

                z[i, j, k] = mesh_z[k_ind, j_ind, i_ind]  # mesh_pad[i_ind,j_ind,k_ind] #Z[k]
    # End for

    cellData = {}
    for index, cell_array in enumerate(cell_arrays):
        zone = np.zeros((nx - 1, ny - 1, nz - 1))

        for k in range(nz - 1):
            for j in range(ny - 1):
                for i in range(nx - 1):
                    zone[i, j, k] = cell_array[k, j, i]

        cellData[cell_array_names[index]] = zone
    # End for

    pointData = {}
    for index, point_array in enumerate(point_arrays):
        if point_array_names[index] == "z_elev":
            pointData[point_array_names[index]] = z
        else:
            pointData[point_array_names[index]] = point_array
        # End if
    # End for

    gridToVTK(os.path.join(out_path, vtk_out), x, y, z, cellData=cellData, pointData=pointData)


if __name__ == "__main__":
    pass
