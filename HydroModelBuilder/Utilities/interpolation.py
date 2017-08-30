"""
Created on Wed May 11 21:21:51 2016

@author: part0075
"""

import numpy as np
from scipy.interpolate import Rbf
from scipy.interpolate import griddata


def GridData(points, values, xi, method='nearest'):
    """
    Options for griddata method are ['nearest', 'linear', 'cubic']

    This has been modded to use nearest method outside of interpolation points
    convex hull.
    """
    if method in ['linear', 'cubic']:
        z1 = griddata(points, values, xi, method=method)
        z0 = griddata(points, values, xi, method='nearest')
        z2 = z1
        z2[np.isnan(z1)] = z0[np.isnan(z1)]
    else:
        z2 = griddata(points, values, xi, method=method)
    # End if

    return z2
# End GridData()


def RadialBasisFunctions2D(points, values, xi, function='multiquadric', epsilon=2):
    """
    Available types of interpolation are:
    ['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate']
    """
    rbf = Rbf(np.array(points)[:, 0], np.array(points)[:, 1], np.array(values), epsilon=epsilon, function=function)
    ZI = rbf(xi[0], xi[1])
    return ZI
# End RadialBasisFunctions2D()


def RadialBasisFunctions3D(points, values, xi, function='multiquadric', epsilon=2):
    """
    Available types of interpolation are:
    ['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate']
    """
    rbf = Rbf(points[0], points[1], points[2], values, epsilon=epsilon, function=function)
    ValuesI = rbf(xi[0], xi[1], xi[2])
    return ValuesI
# End RadialBasisFunctions3D()


def Interpolator(mesh_type, points, values, to_array, method='linear', use='griddata', function='multiquadric', epsilon=2):
    if mesh_type == 'structured' and len(points.shape) == 2:
        if use == 'griddata':
            to_values = GridData(points, values, to_array, method=method)
        elif use == 'rbf':
            to_values = RadialBasisFunctions2D(points, values, to_array, function=function, epsilon=epsilon)
        else:
            print 'Interpolation "use" not recognised'
        # End if
    elif mesh_type == 'structured' and len(points.shape) > 2:
        pass
    elif mesh_type == 'unstructured':
        pass
    else:
        print 'mesh type not recognised: ', mesh_type
    # End if

    return to_values
# End Interpolator()


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from matplotlib import cm

    # Testing
    def func(x, y):
        return x * (1 - x) * np.cos(4 * np.pi * x) * np.sin(4 * np.pi * y**2)**2

    grid_x, grid_y = np.mgrid[0:1:50j, 0:1:200j]

    xi = (grid_x, grid_y)

    points = np.random.rand(1000, 2)
    values = func(points[:, 0], points[:, 1])

    values_interp = Interpolator('structured', points, values, xi, method='cubic')

    def func(x, y):
        return x * (1 - x) * np.cos(4 * np.pi * x) * np.sin(4 * np.pi * y**2)**2

    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

    xi = (grid_x, grid_y)

    points = np.random.rand(1000, 2)
    values = func(points[:, 0], points[:, 1])

    grid_z0 = Interpolator('structured', points, values, xi, method='nearest')
    grid_z1 = Interpolator('structured', points, values, xi, method='linear')
    grid_z2 = Interpolator('structured', points, values, xi, method='cubic')

    plt.subplot(221)
    plt.imshow(func(grid_x, grid_y).T, extent=(0, 1, 0, 1), origin='lower')
    plt.plot(points[:, 0], points[:, 1], 'k.', ms=1)
    plt.title('Original')
    plt.subplot(222)
    plt.imshow(grid_z0.T, extent=(0, 1, 0, 1), origin='lower')
    plt.title('Nearest')
    plt.subplot(223)
    plt.imshow(grid_z1.T, extent=(0, 1, 0, 1), origin='lower')
    plt.title('Linear')
    plt.subplot(224)
    plt.imshow(grid_z2.T, extent=(0, 1, 0, 1), origin='lower')
    plt.title('Cubic')
    plt.gcf().set_size_inches(6, 6)
    plt.show()

    # Test for radial basis functions

    x = np.random.rand(100) * 4.0 - 2.0
    y = np.random.rand(100) * 4.0 - 2.0
    z = x * np.exp(-x**2 - y**2)

    x = np.array(points)[:, 0]
    y = np.array(points)[:, 1]

    # 2-d tests - setup scattered data
    ti = np.linspace(-2.0, 2.0, 100)
    XI, YI = np.meshgrid(ti, ti)
    epsilon = 0.01

    ZI = Interpolator('structured', points, values, xi, use='rbf', function='multiquadric', epsilon=epsilon)

    plt.figure()
    plt.subplot(1, 1, 1)
    #plt.pcolor(XI, YI, ZI, cmap=cm.jet)
    plt.pcolor(grid_x, grid_y, ZI, cmap=cm.jet)
    #plt.scatter(x, y, 100, z, cmap=cm.jet)
    plt.scatter(x, y, 40, values, cmap=cm.jet)
    plt.title('RBF interpolation - multiquadrics')
    #plt.xlim(-2, 2)
    #plt.ylim(-2, 2)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.colorbar()
