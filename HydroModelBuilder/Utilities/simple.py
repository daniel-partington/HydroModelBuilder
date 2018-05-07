

def rect_interpolate(interp_point, known_points):
    """
    :param interp_point:
    :param known_points:
    """

    # Function to linearly interpolate data between 4 points on a rectangular grid
    # where interp point is an (x,y) pair for which we want a z, and all other points are (x,y,z) triplets
    #
    #  point1     point2
    #  .    |       .
    #
    #
    #       + interp_point
    #
    #  .    |       .
    #  point3     point4
    #
    #  Formula for linear interpolation between two points is used:
    a1 = known_points[0][2] + (interp_point[0] - known_points[0][0]) * (known_points[1][2] -
                                                                        known_points[0][2]) / (known_points[1][0] - known_points[0][0])
    a2 = known_points[2][2] + (interp_point[0] - known_points[2][0]) * (known_points[3][2] -
                                                                        known_points[2][2]) / (known_points[3][0] - known_points[2][0])
    z = a1 + (interp_point[1] - known_points[0][1]) * (a2 - a1) / (known_points[2][1] - known_points[0][1])

    return z


def get_four_points(x_loc, y_loc, grid):
    """
    :param x_loc:
    :param y_loc:
    :param grid:
    """

    point = 4 * [None]    # Variables that returns four points bounding x_loc and y_loc

    # Find indices and values of x coordinates bounding x_loc
    if ((x_loc >= grid[0][0]) & (x_loc <= grid[0][-1])):
        for index, x in enumerate(grid[0], start=0):
            if x_loc <= x:
                x_left = x_loc
                x_right = grid[0][index + 1]
                x_index = index
                break
    else:
        x_left = None

    # Find indices and values of y coordinates bounding y_loc
    if ((grid[1][0] >= y_loc) & (y_loc >= grid[1][-1])):
        for index, y in enumerate(grid[1], start=0):
            if y_loc >= y:
                y_up = y_loc
                y_down = grid[1][index - 1]
                y_index = index
                break
    else:
        y_down = None

    if (x_left == None) or (y_down == None):
        return point, None, None

     # Define four points that bound x_loc and y_loc for interpolation
    point[0] = [x_left, y_up]
    point[1] = [x_right, y_up]
    point[2] = [x_left, y_down]
    point[3] = [x_right, y_down]

    return point, x_index, y_index


def get_four_elevations(point, x_index, y_index, z):
    """
    :param point:
    :param x_index:
    :param y_index:
    :param z:
    """

    layer_point = 4 * [None]
    # for index, points in enumerate(point, start=0):
    layer_point[0] = [point[0][0], point[0][1], z[y_index][x_index]]
    layer_point[1] = [point[1][0], point[1][1], z[y_index][x_index + 1]]
    layer_point[2] = [point[2][0], point[2][1], z[y_index - 1][x_index]]
    layer_point[3] = [point[3][0], point[3][1], z[y_index - 1][x_index + 1]]

    return layer_point
