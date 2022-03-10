import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
np.set_printoptions(threshold=sys.maxsize)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)

# ux(0, y) = 0
# u(1, y) = 1 - y^2
# uy(x, 0) = 0
# u(x, 1) = x^2 - 1

def ExactSolution(x, y):
    result = x ** 2 - y ** 2

    return result

def uy(y):
    result = 1 - y * y

    return result

def ux(x):
    result = x * x - 1

    return result

def Norm(curr_grid, prev_grid):
    max = 0
    for i in range(curr_grid.shape[0]):
        for j in range(curr_grid.shape[1]):
            if abs(curr_grid[i, j] - prev_grid[i, j]) > max:
                max = abs(curr_grid[i, j] - prev_grid[i, j])

    return max

def Liebman(x, y, eps):
    curr_grid = np.zeros((len(y), len(x)))
    curr_grid[-1, :] = ux(x)
    curr_grid[:, -1] = uy(y)

    prev_grid = np.zeros((len(y), len(x)))
    count = 0
    h = x[1] - x[0]

    while Norm(curr_grid, prev_grid) > eps:
        count += 1
        prev_grid = np.copy(curr_grid)
        curr_grid[0, :] = curr_grid[1, :]
        curr_grid[:, 0] = curr_grid[:, 1]
        for i in range(1, len(y) - 1):
            for j in range(1, len(x) - 1):
                curr_grid[i][j] = (prev_grid[i + 1][j] + prev_grid[i - 1][j] + prev_grid[i][j - 1] + prev_grid[i][j + 1]) / 4

    np.savetxt("liebman.txt", curr_grid, fmt='%.5f')

    return curr_grid

def Relaxation(x, y, eps, c):
    curr_grid = np.zeros((len(y), len(x)))
    curr_grid[-1, :] = ux(x)
    curr_grid[:, -1] = uy(y)

    prev_grid = np.zeros((len(y), len(x)))
    count = 0
    h = x[1] - x[0]

    while Norm(curr_grid, prev_grid) > eps:
        count += 1
        prev_grid = np.copy(curr_grid)
        curr_grid[0, :] = curr_grid[1, :]
        curr_grid[:, 0] = curr_grid[:, 1]
        for i in range (1, len(y) - 1):
            for j in range(1, len(x) - 1):
                curr_grid[i][j] = (1 - c) * prev_grid[i][j] + c * (prev_grid[i + 1][j] + prev_grid[i - 1][j] + prev_grid[i][j - 1] + prev_grid[i][j + 1]) / 4

    np.savetxt("relaxation.txt", curr_grid, fmt='%.5f')

    return curr_grid

def Zeidel(x, y, eps, c):
    curr_grid = np.zeros((len(y), len(x)))
    curr_grid[-1, :] = ux(x)
    curr_grid[:, -1] = uy(y)

    prev_grid = np.zeros((len(y), len(x)))
    count = 0
    h = x[1] - x[0]
    
    while Norm(curr_grid, prev_grid) > eps:
        count += 1
        prev_grid = np.copy(curr_grid)
        curr_grid[0, :] = curr_grid[1, :]
        curr_grid[:, 0] = curr_grid[:, 1]
        for i in range(1, len(y) - 1):
            for j in range(1, len(x) - 1):
                curr_grid[i][j] = (1 - c) * prev_grid[i][j] + c * (prev_grid[i + 1][j] + curr_grid[i - 1][j] + curr_grid[i][j - 1] + prev_grid[i][j + 1]) / 4

    np.savetxt("zeidel.txt", curr_grid, fmt='%.5f')

    return curr_grid

def main():

    dx = 0.05
    dy = 0.05
    
    ly = 1
    lx = 1

    x = np.arange(0, lx + dx, dx)
    y = np.arange(0, ly + dy, dy)

    print(x)
    print(y)

    eps = 0.000001
    print(eps)

    c = 0.5

    liebman = Liebman(x, y, eps)

    relaxation = Relaxation(x, y, eps, c)

    zeidel = Zeidel(x, y, eps, c)

    X, Y = np.meshgrid(x, y)
    analitical = ExactSolution(X, Y)
    error_x_l = []
    error_y_l = []
    error_x_r = []
    error_y_r = []
    error_x_z = []
    error_y_z = []
    for i in range(len(y)):
        # error_x_l.append(max(abs(analitical[:, i] - liebman[:, i])))
        # error_y_l.append(max(abs(analitical[i, :] - liebman[i, :])))
        # error_x_r.append(max(abs(analitical[:, i] - relaxation[:, i])))
        # error_y_r.append(max(abs(analitical[i, :] - relaxation[i, :])))
        error_x_z.append(max(abs(analitical[:, i] - zeidel[:, i])))
        error_y_z.append(max(abs(analitical[i, :] - zeidel[i, :])))
    # plt.plot(x, error_y_l, label = "Fixed x Liebman Error")
    # plt.plot(x, error_x_l, label = "Fixed y Liebman Error")
    # plt.plot(x, error_y_r, label = "Fixed x Relaxation Error")
    # plt.plot(x, error_x_r, label = "Fixed y Relaxation Error")
    plt.plot(x, error_y_z, label = "Fixed x Zeidel Error")
    plt.plot(x, error_x_z, label = "Fixed y Zeidel Error")
    plt.legend()
    figer = plt.figure()
    axx = figer.add_subplot(111, projection="3d")
    # axx.plot_surface(X, Y, abs(analitical - liebman))
    # axx.plot_surface(X, Y, abs(analitical - relaxation))
    # axx.plot_surface(X, Y, abs(analitical - zeidel))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, analitical, label="Exact")
    ax.plot_surface(X, Y, liebman, label="Liebman")
    ax.plot_surface(X, Y, relaxation, label="Relaxation")
    ax.plot_surface(X, Y, zeidel, label="Zeidel")
    plt.legend()
    plt.show()

    return


main()