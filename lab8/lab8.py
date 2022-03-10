import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
np.set_printoptions(threshold=sys.maxsize)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)

def ux(y, t, mu1, mu2, a):
    result = np.cos(mu2 * y) * np.exp(-(mu1 ** 2 + mu2 ** 2) * a * t)

    return result

def uy(x, t, mu1, mu2, a):
    result = np.cos(mu1 * x) * np.exp(-(mu1 ** 2 + mu2 ** 2) * a * t)

    return result

def ut(x, y, mu1, mu2):
    result = np.cos(mu1 * x) * np.cos(mu2 * y)

    return result

def ExactSolution(x, y, t, mu1, mu2, a):
    result = np.zeros((len(x), len(y)))

    result = np.cos(mu1 * x) * np.cos(mu2 * y) * np.exp(-(mu1 ** 2 + mu2 ** 2) * a * t)

    return result

def TridiagSolve(a_vec, b_vec, c_vec, d_vec):

    n = len(d_vec)
    roots = np.zeros(len(d_vec))
    p = np.zeros(n)
    q = np.zeros(n)

    p[0] = -c_vec[0] / b_vec[0]
    q[0] = d_vec[0] / b_vec[0]

    for i in range(1, n - 1):
        temp = (b_vec[i] + a_vec[i] * p[i - 1])
        p[i] = - c_vec[i] / temp
        q[i] = (d_vec[i] - a_vec[i] * q[i - 1]) / temp

    q[-1] = (d_vec[-1] - a_vec[-1] * q[-2]) / (b_vec[-1] + a_vec[-1] * p[-2])

    roots[-1] = q[-1]
    for i in range(len(q) - 2, -1, -1):
        roots[i] = p[i] * roots[i + 1] + q[i]

    return roots

def Norm(curr_grid, prev_grid):
    max = 0
    for i in range(curr_grid.shape[0]):
        for j in range(curr_grid.shape[1]):
            if abs(curr_grid[i, j] - prev_grid[i, j]) > max:
                max = abs(curr_grid[i, j] - prev_grid[i, j])

    return max

def AlternatingDir(x, y, t, mu1, mu2, a):
    grid = np.zeros((len(t), len(x), len(y)))
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = t[1] - t[0]

    sigma1 = a * dt / (dx**2)
    sigma2 = a * dt / (dy**2)
    print("sigma1 =", sigma1)
    print("sigma2 =", sigma2)

    for j in range(len(y)):
        for i in range(len(x)):
            grid[0,i,j] = ut(x[i], y[j], mu1, mu2)

    for k in range(1, len(t)):
        for j in range(len(y)):
            grid[k, 0, j] = ux(y[j], t[k], mu1, mu2, a)
        for i in range(len(x)):
            grid[k, i, 0] = uy(x[i], t[k], mu1, mu2, a)

    d1_vec = np.zeros(len(y))
    a1_vec = np.zeros(len(d1_vec))
    c1_vec = np.zeros(len(d1_vec))
    b1_vec = np.zeros(len(d1_vec))
    a1_vec.fill(sigma1 / 2)
    c1_vec.fill(sigma1 / 2)
    b1_vec.fill(-1 - sigma1)
    a1_vec[-1] = 0
    c1_vec[0] = 0
    b1_vec[0] = 1
    b1_vec[-1] = 1

    d2_vec = np.zeros(len(x))
    a2_vec = np.zeros(len(d2_vec))
    c2_vec = np.zeros(len(d2_vec))
    b2_vec = np.zeros(len(d2_vec))
    a2_vec.fill(sigma2 / 2)
    c2_vec.fill(sigma2 / 2)
    b2_vec.fill(-1 - sigma2)
    a2_vec[-1] = 0
    c2_vec[0] = 0
    b2_vec[0] = 1
    b2_vec[-1] = 1

    for k in range(1, len(t)):
        temp_grid = np.zeros((len(x), len(y)))
        for i in range(1, len(x) - 1):
            for j in range(1, len(y) - 1):
                d1_vec[j] = (-sigma2 / 2) * (grid[k - 1, i, j + 1] - 2 * grid[k - 1, i, j] + grid[k - 1, i, j - 1]) - grid[k - 1, i, j] 
            d1_vec[0] = uy(x[i], t[k] - dt / 2, mu1, mu2, a)
            d1_vec[-1] = 0
            temp_grid[i] = TridiagSolve(a1_vec, b1_vec, c1_vec, d1_vec)
        temp_grid[0, :] = ux(y, t[k] - dt / 2, mu1, mu2, a)
        temp_grid[:, 0] = uy(x, t[k] - dt / 2, mu1, mu2, a)
        temp_grid[-1, :] = 0
        temp_grid[:, -1] = 0
        with open("temp_grid.txt", "a") as tempfile:
            np.savetxt(tempfile,temp_grid, fmt='%.5f')
            tempfile.write("\n")
        
        for j in range(1, len(y) - 1):
            for i in range(1, len(x) - 1):
                d2_vec[i] = (-sigma1 / 2) * (temp_grid[i][j + 1] - 2 * temp_grid[i][j] + temp_grid[i][j - 1]) - temp_grid[i][j]
            d2_vec[0] = ux(y[j], t[k], mu1, mu2, a)
            d2_vec[-1] = 0
            grid[k, :, j] = TridiagSolve(a2_vec, b2_vec, c2_vec, d2_vec)
        grid[k, 0, :] = ux(y, t[k], mu1, mu2, a)
        grid[k, :, 0] = uy(x, t[k], mu1, mu2, a)
        grid[k, :, -1] = 0
        grid[k, -1, :] = 0

    with open("sd_grid.txt", "w") as outfile:
        for i in range(len(t)):
            np.savetxt(outfile, grid[i], fmt='%.5f')
            outfile.write("\n")

    return grid.T

def FractionSteps(x, y, t, mu1, mu2, a):
    grid = np.zeros((len(t), len(x), len(y)))
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = t[1] - t[0]

    sigma1 = a * dt / (dx**2)
    sigma2 = a * dt / (dy**2)
    print("sigma1 =", sigma1)
    print("sigma2 =", sigma2)

    for j in range(len(y)):
        for i in range(len(x)):
            grid[0,i,j] = ut(x[i], y[j], mu1, mu2)

    for k in range(1, len(t)):
        for j in range(len(y)):
            grid[k, 0, j] = ux(y[j], t[k], mu1, mu2, a)
        for i in range(len(x)):
            grid[k, i, 0] = uy(x[i], t[k], mu1, mu2, a)

    d1_vec = np.zeros(len(y))
    a1_vec = np.zeros(len(d1_vec))
    c1_vec = np.zeros(len(d1_vec))
    b1_vec = np.zeros(len(d1_vec))
    a1_vec.fill(sigma1)
    c1_vec.fill(sigma1)
    b1_vec.fill(-1 - 2*sigma1)
    a1_vec[-1] = 0
    c1_vec[0] = 0
    b1_vec[0] = 1
    b1_vec[-1] = 1

    d2_vec = np.zeros(len(x))
    a2_vec = np.zeros(len(d2_vec))
    c2_vec = np.zeros(len(d2_vec))
    b2_vec = np.zeros(len(d2_vec))
    a2_vec.fill(sigma2)
    c2_vec.fill(sigma2)
    b2_vec.fill(-1 - 2 * sigma2)
    a2_vec[-1] = 0
    c2_vec[0] = 0
    b2_vec[0] = 1
    b2_vec[-1] = 1

    for k in range(1, len(t)):
        temp_grid = np.zeros((len(x), len(y)))
        for i in range(1, len(x) - 1):
            for j in range(1, len(y) - 1):
                d1_vec[j] = - grid[k - 1, i, j] 
            d1_vec[0] = uy(x[i], t[k] - dt / 2, mu1, mu2, a)
            d1_vec[-1] = 0
            temp_grid[i] = TridiagSolve(a1_vec, b1_vec, c1_vec, d1_vec)
        temp_grid[0, :] = ux(y, t[k] - dt / 2, mu1, mu2, a)
        temp_grid[:, 0] = uy(x, t[k] - dt / 2, mu1, mu2, a)
        temp_grid[-1, :] = 0
        temp_grid[:, -1] = 0
        with open("temp_grid.txt", "a") as tempfile:
            np.savetxt(tempfile,temp_grid, fmt='%.5f')
            tempfile.write("\n")
        
        for j in range(1, len(y) - 1):
            for i in range(1, len(x) - 1):
                d2_vec[i] = - temp_grid[i][j]
            d2_vec[0] = ux(y[j], t[k], mu1, mu2, a)
            d2_vec[-1] = 0
            grid[k, :, j] = TridiagSolve(a2_vec, b2_vec, c2_vec, d2_vec)
        grid[k, 0, :] = ux(y, t[k], mu1, mu2, a)
        grid[k, :, 0] = uy(x, t[k], mu1, mu2, a)
        grid[k, :, -1] = 0
        grid[k, -1, :] = 0

    with open("fs_grid.txt", "w") as outfile:
        for i in range(len(t)):
            np.savetxt(outfile, grid[i], fmt='%.5f')
            outfile.write("\n")

    return grid.T

def main():
    # mu1 = 1
    # mu2 = 1
    a = 1

    mu1 = 1
    mu2 = 2

    # mu1 = 2
    # mu2 = 1

    dx = np.pi / 32
    dy = np.pi / 32
    dt = 0.01
    
    ly = mu2 * np.pi / 2
    lx = mu1 * np.pi / 2
    lt = 0.1

    x = np.arange(0, lx + dx, dx)
    y = np.arange(0, ly + dy, dy)
    t = np.arange(0, lt + dt, dt)

    print(len(x))
    print(len(y))
    print(len(t))
    print(x[-1])
    print(np.cos(x[-1]))

    eps = 0.0001
    print("eps =", eps)

    c = 0.5

    file = open("temp_grid.txt", "w")
    file.close()
    alternating_dir = AlternatingDir(x, y, t, mu1, mu2, a)
    
    fractionsteps = FractionSteps(x, y, t, mu1, mu2, a)
    X, Y = np.meshgrid(x, y)
    analitical = ExactSolution(X, Y, t[-1], mu1, mu2, a)
    np.savetxt("analyt.txt", analitical, fmt='%.5f')
    error_x = []
    error_y = []
    error_t = []
    for i in range(len(x)):
        #error_y.append(max(abs(analitical[:, i] - alternating_dir[:, i, -1])))
        error_y.append(max(abs(analitical[:, i] - fractionsteps[:, i, -1])))
    for j in range(len(y)):
        #error_x.append(max(abs(analitical[j, :] - alternating_dir[j, :, -1])))
        error_x.append(max(abs(analitical[j, :] - fractionsteps[j, :, -1])))
    for k in range(len(t)):
        #error_t.append(Norm(ExactSolution(X, Y, t[k], mu1, mu2, a), alternating_dir[:, :, k]))
        error_t.append(Norm(ExactSolution(X, Y, t[k], mu1, mu2, a), fractionsteps[:, :, k]))
    plt.title("График ошибок")
    plt.plot(x, error_y, label = "При фиксированном x в выбранный момент времени", color = "red")
    plt.plot(y, error_x, label = "При фиксированном y в выбранный момент времени", color = "blue")
    plt.plot(t, error_t, label = "По 'x' и 'y' во всех временных промежутках", color = "green")
    plt.xlabel("x, y, t")
    plt.ylabel("error")
    plt.grid()
    plt.legend()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, analitical, label="Exact", color="red")
    # ax.plot_surface(X, Y, alternating_dir[:,:,-1], label="Alternating Directions", color="green")
    ax.plot_surface(X, Y, fractionsteps[:,:,-1], label="Fractional Steps")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("U")
    plt.show()

    return

main()