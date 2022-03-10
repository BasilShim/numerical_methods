import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)

# ux(0, t) = u(0, t)
# ux(pi, t) = u(pi, t)
# ux(0, 0) = u(0, 0) = sin(0) + cos(0) = 1
# ux(pi, 0) = u(pi, 0) = sin(pi) + cos(pi) = -1

def ExactSolution(x, t, a):
    result = np.zeros(len(x))
    for i in range(len(x)):
        temp = np.sin(x[i] - a * t[-1]) + np.cos(x[i] + a * t[-1])
        result[i] = temp

    return result

def u(x):
    result = np.sin(x) + np.cos(x)
    
    return result

def ut(x, a):
    result = - a * (np.sin(x) + np.cos(x))

    return result

def FirstOrderApprox(x, t, h):
    result = np.zeros(len(t))
    result[0] = u(x)

    for i in range(1, len(result) - 1):
        result[i] = h * result[i - 1] + result[i - 1]

    return result

def TridiagSolve(a_vec, b_vec, c_vec, d_vec):

    n = len(d_vec)
    roots = np.zeros(len(d_vec))
    p = np.zeros(n)
    q = np.zeros(n)

    p[0] = -c_vec[0] / b_vec[0]
    q[0] = d_vec[0] / b_vec[0]

    for i in range(1, n - 1):
        temp = (b_vec[i] + a_vec[i - 1] * p[i - 1])
        p[i] = - c_vec[i] / temp
        q[i] = (d_vec[i] - a_vec[i - 1] * q[i - 1]) / temp

    q[-1] = (d_vec[-1] - a_vec[-1] * q[-2]) / (b_vec[-1] + a_vec[-1] * p[-2])

    roots[-1] = q[-1]
    for i in range(len(q) - 2, -1, -1):
        roots[i] = p[i] * roots[i + 1] + q[i]

    return roots

def ExFinDiff(x, t, sigma, tau, a, h, order, approx):
    grid = np.zeros((len(t), len(x)))
    grid[0] = u(x)

    if (order == "first"):
        grid[1, 1:-1] = tau * ut(x[1:-1], a) + grid[0][1:-1]
    if (order == "second"):
        for i in range (1, len(grid[1]) - 1):
            grid[1][i] = grid[0][i] + tau * ut(x[i], a) + sigma * (grid[0][i + 1] - 2 * grid[0][i] + grid[0][i - 1]) / 2
    if(approx == "first"): 
        grid[1][0] = grid[1][1] / (1 + h)
        grid[1][-1] = grid[1][-2] / (1 - h)
    if(approx == "second"):
        grid[1][0] = (4 * grid[1][1] - grid[1][2]) / (2 * h + 3)
        grid[1][-1] = (4 * grid[1][-2] - grid[1][-3]) / (3 - 2 * h)

    for i in range(2, len(t)):
        for j in range(1, len(grid[i]) - 1):
            grid[i][j] = sigma * (grid[i - 1][j + 1] - 2 * grid[i - 1][j] + grid[i - 1][j - 1]) - grid[i - 2][j] + 2 * grid[i - 1][j]
        if(approx == "first"):
            grid[i][0] = grid[i][1] / (1 + h)
            grid[i][-1] = grid[i][-2] / (1 - h)
        if(approx == "second"):
            grid[i][0] = (4 * grid[i][1] - grid[i][2]) / (2 * h + 3)
            grid[i][-1] = (4 * grid[i][-2] - grid[i][-3]) / (3 - 2 * h)

    np.savetxt("exp_grid.txt", grid, fmt='%.5f')
    return grid

def ImFinDiff(x, t, sigma, tau, a, h, order, approx):
    grid = np.zeros((len(t), len(x)))
    grid[0] = u(x)

    if (order == "first"):
        grid[1, 1:-1] = tau * ut(x[1:-1], a) + grid[0][1:-1]
    if (order == "second"):
        for i in range (1, len(grid[1]) - 1):
            grid[1][i] = grid[0][i] + tau * ut(x[i], a) + sigma * (grid[0][i + 1] - 2 * grid[0][i] + grid[0][i - 1]) / 2
    if(approx == "first"): 
        grid[1][0] = grid[1][1] / (1 + h)
        grid[1][-1] = grid[1][-2] / (1 - h)
    if(approx == "second"):
        grid[1][0] = (4 * grid[1][1] - grid[1][2]) / (2 * h + 3)
        grid[1][-1] = (4 * grid[1][-2] - grid[1][-3]) / (3 - 2 * h)

    d_vec = - grid[0, 1:-1] + 2 * grid[1, 1:-1]

    a_vec = np.zeros(len(d_vec) - 1)
    c_vec = np.zeros(len(d_vec) - 1)
    b_vec = np.zeros(len(d_vec))
    a_vec.fill(-sigma)
    c_vec.fill(-sigma)
    b_vec.fill(1 + 2 * sigma)

    for i in range(2, len(t)):
        grid[i, 1:-1] = TridiagSolve(a_vec, b_vec, c_vec, d_vec)
        if(approx == "first"):
            grid[i][0] = grid[i][1] / (1 + h)
            grid[i][-1] = grid[i][-2] / (1 - h)
        if(approx == "second"):
            grid[i][0] = (4 * grid[i][1] - grid[i][2]) / (2 * h + 3)
            grid[i][-1] = (4 * grid[i][-2] - grid[i][-3]) / (3 - 2 * h)

    np.savetxt("imp_grid.txt", grid, fmt='%.5f')

    return grid

def main():
    T = 0.01
    N = 25
    K = 25
    a = 1

    l = np.pi
    tau = T / K
    h = l / N

    sigma = (a * a * tau) / (h * h)
    print("Sigma:\n", sigma)

    x = np.arange(0, l + h, h)
    t = np.arange(0, T + tau, tau)

    # print(x)
    # print(t)

    analitical = ExactSolution(x, t, a)
    print("Analitical Solution:\n", analitical)

    sec_explicit_fd = ExFinDiff(x, t, sigma, tau, a, h, order="second", approx="first")[-1]
    print("Explicit Finite Difference:\n", sec_explicit_fd)

    sec_implicit_fd = ImFinDiff(x, t, sigma, tau, a, h, order="second", approx="first")[-1]
    print("Implicit Difference:\n", sec_implicit_fd)

    explicit_fd = ExFinDiff(x, t, sigma, tau, a, h, order="second", approx="second")[-1]
    print("Explicit Finite Difference:\n", explicit_fd)

    implicit_fd = ImFinDiff(x, t, sigma, tau, a, h, order="second", approx="second")[-1]
    print("Implicit Difference:\n", implicit_fd)

    plt.plot(x, analitical, label = "Analitical Solution", color="lime")
    plt.plot(x, explicit_fd, label = "Explicit Solution, First Order", color="red")
    plt.plot(x, implicit_fd, label = "Implicit Solution, First Order", color="blue")
    plt.plot(x, sec_explicit_fd, label = "Explicit Solution, Second Order", color="orange")
    plt.plot(x, sec_implicit_fd, label = "Implicit Solution, Second Order", color="black")
    plt.xlabel("x")
    plt.ylabel("U(x,t)")
    plt.xticks(np.arange(x[0], x[-1], 0.2))
    plt.autoscale(True)
    plt.legend()
    plt.show()

    return


main()